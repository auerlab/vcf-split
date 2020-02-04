/***************************************************************************
 *  Description:
 *      Generate single-sample VCF files suitable for haplohseq from
 *      dbGaP multisample BCFs.
 *
 *      Input is a VCF stream, usually piped in from bcftools view,
 *      containing multiple samples per call.  vcf-split outputs a
 *      single-sample VCF file for each sample in the columms specified
 *      on the command-line.  Output filenames are based on the sample
 *      ID for each column.
 *
 *  Arguments:
 *      first-column
 *      last-column
 *
 *  Returns:
 *      See "man sysexits"
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Wayne BaconBegin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>     // PATH_MAX
#include <stdint.h>     // SIZE_MAX
#include <ctype.h>
#include <stdbool.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include "vcf-split.h"
#include "tsvio.h"
#include "vcfio.h"

int     main(int argc,const char *argv[])

{
    char        *eos;
    const char  *outfile_prefix,
		*selected_samples_file;
    id_list_t   *selected_sample_ids;
    size_t      first_col,
		last_col,
		max_calls = SIZE_MAX,
		next_arg = 1;
    flag_t      flags = 0;
    
    if ( (argc < 6) || (argc > 10) )
	usage(argv);
    
    next_arg = 1;
    while ( argv[next_arg][0] == '-' )
    {
	if ( strcmp(argv[next_arg],"--max-calls") == 0 )
	{
	    max_calls = strtoul(argv[++next_arg], &eos, 10);
	    if ( *eos != '\0' )
	    {
		fprintf(stderr, "%s: %s: Max calls must be a positive integer.\n",
			argv[0], argv[next_arg]);
		exit(EX_DATAERR);
	    }
	    ++next_arg;
	}
	
	if ( strcmp(argv[next_arg],"--sample-id-file") == 0 )
	{
	    selected_samples_file = argv[++next_arg];
	    ++next_arg;
	    selected_sample_ids =
		read_selected_sample_ids(argv, selected_samples_file);
	    fprintf(stderr, "%zu selected samples.\n", selected_sample_ids->count);
	}
	
	/*
	 *  Adding filter flags like --genotype het to bcftools can slow
	 *  it down significantly.  Choosing flags that maximize the speed
	 *  of bcftools and filtering for heterozygous sites here will make
	 *  the analysis run much faster.
	 */
	
	if ( strcmp(argv[next_arg], "--het-only") == 0 )
	{
	    flags |= FLAG_HET;
	    ++next_arg;
	}

	/* Find a way to do this without thousands of xz processes
	if ( strcmp(argv[next_arg], "--xz") == 0 )
	{
	    flags |= FLAG_XZ;
	    ++next_arg;
	}
	*/
    }
    
    outfile_prefix = argv[next_arg++];
    
    first_col = strtoul(argv[next_arg], &eos, 10);
    if ( *eos != '\0' )
    {
	fprintf(stderr, "%s: %s: Column must be a positive integer.\n",
		argv[0], argv[next_arg]);
	exit(EX_DATAERR);
    }
    
    last_col = strtoul(argv[++next_arg], &eos, 10);
    if ( *eos != '\0' )
    {
	fprintf(stderr, "%s: %s: Column must be a positive integer.\n",
		argv[0], argv[next_arg]);
	exit(EX_DATAERR);
    }
    
    if ( (first_col < 1) || (last_col < 1) || (first_col > last_col) )
    {
	fprintf(stderr, "%s: Columns must be >= 1 and last_col >= first_col.\n", argv[0]);
	exit(EX_DATAERR);
    }
    
    if ( last_col - first_col >= MAX_OUTFILES )
    {
	fprintf(stderr, "%s: Maximum columns is %u\n", argv[0], MAX_OUTFILES);
	exit(EX_USAGE);
    }

    return vcf_split(argv, stdin, outfile_prefix, first_col, last_col,
		     selected_sample_ids, max_calls, flags);
}


/***************************************************************************
 *  Description:
 *      Split a multisample VCF stream into single-sample files.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

int     vcf_split(const char *argv[], FILE *vcf_infile,
		  const char *outfile_prefix,
		  size_t first_col, size_t last_col,
		  id_list_t *selected_sample_ids, size_t max_calls,
		  flag_t flags)

{
    char    inbuf[BUFF_SIZE + 1],
	    *sample_ids[last_col - first_col + 1];
    bool    selected[last_col - first_col + 1];
    
    // Input is likely to come from "bcftools view" stdout.
    // What is optimal buffering for a Unix pipe?  Benchmark several values.
    setvbuf(vcf_infile, inbuf, _IOFBF, BUFF_SIZE);
    vcf_skip_header(argv, vcf_infile);
    vcf_get_sample_ids(argv, vcf_infile, sample_ids, first_col, last_col);
    tag_selected_columns(sample_ids, selected_sample_ids, selected,
			 first_col, last_col);
    write_output_files(argv, vcf_infile, (const char **)sample_ids, selected,
		       outfile_prefix, first_col, last_col, max_calls, flags);
    
    return EX_OK;
}


/***************************************************************************
 *  Description:
 *      Write individual VCF files for each sample.
 *
 *  Arguments:
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

void    write_output_files(const char *argv[], FILE *vcf_infile,
			    const char *sample_ids[], bool selected[],
			    const char *outfile_prefix,
			    size_t first_col, size_t last_col,
			    size_t max_calls, flag_t flags)

{
    size_t  columns = last_col - first_col + 1,
	    c;
    int     fd;
    FILE    *vcf_outfiles[columns];
    char    filename[PATH_MAX + 1],
	    cmd[CMD_MAX + 1];
    
    // Open all output streams
    for (c = 0; c < columns; ++c)
    {
	if ( selected[c] )
	{
	    snprintf(filename, PATH_MAX, "%s%s.vcf",
		     outfile_prefix, sample_ids[c]);
	    /* Find a way to do this without thousands of xz processes
	    if ( flags & FLAG_XZ )
	    {
		snprintf(cmd, CMD_MAX, "xz --stdout > %s.xz",filename);
		vcf_outfiles[c] = popen(cmd, "w");
	    }
	    else
	    */
	    vcf_outfiles[c] = fopen(filename, "w");
	    if ( vcf_outfiles[c] == NULL )
	    {
		fprintf(stderr, "%s: Cannot create %s: %s.\n",
			argv[0], filename, strerror(errno));
		exit(EX_CANTCREAT);
	    }
	    tsv_skip_rest_of_line(argv, vcf_infile);
	}
    }

    // Temporary hack for testing.  Remove limit.
    c = 0;
    while ( (c < max_calls) && split_line(argv, vcf_infile, vcf_outfiles,
	    sample_ids, selected, first_col, last_col, flags) )
	++c;
    
    // Close all output streams
    for (c = 0; c < columns; ++c)
    {
	if ( selected[c] )
	{
	    /*
	    if ( flags & FLAG_XZ )
		pclose(vcf_outfiles[c]);
	    else
	    */
	    fclose(vcf_outfiles[c]);
	    snprintf(filename, PATH_MAX, "%s%s.vcf.done",
		     outfile_prefix, sample_ids[c]);
	    
	    /*
	     *  Touch a .done file to indicate completion.  Another script
	     *  can use this to determine which .vcf files are ready for
	     *  compression.
	     */
	    if ( (fd = open(filename, O_CREAT | O_TRUNC | O_WRONLY)) != -1 )
		close(fd);
	    else
		fprintf(stderr, "%s: Warning: Could not create %s: %s.\n",
			argv[0], filename, strerror(errno));
	}
    }
}


/***************************************************************************
 *  Description:
 *      Split a single input line from the multisample VCF.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

int     split_line(const char *argv[], FILE *vcf_infile, FILE *vcf_outfiles[],
		   const char *sample_ids[], bool selected[],
		   size_t first_col, size_t last_col, flag_t flags)

{
    size_t      c;
    vcf_call_t  vcf_call;
    char        genotype[VCF_SAMPLE_MAX_CHARS + 1];
    
    /*
     *  Read in VCF fields
     */
    
    if ( vcf_read_static_fields(argv, vcf_infile, &vcf_call) )
    {
	//fprintf(stderr, "POS = %s\n", vcf_call.pos_str);
	// Skip columns before first_col
	for (c = 1; c < first_col; ++c)
	    tsv_skip_field(argv, vcf_infile);
	
	for (; c <= last_col; ++c)
	{
	    tsv_read_field(argv, vcf_infile, genotype, VCF_SAMPLE_MAX_CHARS);
	    if ( selected[c - first_col] )
	    {
		if ( !(flags & FLAG_HET) || (genotype[0] != genotype[2]) )
		{
#if DEBUG
		    fprintf(stderr, "%zu %s:\n", c, sample_ids[c - first_col]);
		    fprintf(stderr, "%s\t%s\t.\t%s\t%s\t.\t.\t.\t%s\t%s\n",
			vcf_call.chromosome, vcf_call.pos_str,
			vcf_call.ref, vcf_call.alt,
			vcf_call.format, genotype);
#endif
		    // FIXME: Use vcf_write_ss_call()
		    fprintf(vcf_outfiles[c - first_col],
			"%s\t%s\t.\t%s\t%s\t.\t.\t.\t%s\t%s\n",
			vcf_call.chromosome, vcf_call.pos_str,
			vcf_call.ref, vcf_call.alt, 
			vcf_call.format, genotype);
		}
	    }
	}
	
	tsv_skip_rest_of_line(argv, vcf_infile);
	return 1;
    }
    else
	return 0;
}


/***************************************************************************
 *  Description:
 *      Read list of sample IDs from a file.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-20  Jason Bacon Begin
 ***************************************************************************/

id_list_t   *read_selected_sample_ids(const char *argv[], const char *samples_file)

{
    id_list_t   *list;
    size_t      c;
    char        temp_id[SAMPLE_ID_MAX + 1];
    FILE        *fp;
    
    if ( (fp = fopen(samples_file, "r")) == NULL )
    {
	fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], samples_file,
		strerror(errno));
	exit(EX_NOINPUT);
    }
    
    if ( (list = (id_list_t *)malloc(sizeof(id_list_t))) == NULL )
    {
	fprintf(stderr, "%s: Cannot allocate sample_ids structure.\n", argv[0]);
	exit(EX_UNAVAILABLE);
    }
    
    /*
     *  With such a small file, the cleanest approach is read to count the
     *  IDS, do one malloc(), rewind and read again.
     */

    for (list->count = 0; read_string(fp, temp_id, SAMPLE_ID_MAX) > 0; ++list->count )
	;
    
    if ( (list->ids = (char **)malloc(list->count * sizeof(char **))) == NULL )
    {
	fprintf(stderr, "%s: Cannot allocate sample_ids list.\n", argv[0]);
	exit(EX_UNAVAILABLE);
    }
    
    rewind(fp);
    for (c = 0; c < list->count; ++c)
    {
	read_string(fp, temp_id, SAMPLE_ID_MAX);
	if ( (list->ids[c] = strdup(temp_id)) == NULL )
	{
	    fprintf(stderr, "%s: Cannot allocate sample id %zu.\n", argv[0], c);
	    exit(EX_UNAVAILABLE);
	}
    }
    fclose(fp);
    
    // List may already be sorted, so stay away from qsort, which has
    // worst-case performance on presorted data
    // List should not be huge, so minor issue here
    // heapsort() is not portable, however
    qsort(list->ids, list->count, sizeof(char *),
	(int (*)(const void *, const void *))strptrcmp);
    return list;
}


/***************************************************************************
 *  Description:
 *      Read a whitespace-separated string.
 *
 *  Returns:
 *      String length.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-23  Jason Bacon Begin
 ***************************************************************************/

size_t  read_string(FILE *fp, char *buff, size_t maxlen)

{
    size_t  c = 0;
    int     ch;
    
    while ( isspace(ch = getc(fp)) )
	;
    
    while ( (c < maxlen) && ! isspace(ch) && (ch != EOF) )
    {
	buff[c++] = ch;
	ch = getc(fp);
    }
    buff[c] = '\0';
    return c;
}


void    usage(const char *argv[])

{
    // Find a way to do this without thousands of xz processes
    // fprintf(stderr, "Usage: %s: [--xz] [--het-only] [--max-calls N] [--sample-id-file file] output-file-prefix first-column last-column\n", argv[0]);
    fprintf(stderr, "Usage: %s: [--het-only] [--max-calls N] [--sample-id-file file] output-file-prefix first-column last-column\n", argv[0]);
    exit(EX_USAGE);
}


void    tag_selected_columns(char *sample_ids[], id_list_t *selected_sample_ids,
			     bool selected[], size_t first_col, size_t last_col)

{
    size_t  c;
    
    for (c = first_col; c <= last_col; ++c)
    {
	selected[c - first_col] =
	    (bsearch(sample_ids + c - first_col, selected_sample_ids->ids,
		selected_sample_ids->count, sizeof(char *),
		(int (*)(const void *, const void *))strptrcmp)
		!= NULL);
    }
}


int     strptrcmp(const char **p1, const char **p2)

{
    return strcmp(*p1, *p2);
}
