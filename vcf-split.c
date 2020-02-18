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
#include "tsvio.h"
#include "vcfio.h"
#include "vcf-split.h"

int     main(int argc,const char *argv[])

{
    char        *eos;
    const char  *outfile_prefix,
		*selected_samples_file = NULL;
    id_list_t   *selected_sample_ids = NULL;
    size_t      first_col,
		last_col,
		max_calls = SIZE_MAX,
		next_arg = 1;
    flag_t      flags = 0;
    
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
	
	else if ( strcmp(argv[next_arg],"--sample-id-file") == 0 )
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
	
	else if ( strcmp(argv[next_arg], "--het-only") == 0 )
	{
	    flags |= FLAG_HET;
	    ++next_arg;
	}

	/* Find a way to do this without thousands of xz processes
	else if ( strcmp(argv[next_arg], "--xz") == 0 )
	{
	    flags |= FLAG_XZ;
	    ++next_arg;
	}
	*/
	
	else
	    usage(argv);
    }
    
    outfile_prefix = argv[next_arg++];
    
    first_col = strtoul(argv[next_arg], &eos, 10);
    if ( *eos != '\0' )
    {
	fprintf(stderr, "%s: %s: Column must be a positive integer.\n",
		argv[0], argv[next_arg]);
	usage(argv);
    }
    
    last_col = strtoul(argv[++next_arg], &eos, 10);
    if ( *eos != '\0' )
    {
	fprintf(stderr, "%s: %s: Column must be a positive integer.\n",
		argv[0], argv[next_arg]);
	usage(argv);
    }
    
    if ( (first_col < 1) || (last_col < 1) || (first_col > last_col) )
    {
	fprintf(stderr, "%s: Columns must be >= 1 and last_col >= first_col.\n", argv[0]);
	usage(argv);
    }
    
    if ( last_col - first_col >= MAX_OUTFILES )
    {
	fprintf(stderr, "%s: Maximum columns is %u\n", argv[0], MAX_OUTFILES);
	usage(argv);
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
	    *all_sample_ids[last_col - first_col + 1];
    bool    selected[last_col - first_col + 1];
    size_t  c;
    
    // Input is likely to come from "bcftools view" stdout.
    // What is optimal buffering for a Unix pipe?  Benchmark several values.
    setvbuf(vcf_infile, inbuf, _IOFBF, BUFF_SIZE);
    vcf_skip_header(argv, vcf_infile);
    vcf_get_sample_ids(argv, vcf_infile, all_sample_ids, first_col, last_col);

    /*
    fputs("All sample IDs:", stderr);
    for (c = 0; c < last_col - first_col + 1; ++c)
	fprintf(stderr, " %s", all_sample_ids[c]);
    putc('\n', stderr);
    */
    
    tag_selected_columns(all_sample_ids, selected_sample_ids, selected,
			 first_col, last_col);
    
    fputs("Tagged sample IDs:", stderr);
    for (c = 0; c < selected_sample_ids->count; ++c)
	fprintf(stderr, " %s", selected_sample_ids->ids[c]);
    fputc('\n', stderr);
    
    write_output_files(argv, vcf_infile, (const char **)all_sample_ids, selected,
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
			    const char *all_sample_ids[], bool selected[],
			    const char *outfile_prefix,
			    size_t first_col, size_t last_col,
			    size_t max_calls, flag_t flags)

{
    size_t  columns = last_col - first_col + 1,
	    c;
    int     fd;
    FILE    *vcf_outfiles[columns];
    char    filename[PATH_MAX + 1];
    
    // Open all output streams
    for (c = 0; c < columns; ++c)
    {
	if ( selected[c] )
	{
	    snprintf(filename, PATH_MAX, "%s%s.vcf",
		     outfile_prefix, all_sample_ids[c]);
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
	}
    }

    // Heart of the program, split each VCF line across multiple files
    for (c = 0; (c < max_calls) &&
		split_line(argv, vcf_infile, vcf_outfiles, all_sample_ids,
			   selected, first_col, last_col, flags); ++c)
	;
    
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
	    /*
	     *  Touch a .done file to indicate completion.  Another script
	     *  can use this to determine which .vcf files are ready for
	     *  compression.
	     */
	    snprintf(filename, PATH_MAX, "%s%s.vcf.done",
		     outfile_prefix, all_sample_ids[c]);
	    if ( (fd = open(filename, O_CREAT | O_TRUNC | O_WRONLY, 0644)) != -1 )
		close(fd);
	    else
		fprintf(stderr, "%s: Warning: Could not create %s: %s.\n",
			argv[0], filename, strerror(errno));
	}
    }
    fprintf(stderr, "%s completed successfully.\n", argv[0]);
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
		   const char *all_sample_ids[], bool selected[],
		   size_t first_col, size_t last_col, flag_t flags)

{
    static size_t   line_count = 0;
    size_t      c,
		field_len;
    vcf_call_t  vcf_call;
    char        genotype[VCF_SAMPLE_MAX_CHARS + 1];
    
    /*
     *  Read in VCF fields
     */
    
    if ( vcf_read_static_fields(argv, vcf_infile, &vcf_call) )
    {
	if ( ++line_count % 100 == 0 )
	    fprintf(stderr, "%zu\r", line_count);
	
	//fprintf(stderr, "POS = %s\n", vcf_call.pos_str);
	// Skip columns before first_col
	for (c = 1; c < first_col; ++c)
	    if ( tsv_skip_field(argv, vcf_infile) == EOF )
	    {
		dump_line(argv, "split_line(): Hit EOF skipping fields before first_col.\n",
			  &vcf_call, line_count, c, first_col, all_sample_ids,
			  "Genotype not relevant");
		exit(EX_DATAERR);
	    }
	
	for (; c <= last_col; ++c)
	{
	    if ( tsv_read_field(argv, vcf_infile, genotype,
				VCF_SAMPLE_MAX_CHARS, &field_len) == EOF )
	    {
		dump_line(argv, "split_line(): Encountered EOF while reading genotype fields.\n",
			  &vcf_call, line_count, c, first_col, all_sample_ids,
			  genotype);
		exit(EX_DATAERR);
	    }
	    if ( selected[c - first_col] )
	    {
		if ( !(flags & FLAG_HET) || (genotype[0] != genotype[2]) )
		{
		    // FIXME: Use vcf_write_ss_call()
		    fprintf(vcf_outfiles[c - first_col],
			"%s\t%s\t.\t%s\t%s\t.\t.\t.\t%s\t%s\n",
			vcf_call.chromosome, vcf_call.pos_str,
			vcf_call.ref, vcf_call.alt, 
			vcf_call.format, genotype);
		}
	    }
	}
	
	if ( tsv_skip_rest_of_line(argv, vcf_infile) == EOF )
	{
	    dump_line(argv, "split_line(): Encountered EOF skipping fields after last_col.\n",
		      &vcf_call, line_count, c, first_col, all_sample_ids,
		      "Genotype not relevant");
	    exit(EX_DATAERR);
	}
	else
	    return 1;
    }
    else
    {
	fprintf(stderr, "%s: split_line(): No more VCF calls.\n", argv[0]);
	fprintf(stderr, "Processed %zu multi-sample VCF calls.\n", line_count);
	return 0;
    }
}


void    dump_line(const char *argv[], const char *message, 
		  vcf_call_t *vcf_call, size_t line_count, size_t col,
		  size_t first_col, const char *all_sample_ids[], char *genotype)

{
    fprintf(stderr, "%s: %s: ", argv[0], message);
    fprintf(stderr, "Line #%zu\n", line_count);
    fprintf(stderr, "Column: %zu Sample ID: %s:\n",
	    col, all_sample_ids[col - first_col]);
    fprintf(stderr, "SS VCF: %s\t%s\t.\t%s\t%s\t.\t.\t.\t%s\t%s\n",
	    vcf_call->chromosome, vcf_call->pos_str,
	    vcf_call->ref, vcf_call->alt,
	    vcf_call->format, genotype);
}


/***************************************************************************
 *  Description:
 *      Read list of sample IDs from a file.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-20  Jason Bacon Begin
 ***************************************************************************/

id_list_t   *read_selected_sample_ids(const char *argv[],
				      const char *samples_file)

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
    fprintf(stderr, "Usage: %s [--het-only] [--max-calls N] [--sample-id-file file] output-file-prefix first-column last-column\n", argv[0]);
    exit(EX_USAGE);
}


void    tag_selected_columns(char *all_sample_ids[],
			     id_list_t *selected_sample_ids, bool selected[],
			     size_t first_col, size_t last_col)

{
    size_t  c;
    
    // No --sample-id-file, select all samples
    if ( selected_sample_ids == NULL )
	for (c = first_col; c <= last_col; ++c)
	    selected[c - first_col] = true;
    else
	for (c = first_col; c <= last_col; ++c)
	{
	    selected[c - first_col] =
		(bsearch(all_sample_ids + c - first_col, selected_sample_ids->ids,
		    selected_sample_ids->count, sizeof(char *),
		    (int (*)(const void *, const void *))strptrcmp)
		    != NULL);
	}
}


int     strptrcmp(const char **p1, const char **p2)

{
    return strcmp(*p1, *p2);
}
