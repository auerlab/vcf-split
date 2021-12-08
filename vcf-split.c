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
#include <xtend/dsv.h>
#include <xtend/string.h>
#include <biolibc/vcf.h>
#include "vcf-split.h"

int     main(int argc, char *argv[])

{
    char        *field_spec,
		*eos;
    const char  *outfile_prefix,
		*selected_samples_file = NULL;
    id_list_t   *selected_sample_ids = NULL;
    size_t      first_col,
		last_col,
		max_calls = SIZE_MAX;
    int         next_arg = 1;
    flag_t      flags = 0;
    // Overridden if specified on command line
    vcf_field_mask_t    field_mask = BL_VCF_FIELD_ALL;
    
    if ( (argc == 2) && (strcmp(argv[1],"--version")) == 0 )
    {
	printf("vcf-split %s\n", VERSION);
	return EX_OK;
    }
    
    next_arg = 1;
    while ( (next_arg < argc ) && (argv[next_arg][0] == '-') )
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

	else if ( strcmp(argv[next_arg], "--alt-only") == 0 )
	{
	    if ( flags & FLAG_HET )
	    {
		fprintf(stderr,
		    "%s: --het-only and --alt-only are mutually exclusive.\n",
		    argv[0]);
		usage(argv);
	    }
	    flags |= FLAG_ALT;
	    ++next_arg;
	}

	else if ( strcmp(argv[next_arg], "--fields") == 0 )
	{
	    ++next_arg;
	    field_spec = argv[next_arg++];
	    if ( (field_mask = bl_vcf_parse_field_spec(field_spec))
		    == BL_VCF_FIELD_ERROR )
		usage(argv);
	}
	
	else
	    usage(argv);
    }

    if ( next_arg > argc - 3 )
	usage(argv);

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
    
    return vcf_split(argv, stdin, outfile_prefix, first_col, last_col,
		     selected_sample_ids, max_calls, flags, field_mask);
}


/***************************************************************************
 *  Description:
 *      Split a multisample VCF stream into single-sample files.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

int     vcf_split(char *argv[], FILE *vcf_infile,
		  const char *outfile_prefix,
		  size_t first_col, size_t last_col,
		  id_list_t *selected_sample_ids, size_t max_calls,
		  flag_t flags, vcf_field_mask_t field_mask)

{
    char    inbuf[BUFF_SIZE + 1],
	    *all_sample_ids[last_col - first_col + 1];
    bool    selected[last_col - first_col + 1];
    size_t  c;
    FILE    *meta_stream;
    
    // Input is likely to come from "bcftools view" stdout.
    // What is optimal buffering for a Unix pipe?  Benchmark several values.
    setvbuf(vcf_infile, inbuf, _IOFBF, BUFF_SIZE);
    if ( bl_vcf_skip_meta_data(vcf_infile, &meta_stream) != BL_READ_OK )
	exit(EX_DATAERR);
    bl_vcf_get_sample_ids(vcf_infile, all_sample_ids, first_col, last_col);

    /*
    fputs("All sample IDs:", stderr);
    for (c = 0; c < last_col - first_col + 1; ++c)
	fprintf(stderr, " %s", all_sample_ids[c]);
    putc('\n', stderr);
    */
    
    tag_selected_columns(all_sample_ids, selected_sample_ids, selected,
			 first_col, last_col);

    if ( selected_sample_ids != NULL )
    {
	fputs("Tagged sample IDs:", stderr);
	for (c = 0; c < selected_sample_ids->count; ++c)
	    fprintf(stderr, " %s", selected_sample_ids->ids[c]);
	fputc('\n', stderr);
    }
    
    write_output_files(argv, vcf_infile, meta_stream, 
		       (const char **)all_sample_ids,
		       selected, outfile_prefix,
		       first_col, last_col, max_calls, flags, field_mask);
    
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

void    write_output_files(char *argv[], FILE *vcf_infile, FILE *header,
			    const char *all_sample_ids[], bool selected[],
			    const char *outfile_prefix,
			    size_t first_col, size_t last_col,
			    size_t max_calls, flag_t flags,
			    vcf_field_mask_t field_mask)

{
    size_t  columns = last_col - first_col + 1,
	    c;
    int     fd;
    FILE    *vcf_outfiles[columns];
    char    filename[PATH_MAX + 1],
	    file_format[129];
    
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
	    
	    /*
	     *  Add basic header
	     *  FIXME: Add option to copy all/part of source header
	     */
	    
	    rewind(header);
	    fgets(file_format, 128, header);
	    if ( memcmp(file_format, "##fileformat", 12) == 0 )
		fputs(file_format, vcf_outfiles[c]);
	    fputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
		  vcf_outfiles[c]);
	}
    }

    // Heart of the program, split each VCF line across multiple files
    for (c = 0; split_line(argv, vcf_infile, vcf_outfiles, all_sample_ids,
			   selected, first_col, last_col, max_calls, flags,
			   field_mask);
			   ++c)
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

int     split_line(char *argv[], FILE *vcf_infile, FILE *vcf_outfiles[],
		   const char *all_sample_ids[], bool selected[],
		   size_t first_col, size_t last_col, size_t max_calls,
		   flag_t flags, vcf_field_mask_t field_mask)

{
    static size_t   line_count = 0,
		    max_info_len = 0;
    size_t          c,
		    field_len;
    int             delimiter = ' ';    // Silence false uninit warning
    static bl_vcf_t vcf_call;  // Use bl_vcf_init() function to initialize
    char            *genotype;
    
    /*
     *  Read in VCF fields
     */
    
    /* Declared as static: Allocate only once and reuse */
    if ( vcf_call.single_sample == NULL )
	bl_vcf_init(&vcf_call, VCF_INFO_MAX_CHARS, VCF_FORMAT_MAX_CHARS,
		  VCF_SAMPLE_MAX_CHARS);
    genotype = vcf_call.single_sample;
    
    // Check max_calls here rather than outside in order to print the
    // end-of-run report below
    if ( (line_count < max_calls) && 
	 (bl_vcf_read_static_fields(vcf_infile, &vcf_call, field_mask) == BL_READ_OK) )
    {
	if ( (++line_count % 100 == 0) && isatty(fileno(stderr)) )
	    fprintf(stderr, "%zu\r", line_count);
	
	if ( BL_VCF_INFO_LEN(&vcf_call) > max_info_len )
	    max_info_len = BL_VCF_INFO_LEN(&vcf_call);
	
	// Skip columns before first_col
	for (c = 1; c < first_col; ++c)
	{
	    delimiter = tsv_skip_field(vcf_infile);
	    if ( delimiter == EOF )
	    {
		dump_line(argv, "split_line(): Hit EOF skipping fields before first_col.\n",
			  &vcf_call, line_count, c, first_col, all_sample_ids,
			  "Genotype not relevant");
		exit(EX_DATAERR);
	    }
	    else if ( delimiter == '\n' )
	    {
		fprintf(stderr, "%s: split_line(): Reached EOL before first_col.\n", argv[0]);
		fprintf(stderr, "Does your input really have %zu samples?\n",
			first_col);
		usage(argv);
	    }
	}
	
	for (; (c <= last_col) && 
	       (delimiter = tsv_read_field(vcf_infile, genotype,
				VCF_SAMPLE_MAX_CHARS, &field_len)) != EOF;
				++c)
	{
	    if ( selected[c - first_col] )
	    {
		if ( (flags == FLAG_NONE) ||
		     ((flags == FLAG_HET) && (genotype[0] != genotype[2])) ||
		     ((flags == FLAG_ALT) &&
		      ((genotype[0] == '1') || (genotype[2] == '1'))) )
		{
		    bl_vcf_write_ss_call(vcf_outfiles[c - first_col],
					    &vcf_call, BL_VCF_FIELD_ALL);
		    /*
		    fprintf(vcf_outfiles[c - first_col],
			    "%s\t%s\t%s\t%s\t%s\t.\t.\t.\t%s\t%s\n",
			    BL_VCF_CHROM(&vcf_call), BL_VCF_POS_STR(&vcf_call),
			    BL_VCF_ID(&vcf_call), BL_VCF_REF(&vcf_call),
			    BL_VCF_ALT(&vcf_call), BL_VCF_FORMAT(&vcf_call),
			    genotype);
		    */
		}
	    }
	}
	
	if ( delimiter == EOF )
	{
	    dump_line(argv, "split_line(): Encountered EOF while reading genotype fields.\n",
		      &vcf_call, line_count, c, first_col, all_sample_ids,
		      genotype);
	    exit(EX_DATAERR);
	}

	if ( (delimiter == '\n') && (c < last_col) )
	{
	    fprintf(stderr, "%s: split_line(): Reached EOL before last_col.\n", argv[0]);
	    fprintf(stderr, "Does your input really have %zu samples?\n", last_col);
	    usage(argv);
	}

	// If this wasn't the last sample in the line, skip to EOL
	if ( delimiter != '\n' )
	{
	    if ( tsv_skip_rest_of_line(vcf_infile) == EOF )
	    {
		dump_line(argv, "split_line(): Encountered EOF skipping fields after last_col.\n",
			  &vcf_call, line_count, c, first_col, all_sample_ids,
			  "Genotype not relevant");
		exit(EX_DATAERR);
	    }
	}
	return 1;
    }
    else
    {
	fprintf(stderr, "%s: split_line(): No more VCF calls.\n", argv[0]);
	fprintf(stderr, "Processed %zu multi-sample VCF calls.\n", line_count);
	fprintf(stderr, "Max info_len = %zu.\n", max_info_len);
	return 0;
    }
}


void    dump_line(char *argv[], const char *message, 
		  bl_vcf_t *vcf_call, size_t line_count, size_t col,
		  size_t first_col, const char *all_sample_ids[], char *genotype)

{
    fprintf(stderr, "%s: %s\n", argv[0], message);
    fprintf(stderr, "Input line: %zu\n", line_count);
    fprintf(stderr, "Column: %zu Sample ID: %s:\n",
	    col, all_sample_ids[col - first_col]);
    fprintf(stderr, "SS VCF: %s\t%" PRIu64 "\t.\t%s\t%s\t.\t.\t.\t%s\t%s\n",
	    BL_VCF_CHROM(vcf_call), BL_VCF_POS(vcf_call),
	    BL_VCF_REF(vcf_call), BL_VCF_ALT(vcf_call),
	    BL_VCF_FORMAT(vcf_call), genotype);
}


/***************************************************************************
 *  Description:
 *      Read list of sample IDs from a file.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-20  Jason Bacon Begin
 ***************************************************************************/

id_list_t   *read_selected_sample_ids(char *argv[],
				      const char *samples_file)

{
    id_list_t   *list;
    size_t      c;
    char        temp_id[BL_VCF_ID_MAX_CHARS + 1];
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

    for (list->count = 0; read_string(fp, temp_id, BL_VCF_ID_MAX_CHARS) > 0; ++list->count )
	;
    
    if ( (list->ids = (char **)malloc(list->count * sizeof(char **))) == NULL )
    {
	fprintf(stderr, "%s: Cannot allocate sample_ids list.\n", argv[0]);
	exit(EX_UNAVAILABLE);
    }
    
    rewind(fp);
    for (c = 0; c < list->count; ++c)
    {
	read_string(fp, temp_id, BL_VCF_ID_MAX_CHARS);
	if ( (list->ids[c] = strdup(temp_id)) == NULL )
	{
	    fprintf(stderr, "%s: Cannot allocate sample id %zu.\n", argv[0], c);
	    exit(EX_UNAVAILABLE);
	}
    }
    fclose(fp);
    
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


void    usage(char *argv[])

{
    fprintf(stderr, "\nUsage: %s\n\t[--version]\n", argv[0]);
    fprintf(stderr, "\nUsage: %s\n\t[--het-only]\n\t[--alt-only]\n\t"
		    "[--max-calls N]\n\t[--sample-id-file file]\n\t"
		    "[--fields field-spec]\n\toutput-file-prefix\n\t"
		    "first-column\n\tlast-column\n\n", argv[0]);
    fputs("Press return to continue...", stderr);
    getchar();
    fprintf(stderr, "\n--het-only indicates that only heterozygous fields are output.\n"
		    "Allowing vcf-split to perform this filtering is faster than doing\n"
		    "it in bcftools in some cases.\n\n"
		    "--alt-only indicates that only fields with at least one alt allele are output.\n\n"
		    "--max-calls limits the number of calls processed (for testing purposes).\n\n"
		    "--sample-id-file indicates a list of samples to extract.  Names must\n"
		    "match the column header in the input VCF.\n\n"
		    "field-spec is a comma-separated list of fields to include in the output\n"
		    "including one or more of\n\n"
		    "chrom,pos,id,ref,alt,qual,filter,info\n\n"
		    "fields not indicated are replaced with a '.'.\n\n"
		    "output-file-prefix is prepended to .vcf\n"
		    "It may include one or more subdirectories\n\n"
		    "first-column and last column indicate the range of samples to process\n"
		    "Output is the intersection of this range and --sample-id-file.\n\n"
		    );
    exit(EX_USAGE);
}


void    tag_selected_columns(char *all_sample_ids[],
			     id_list_t *selected_sample_ids, bool selected[],
			     size_t first_col, size_t last_col)

{
    size_t  c,
	    total_selected;
    
    // No --sample-id-file, select all samples
    if ( selected_sample_ids == NULL )
    {
	for (c = first_col; c <= last_col; ++c)
	    selected[c - first_col] = true;
	total_selected = last_col - first_col + 1;
    }
    else
    {
	for (c = first_col, total_selected = 0; c <= last_col; ++c)
	{
	    selected[c - first_col] =
		(bsearch(all_sample_ids + c - first_col, selected_sample_ids->ids,
		    selected_sample_ids->count, sizeof(char *),
		    (int (*)(const void *, const void *))strptrcmp)
		    != NULL);
	    if ( selected[c - first_col] )
		++total_selected;
	}
    }

    if ( total_selected > MAX_OUTFILES )
    {
	fprintf(stderr,
		"vcf-split:  Maximum output files = %u.  You have selected %zu.\n"
		"Increase MAX_OUTFILES and recompile if a %ux speedup is really not enough\n"
		"and you want to bludgeon your system.  Otherwise, you can process the\n"
		"multi-sample input file in %u-sample chunks.\n",
		MAX_OUTFILES, total_selected, MAX_OUTFILES, MAX_OUTFILES);
	exit(EX_SOFTWARE);
    }
}
