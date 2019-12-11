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
#include "vcf-split.h"
#include "tsvio.h"
#include "vcfio.h"

int     main(int argc,const char *argv[])

{
    char        *eos;
    const char  *prefix;
    size_t      first_col,
		last_col,
		max_calls = SIZE_MAX,
		next_arg = 1;
    
    if ( (argc != 4) && (argc != 6) )
	usage(argv);
    
    if ( strcmp(argv[1],"--max-calls") == 0 )
    {
	max_calls = strtoul(argv[2], &eos, 10);
	if ( *eos != '\0' )
	{
	    fprintf(stderr, "%s: %s: Max calls must be a positive integer.\n",
		    argv[0], argv[2]);
	    exit(EX_DATAERR);
	}
	next_arg = 3;
    }
    
    prefix = argv[next_arg++];
    
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
    
    if ( last_col - first_col >= 10000 )
    {
	fprintf(stderr, "%s: Maximum columns is 10,000\n", argv[0]);
	exit(EX_USAGE);
    }

    return split_vcf(argv, stdin, prefix, first_col, last_col, max_calls);
}


void    usage(const char *argv[])

{
    fprintf(stderr, "Usage: %s: [--max-calls N] base-name first-column last-column\n", argv[0]);
    exit(EX_USAGE);
}


/***************************************************************************
 *  Description:
 *      Split a multisample VCF stream into single-sample files.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

int     split_vcf(const char *argv[], FILE *vcf_infile, const char *prefix,
		    size_t first_col, size_t last_col, size_t max_calls)

{
    char    inbuf[BUFF_SIZE + 1],
	    *sample_ids[last_col - first_col + 1];
    
    // Input is likely to come from "bcftools view" stdout.
    // What is optimal buffering for a Unix pipe?  Benchmark several values.
    setvbuf(vcf_infile, inbuf, _IOFBF, BUFF_SIZE);
    vcf_skip_header(argv, vcf_infile);
    vcf_get_sample_ids(argv, vcf_infile, sample_ids, first_col, last_col);
    write_output_files(argv, vcf_infile, (const char **)sample_ids,
			prefix, first_col, last_col, max_calls);
    
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
			    const char *sample_ids[], const char *prefix,
			    size_t first_col, size_t last_col,
			    size_t max_calls)

{
    size_t  columns = last_col - first_col + 1,
	    c;
    FILE    *vcf_outfiles[columns];
    char    filename[PATH_MAX + 1];
    
    // Open all output streams
    for (c = 0; c < columns; ++c)
    {
	snprintf(filename, PATH_MAX, "%s%s.vcf", prefix, sample_ids[c]);
	if ( (vcf_outfiles[c] = fopen(filename, "w")) == NULL )
	{
	    fprintf(stderr, "%s: Error: Cannot create %s.\n", argv[0], filename);
	    exit(EX_CANTCREAT);
	}
	tsv_skip_rest_of_line(argv, vcf_infile);
    }

    // Temporary hack for testing.  Remove limit.
    c = 0;
    while ( (c < max_calls) && split_line(argv, vcf_infile, vcf_outfiles, sample_ids, first_col, last_col) )
	++c;
    
    // Close all output streams
    for (c = 0; c < columns; ++c)
	fclose(vcf_outfiles[c]);
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
		    const char *sample_ids[], size_t first_col, size_t last_col)

{
    size_t      c;
    vcf_call_t  vcf_call;
    char        genotype[VCF_GENOTYPE_NAME_MAX + 1];
    
    /*
     *  Read in VCF fields
     */
    
    if ( vcf_read_call(argv, vcf_infile, &vcf_call) )
    {
	// Skip columns before first_col
	for (c = 1; c < first_col; ++c)
	    tsv_skip_field(argv, vcf_infile);
	
	for (; c <= last_col; ++c)
	{
	    tsv_read_field(argv, vcf_infile, genotype, VCF_GENOTYPE_NAME_MAX);
	    if ( genotype[0] != genotype[2] )
	    {
#if DEBUG
		fprintf(stderr, "%zu %s:\n", c, sample_ids[c-1]);
		fprintf(stderr, "%s\t%s\t.\t%s\t%s\t.\t.\t.\t%s\t%s\n",
		    vcf_call.chromosome, vcf_call.pos_str,
		    vcf_call.ref, vcf_call.alt,
		    vcf_call.format, genotype);
#endif
		fprintf(vcf_outfiles[c-1], "%s\t%s\t.\t%s\t%s\t.\t.\t.\t%s\t%s\n",
		    vcf_call.chromosome, vcf_call.pos_str,
		    vcf_call.ref, vcf_call.alt, 
		    vcf_call.format, genotype);
	    }
	}
	
	tsv_skip_rest_of_line(argv, vcf_infile);
	return 1;
    }
    else
	return 0;
}

