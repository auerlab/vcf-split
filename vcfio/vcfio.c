#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>
#include "tsvio.h"
#include "vcfio.h"

/***************************************************************************
 *  Description:
 *      Skip over header lines in VCF input stream.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

void    vcf_skip_header(const char *argv[], FILE *infile)

{
    char    start[7] = "xxxxxx";
    size_t  count;

    while ( ((count=fread(start, 6, 1, infile)) == 1) && 
	    (memcmp(start, "#CHROM", 6) != 0) )
	tsv_skip_rest_of_line(argv, infile);
    
    // puts(start);
    if ( count == 0 )
    {
	fprintf(stderr, "%s: Error: No #CHROM header found.\n", argv[0]);
	exit(EX_DATAERR);
    }
}


/***************************************************************************
 *  Description:
 *      Extract sample IDs from input header line.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

void    vcf_get_sample_ids(const char *argv[], FILE *infile,
			   char *sample_ids[],
			   size_t first_col, size_t last_col)

{
    size_t  c;
    char    temp_id[VCF_ID_MAX_LEN + 1];
    
    // Skip standard header tags to get to sample IDs
    for (c = 0; c < 9; ++c)
	tsv_skip_field(argv, infile);
    
    // Skip sample IDs before first_col
    for (c = 1; c < first_col; ++c)
	tsv_skip_field(argv, infile);
    
    for (; (c <= last_col) &&
	   tsv_read_field(argv, infile, temp_id, VCF_ID_MAX_LEN) != 0; ++c)
    {
	sample_ids[c - first_col] = strdup(temp_id);
	// fprintf(stderr, "'%s'\n", temp_id);
    }
}


/***************************************************************************
 *  Description:
 *      Read static fields from one line of a single-entry VCF file.
 *      Does not read sample data.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-08  Jason Bacon Begin
 ***************************************************************************/

int     vcf_read_static_fields(const char *argv[],
		      FILE *vcf_stream, vcf_call_t *vcf_call)

{
    char    temp_chromosome[VCF_CHROMOSOME_NAME_MAX + 1],
	    *end;
    
    if ( tsv_read_field(argv, vcf_stream, temp_chromosome, VCF_CHROMOSOME_NAME_MAX) )
    {
	strlcpy(vcf_call->chromosome, temp_chromosome, VCF_CHROMOSOME_NAME_MAX);
	// Call position
	tsv_read_field(argv, vcf_stream, vcf_call->pos_str, VCF_POSITION_MAX_DIGITS);
	vcf_call->pos = strtoul(vcf_call->pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "%s: Invalid call position: %s\n",
		    argv[0], vcf_call->pos_str);
	    exit(EX_DATAERR);
	}
	
	// ID
	tsv_skip_field(argv, vcf_stream);
	
	// Ref
	tsv_read_field(argv, vcf_stream, vcf_call->ref, VCF_REF_NAME_MAX);
	
	// Alt
	tsv_read_field(argv, vcf_stream, vcf_call->alt, VCF_ALT_NAME_MAX);

	// Qual
	tsv_skip_field(argv, vcf_stream);
	
	// Filter
	tsv_skip_field(argv, vcf_stream);
	
	// Info
	tsv_skip_field(argv, vcf_stream);
	
	// Format
	tsv_read_field(argv, vcf_stream, vcf_call->format, VCF_FORMAT_MAX);

#if 0
	fprintf(stderr, "%s %s %s %s %s %s\n",
	    vcf_call->chromosome,
	    vcf_call->pos_str,
	    vcf_call->ref,
	    vcf_call->alt,
	    vcf_call->format);
#endif
	return 1;
    }
    else
	return 0;
}


/***************************************************************************
 *  Description:
 *      Read a single-sample VCF call including genotype.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-11  Jason Bacon Begin
 ***************************************************************************/

int     vcf_read_ss_call(const char *argv[],
		      FILE *vcf_stream, vcf_call_t *vcf_call)

{
    if ( vcf_read_static_fields(argv, vcf_stream, vcf_call) &&
	 tsv_read_field(argv, vcf_stream,
			vcf_call->genotype, VCF_GENOTYPE_NAME_MAX) )
	return 1;
    else
	return 0;
}


/***************************************************************************
 *  Description:
 *      Read in all consecutive calls with the same position.
 *
 *  Returns:
 *      Array of calls and number of calls.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-11  Jason Bacon Begin
 ***************************************************************************/

size_t  vcf_read_duplicate_calls(const char *argv[], FILE *vcf_stream,
				 vcf_duplicate_call_t *vcf_duplicate_calls)

{
    // Cache the next VCF call after the last one returned
    static vcf_call_t   vcf_call;
    static size_t       buffered_calls = 0; // 0 or 1
    size_t              c;
    
    // Prime cache with the first call
    if ( buffered_calls == 0 )
	if ( (buffered_calls = vcf_read_ss_call(argv, vcf_stream, &vcf_call)) == 0 )
	    return 0;
    
    /*
     *  Read all VCF calls with the same position.  The first one with a
     *  different position is not added to vcf_duplicate_calls and is left in
     *  the static vcf_call for the next invocation of this function.
     */
    c = 0;
    do
    {
	vcf_duplicate_calls->vcf_call[c++] = vcf_call;
	buffered_calls = vcf_read_ss_call(argv, vcf_stream, &vcf_call);
    }   while ( (buffered_calls == 1) &&
		(vcf_call.pos == vcf_duplicate_calls->vcf_call[c-1].pos) &&
		(strcmp(vcf_call.chromosome,
			vcf_duplicate_calls->vcf_call[c-1].chromosome) == 0) );

    // Return the number of calls with the same position
    return vcf_duplicate_calls->count = c;
}


#ifdef __linux__
size_t  strlcpy(char * restrict dest, const char * restrictsrc, size_t len)

{
    char   *save_dest, *end;

    save_dest = dest;
    end = (char *)src + len - 1;
    while ((*src != '\0') && (src < end))
	*dest++ = *src++;
    *dest = '\0';
    return dest - save_dest;
}
#endif
