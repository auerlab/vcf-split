// FIXME: Are there limits defined by the VCF format?
#define VCF_ID_MAX_LEN          256
#define VCF_CHROMOSOME_NAME_MAX 256
#define VCF_POSITION_MAX_DIGITS 32
#define VCF_REF_NAME_MAX        256
#define VCF_ALT_NAME_MAX        256
#define VCF_GENOTYPE_NAME_MAX   256
#define VCF_FORMAT_MAX          4096

typedef struct
{
    char    chromosome[VCF_CHROMOSOME_NAME_MAX + 1],
	    pos_str[VCF_POSITION_MAX_DIGITS + 1],
	    ref[VCF_REF_NAME_MAX + 1],
	    alt[VCF_ALT_NAME_MAX + 1],
	    format[VCF_FORMAT_MAX + 1],
	    genotype[VCF_GENOTYPE_NAME_MAX + 1];    // Only used for SS calls
    size_t  pos;
}   vcf_call_t;

typedef struct
{
    size_t      count;
    vcf_call_t  vcf_call[10];   // FIXME: What's the real maximum?
}   vcf_duplicate_call_t;

#ifdef __linux__
size_t strlcpy(char *restrict dest, const char *restrict src, size_t len);
#endif

#include "vcfio-protos.h"
