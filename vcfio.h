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
	    format[VCF_FORMAT_MAX + 1];
    size_t  pos;
}   vcf_call_t;

#include "vcfio-protos.h"
