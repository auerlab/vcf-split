/*
 *  Maybe reduce I/O frequency from pipe.  Experiment with this value.
 */

#define BUFF_SIZE           1024*1024
#define MAX_OUTFILES        10000
#define SAMPLE_ID_MAX       256         // Use only for temp variables
#define CMD_MAX             128

// Match these with ad2vcf
// Yes, we actually saw a few INFO fields over 512k in some dbGap BCFs
#define BL_VCF_INFO_MAX_CHARS          1048576
#define BL_VCF_FORMAT_MAX_CHARS        4096
#define BL_VCF_SAMPLE_MAX_CHARS        2048

// FIXME: Encapsulate this and move to biolibc?
typedef struct
{
    size_t  count;
    char    **ids;
}   id_list_t;

typedef unsigned char   flag_t;

#define FLAG_NONE   0x0
#define FLAG_HET    0x1
#define FLAG_ALT    0x2

#include "vcf-split-protos.h"
