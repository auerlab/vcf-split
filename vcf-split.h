/*
 *  Maybe reduce I/O frequency from pipe.  Experiment with this value.
 */

#define BUFF_SIZE           1024*1024
#define MAX_OUTFILES        200000
#define MAX_SAMPLES         1000000
#define SAMPLE_ID_MAX       256         // Use only for temp variables
#define CMD_MAX             128

typedef struct
{
    size_t  count;
    char    **ids;
}   id_list_t;

typedef unsigned char   flag_t;

#define FLAG_HET    0x01
#define FLAG_XZ     0x02

#include "vcf-split-protos.h"
