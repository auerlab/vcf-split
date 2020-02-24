#!/bin/sh -e

##########################################################################
#   Script description:
#       Compress VCFs output by vcf-split.  This is safe to do
#       periodically while vcf-split is still running, as it only
#       compresses outputs of completed samples.
#       
#   Arguments:
#       Directory containing uncompressed VCF outputs
#       
#   History:
#   Date        Name        Modification
#   2020-02-24  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 vcf-directory\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi

vcf_dir="$1"
for file in $vcf_dir/*.vcf.done; do
    vcf=${file%.done}
    if [ -e $vcf ]; then
	printf "Compressing $vcf...\n"
	ls -l $vcf
	xz $vcf
	ls -l $vcf.xz
    fi
done
