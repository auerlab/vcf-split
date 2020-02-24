#!/bin/sh -e

##########################################################################
#   Script description:
#       Concatenate single-chromosome VCFs into a full-genome VCF for
#       each sample.  Run this only after all vcf-split jobs have
#       finished.
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
for sample in $(cat 10-samples.txt); do
    # Assumes files list in order of chromosomes.
    # vcf-split generates filenames with chr01, chr02, ... chr10, etc.
    # to ensure this.  If using VCFs from another tool that ouputs
    # chr1, chr2, chr10, ..., this will put 10 before 2
    files="$vcf_dir/chr*$sample*.vcf.xz"
    
    # bcftools requires headers in VCF inputs
    # outfile=combined.$sample.vcf
    # bcftools concat $files --output-type=v --output=$outfile
    
    outfile=$vcf_dir/combined.$sample.vcf.xz
    printf "Concatenating $files to $outfile...\n"
    rm -f $outfile
    xzcat $files | xz -c > $outfile
done
