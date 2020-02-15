#!/bin/sh -e

# Example script for compressing finished vcf-split output files while
# vcf-split is running.

# FIXME: Loop until all output VCFs are compressed
while true; do
    for file in *.vcf.done; do
	vcf=${file%.done}
	if [ ! -e $vcf.xz ]; then
	    printf "Compressing $vcf...\n"
	    xz $vcf
	    rm $file
	fi
    done
done
