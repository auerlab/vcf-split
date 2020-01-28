#!/bin/sh -e

# Example script for compressing finished vcf-split output files while
# vcf-split is running.

while true; do
    for file in *.vcf.done; do
	vcf=${file%.done}
	if [ ! -e $vcf.xz ]; then
	    xz $vcf
	    rm $file
	fi
    fi
done
