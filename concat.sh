#!/bin/sh -e

for sample in $(cat 10-samples.txt); do
    echo $sample
    c=1
    files=""
    
    # Filenames don't sort lexically in line with chromosome number
    # E.g. chr2. > chr10.
    # Construct the argument list in order from 1 to 22
    while [ $c -le 22 ]; do
	files="$files chr$c.$sample.vcf.xz"
	c=$((c + 1))
    done
    
    # bcftools requires headers in VCF inputs
    # outfile=combined.$sample.vcf
    # bcftools concat $files --output-type=v --output=$outfile
    
    outfile=combined.$sample.vcf.xz
    rm -f $outfile
    xzcat $files | xz -c > $outfile
done
