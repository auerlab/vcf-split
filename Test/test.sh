#!/bin/sh -e

pause()
{
    local junk
    
    printf "Press return to continue..."
    read junk
}

cd ..
./cave-man-install.sh
cd Test
../vcf-split test-all-fields- 1 11 < test.vcf
../vcf-split --fields chrom,pos,ref,alt,format test-limited-fields- 1 11 < test.vcf
rm -f *.done

printf "All files should be 12 lines:\n"
wc -l test-*.vcf
pause

printf "There should be no differences shown below:\n"
for col in $(seq 11); do
    diff test-all-fields-$col.vcf correct-all-fields-$col.vcf
    diff test-limited-fields-$col.vcf correct-limited-fields-$col.vcf
done
rm -f test-*.vcf
