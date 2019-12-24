#!/bin/sh -e

##########################################################################
#   Script description:
#       
#   History:
#   Date        Name        Modification
#   2019-12-08  Jason Bacon Begin
##########################################################################

# Notes:
#
# Filtering out calls with no het sites is probably useless with 137k
# samples per call.  Every call is likely to have at least one het site.
# Does not seem to work in conjunction with --samples.
# --samples NWDXXXX --genotype het produces lines with 1|1 or 0|0.
# Maybe checking all samples for het instead of desired samples?
#
# vcftools is much slower than bcftools

# From bcf-bench:
# --min-ac + --genotype about doubles run time
# --min-ax + --samples increases run time by order of magnitude

##########################################################################
# vcf-split is limited by the number of open files your filesystem can
# write at once.  Tens of thousands are usually feasible on a high-end
# server.
# Keep in mind that they add up with parallel jobs accessing the same
# file server. Output initially to compute node local disk and use
#
#   #SBATCH --ntasks-per-node=1
#
# to get past this limit.

# vcf-split [--sample-id-file file] first-col last-col
# vcf-split --max-calls N stops after N calls for quick testing.
# Remove it for real runs.

#SBATCH --ntasks=2

bcftools view --min-ac 2 --exclude-types indels \
    ../phased/freeze.8.chr1.pass_only.phased.bcf \
    | ./vcf-split --sample-id-file LLS_WHIMS_NWDID_TOPMED.csv \
    chr1. 1 137977
