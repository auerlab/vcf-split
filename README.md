# vcf-split
Split combined-sample VCF stream into single-sample VCF files.

Traditional methods for splitting a multi-sample VCF stream into single-sample
files involve a loop or parallel job that rereads the multi-sample input for
each sample.  This can
become a major bottleneck where there are many samples and/or the input
is compressed.  For example, using "bcftools view" with optimal filtering
options to merely decode one human chromosome BCF with
137,977 samples and pipe the VCF output through "wc" took 12 hours on a
fast server using 2 cores.  To split it into 137,977 single-sample VCFs
would therefore require about 137,977 * 12 * 2 = ~3 million core-hours.
This translates to 171 years on a single core or 125 days using 1000 cores
on an HPC cluster.

vcf-split solves this problem by writing a large number of single-sample VCFs
simultaneously during a single read through the multi-sample input.  The
number of parallel output files is theoretically limited only by the open file
limit of your system, which is typically at least in the tens of thousands on
a modern Unix-like system.  However, to prevent inadvertent system overloads,
a limit of 10,000 samples is hard-coded as a constant.  You must edit the
code and recompile to override this limit.

vcf-split is written entirely in C and attempts to optimize CPU, memory,
and disk access.  It does not inhale large amounts of data into RAM, so memory
use is trivial and it runs mostly from cache, making it very fast.

The example BCF file mentioned above can be split in a few days on a single
server using two cores, with three runs of about 45,000 samples each.
Note that if running multiple vcf-split processes in parallel on a cluster
using the same file server, you may need to reduce the number of samples.
We found that splitting all 23 chromosomes simultaneously to the same file
server was limited to about 10,000 samples at once to keep CPU utilization
reasonable.

vcf-split is intended to build cleanly in any POSIX environment.  Please
don't hesitate to open an issue if you encounter problems on any
Unix-like system.

Primary development is done on FreeBSD, but the code is frequently tested on
CentOS, Mac OS X, and NetBSD as well.  MS Windows is not supported, unless
using a POSIX environment such as Cygwin or Windows Subsystem for Linux.

Building and installing:

vcf-split depends on [biolibc](https://github.com/auerlab/biolibc).

Set LOCALBASE to the prefix of lib/libbiolibc.a.  Default is ../local.
(See Makefile).

Set PREFIX to the prefix where you would like to install.  Default is
$LOCALBASE.

Then simply run

```sh
make install
```

vcf-split supports masking off useless fields.  For example, haplohseq only
requires the CHROM, POS, REF, ALT, and FORMAT fields along with the single
sample data.  This can be indicated with

```
vcf-split --output-fields chrom,pos,ref,alt,format
```

All other fields are replaced with a '.', saving considerable space.

vcf-split can also filter for heterozygous sites using --het-only or
sites with at least one alt allele using --alt-only.
