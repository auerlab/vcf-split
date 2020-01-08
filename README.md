# vcf-split
Split combined-sample VCF stream into single-sample VCF files.

Traditional bioinformatics tools need to reread the multi-sample VCF for
each extracted sample, which becomes a major bottleneck for large VCFs,
especially when they are compressed.  For example, using bcftools to simply
decode a single-chromosome BCF file from one of our SRA projects and dumping
the VCF output to /dev/null takes nearly 11 hours:

time bcftools view --no-header freeze.8.chr22.pass_only.phased.bcf > /dev/null

    39486.70 real     39456.34 user        28.55 sys

Repeating this process for the 100,000+ samples in this file would turn it
into a massive and senseless HPC job, wasting cluster resources that can be
put to better use.

vcf-split can write many single-sample VCFs at the same time. It is limited
only by the number of open files supported by your operating system, which
varies by operating system, file system, and hardware specs, but is typically
at least in the tens of thousands.

vcf-split does not inhale data into memory, so memory use is trivial and
it works mostly with cache memory, making it very fast.

Hence, the sample BCF above can be split in about 30 hours on a decent
workstation.

vcf-split is intended to build cleanly in any POSIX environment.  Please
don't hesitate to open an issue if you encounter problems on any
Unix-like system.

Primary development is done on FreeBSD, but the code is frequently tested on
CentOS, Mac OS X, and NetBSD as well.  MS Windows is not supported, unless
using a POSIX environment such as Cygwin or Windows Subsystem for Linux.
