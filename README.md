# vcf-split
Split combined-sample VCF stream into single-sample VCF files.

It is currently suitable for generating inputs for haplohseq, but the
intention is to generalize it in the future.

Traditional methods for splitting a multi-sample VCF stream into single-sample
files involve a loop or parallel job that rereads the multi-sample input for
every sample.  This can
become a major bottleneck where there are many samples and/or the input
is compressed.  For example, using "bcftools view" with optimal filtering
options to merely decode one human chromosome BCF with
137,977 samples and pipe the VCF output through "wc" took 12 hours on a
fast server using 2 cores.  To split it into 137,977 single-sample VCFs
would therefore require about 137,977 * 12 * 2 = ~3 million core-hours.
This translates to 171 years on a single server or 125 days using 1000 cores
on an HPC cluster.

vcf-split solves this problem by writing a large number of single-sample VCFs
simultaneously during a single read of the multi-sample input.  The number
of parallel output files is limited only by the open file limit of your
system, which is typically at least in the tens of thousands on a modern
Unix workstation or server.

vcf-split is written entirely in C and attempts to optimize CPU, memory,
and disk access.  It does not inhale large amounts of data into RAM, so memory
use is trivial and it runs mostly from cache, making it very fast.

The example BCF file mentioned above can be split in a few days on a single
server using two cores, with three runs of about 45,000 output files each.

vcf-split is intended to build cleanly in any POSIX environment.  Please
don't hesitate to open an issue if you encounter problems on any
Unix-like system.

Primary development is done on FreeBSD, but the code is frequently tested on
CentOS, Mac OS X, and NetBSD as well.  MS Windows is not supported, unless
using a POSIX environment such as Cygwin or Windows Subsystem for Linux.
