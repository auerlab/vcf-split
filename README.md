# vcf-split

## Description

vcf-split splits a combined-sample VCF stream into single-sample VCF files.

Prior methods for splitting a multi-sample VCF stream into single-sample
files involve a loop or parallel job that rereads the multi-sample input for
each sample, e.g. using "bcftools view --samples" to extract one sample at a
time.  This can
become a major bottleneck where there are many samples and/or the input
is compressed.  For example, using "bcftools view" with optimal filtering
options to merely decode one human chromosome BCF with
137,977 samples and pipe the VCF output through "wc" took 12 hours on a
fast server using 2 cores.  To split it into 137,977 single-sample VCFs
would therefore require about 137,977 * 12 * 2 = ~3 million core-hours.
This translates to 171 years on a single core or 125 days using 1000 cores
on an HPC cluster.  This is assuming the cluster filesystem can keep up,
which is doubtful.

vcf-split solves this problem by writing a large number of single-sample VCFs
simultaneously during a single read through the multi-sample input.  The
number of parallel output files is theoretically limited only by the open file
limit of your system, which is typically at least in the tens of thousands on
a modern Unix-like system.  However, to prevent inadvertent system overloads,
a limit of 10,000 samples is hard-coded as a constant.  You must edit the
code and recompile to override this limit.  This should only be done on
your own hardware, to avoid impacting other users.  10,000 samples at
a time should be more than adequate for most purposes.  The 137,977-sample
file mentioned above can be split in 14 runs using this limit.

vcf-split is written entirely in C and attempts to optimize CPU, memory,
and disk access.  It does not inhale large amounts of data into RAM, so memory
use is trivial and it runs mostly from cache, making it very fast.

The example BCF file mentioned above can be split in a few days on a single
server using two cores, with three runs of about 45,000 samples each.
Using the default limit of 10,000 output files did not take much
longer, as there was less competition for disk access.  45,000 was
determined to be near optimal by trial-and-error.

Note that if running multiple vcf-split processes in parallel on a cluster
using the same file server, you may need to reduce the number of samples
significantly.
We found that splitting VCF files for
23 chromosomes simultaneously to the same file
server was limited to about 10,000 samples at once to keep CPU utilization
reasonable.

## Design and Implementation

The code is organized following basic object-oriented design principals, but
implemented in C to minimize overhead and keep the source code accessible to
scientists who don't have time to master the complexities of C++.

Structures are treated as classes, with accessor macros and mutator functions
provided, so dependent applications and libraries need not access
structure members directly.  Since the C language cannot enforce this, it's
up to application programmers to exercise self-discipline.

## Building and installing

vcf-split is intended to build cleanly in any POSIX environment.  Please
don't hesitate to open an issue if you encounter problems on any
Unix-like system.

Primary development is done on FreeBSD with clang, but the code is frequently
tested on Linux, MacOS, NetBSD, and OpenIndiana as well.  MS Windows is not supported,
unless using a POSIX environment such as Cygwin or Windows Subsystem for Linux.

The Makefile is designed to be friendly to package managers, such as
[Debian packages](https://www.debian.org/distrib/packages),
[FreeBSD ports](https://www.freebsd.org/ports/),
[MacPorts](https://www.macports.org/), [dreckly](https://github.com/drecklypkg/dreckly), etc.

End users should install using a package manager, to ensure that
dependencies are properly managed.

I maintain a FreeBSD port and a dreckly package, which is sufficient to install
cleanly on virtually any POSIX platform.  If you would like to see a
package in another package manager, please consider creating a package
yourself.  This will be one of the easiest packages in the collection and
hence a good vehicle to learn how to create packages.

Note that dreckly can be used by anyone, on virtually any POSIX operating
system, with or without administrator privileges.

For an overview of available package managers, see the
[Repology website](https://repology.org/).

### Installing vcf-split on FreeBSD:

FreeBSD is a highly underrated platform for scientific computing, with over
2,000 scientific libraries and applications in the FreeBSD ports collection
(of more than 30,000 total), modern clang compiler, fully-integrated ZFS
filesystem, and renowned security, performance, and reliability.
FreeBSD has a somewhat well-earned reputation for being difficult to set up
and manage compared to user-friendly systems like [Ubuntu](https://ubuntu.com/).
However, if you're a little bit Unix-savvy, you can very quickly set up a
workstation, laptop, or VM using
[desktop-installer](http://www.acadix.biz/desktop-installer.php).
[GhostBSD](https://ghostbsd.org/) offers an experience very similar
to Ubuntu, but is built on FreeBSD rather than Debian Linux.  GhostBSD
packages lag behind FreeBSD ports slightly, but this is not generally
an issue and there are workarounds.

To install the binary package on FreeBSD:

```
pkg install vcf-split
```

You can just as easily build and install from source.  This is useful for
FreeBSD ports with special build options, for building with non-portable
optimizations such as -march=native, and for 
[work-in-progress ports](https://github.com/outpaddling/freebsd-ports-wip),
for which binary packages are not yet maintained.

```
cd /usr/ports/biology/vcf-split && env CFLAGS='-march=native -O2' make install
cd /usr/ports/wip/vcf-split && make install
```

### Installing via dreckly

[Dreckly](https://github.com/drecklypkg/dreckly) is a cross-platform package manager that works on any Unix-like
platform. It is derived from pkgsrc, which is part of [NetBSD](https://www.netbsd.org/), and well-supported
on [Illumos](https://illumos.org/), [MacOS](https://www.apple.com/macos/),
[RHEL](https://www.redhat.com)/[CentOS](https://www.centos.org/), and
many other Linux distributions.
Unlike most package managers, using dreckly does not require admin privileges.  You can install a dreckly
tree in any directory to which you have write access and easily install any
of the nearly 20,000 packages in the collection.

The
[auto-dreckly-setup](https://github.com/outpaddling/auto-admin/blob/master/User-scripts/auto-dreckly-setup)
script will help you install dreckly in about 10 minutes.  Just download it
and run

```
sh auto-dreckly-setup
```

Then, assuming you selected current packages and the default prefix

```
source ~/Dreckly/pkg/etc/dreckly.sh   # Or dreckly.csh for csh or tcsh
cd ~/Dreckly/dreckly/biology/vcf-split
sbmake install clean clean-depends
```

## Instructions for packagers

If you would like to add this project to another package manager
rather than use FreeBSD ports or dreckly, basic manual build instructions
for package can be found
[here](https://github.com/outpaddling/Coding-Standards/blob/main/package.md).
Your contribution is greatly appreciated!
