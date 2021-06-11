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

## Design and Implementation

The code is organized following basic object-oriented design principals, but
implemented in C to minimize overhead and keep the source code accessible to
scientists who don't have time to master the complexities of C++.

Structures are treated as classes, with accessor and mutator functions
(or macros) provided, so dependent applications and libraries need not access
structure members directly.  Since the C language cannot enforce this, it's
up to application programmers to exercise self-discipline.

## Building and installing

vcf-split is intended to build cleanly in any POSIX environment.  Please
don't hesitate to open an issue if you encounter problems on any
Unix-like system.

Primary development is done on FreeBSD with clang, but the code is frequently
tested on CentOS, MacOS, and NetBSD as well.  MS Windows is not supported,
unless using a POSIX environment such as Cygwin or Windows Subsystem for Linux.

The Makefile is designed to be friendly to package managers, such as
[Debian packages](https://www.debian.org/distrib/packages),
[FreeBSD ports](https://www.freebsd.org/ports/),
[MacPorts](https://www.macports.org/), [pkgsrc](https://pkgsrc.org/), etc.
End users should install via one of these if at all possible.

I maintain a FreeBSD port and a pkgsrc package.

### Installing vcf-split on FreeBSD:

FreeBSD is a highly underrated platform for scientific computing, with over
1,800 scientific libraries and applications in the FreeBSD ports collection
(of more than 30,000 total), modern clang compiler, fully-integrated ZFS
filesystem, and renowned security, performance, and reliability.
FreeBSD has a somewhat well-earned reputation for being difficult to set up
and manage compared to user-friendly systems like [Ubuntu](https://ubuntu.com/).
However, if you're a little bit Unix-savvy, you can very quickly set up a
workstation, laptop, or VM using
[desktop-installer](http://www.acadix.biz/desktop-installer.php).  If
you're new to Unix, you can also reap the benefits of FreeBSD by running
[GhostBSD](https://ghostbsd.org/), a FreeBSD distribution augmented with a
graphical installer and management tools.  GhostBSD does not offer as many
options as desktop-installer, but it may be more comfortable for Unix novices.

```
pkg install vcf-split
```

### Installing via pkgsrc

pkgsrc is a cross-platform package manager that works on any Unix-like
platform. It is native to [NetBSD](https://www.netbsd.org/) and well-supported
on [Illumos](https://illumos.org/), [MacOS](https://www.apple.com/macos/),
[RHEL](https://www.redhat.com)/[CentOS](https://www.centos.org/), and
many other Linux distributions.
Using pkgsrc does not require admin privileges.  You can install a pkgsrc
tree in any directory to which you have write access and easily install any
of the nearly 20,000 packages in the collection.  The
[auto-pkgsrc-setup](http://netbsd.org/~bacon/) script can assist you with
basic setup.

First bootstrap pkgsrc using auto-pkgsrc-setup or any
other method.  Then run the following commands:

```
cd pkgsrc-dir/biology/vcf-split
bmake install clean
```

There may also be binary packages available for your platform.  If this is
the case, you can install by running:

```
pkgin install vcf-split
```

See the [Joyent Cloud Services Site](https://pkgsrc.joyent.com/) for
available package sets.

### Building vcf-split locally

Below are cave man install instructions for development purposes, not
recommended for regular use.

vcf-split depends on [biolibc](https://github.com/auerlab/biolibc).
Install biolibc before attempting to build vcf-split.

1. Clone the repository
2. Run "make depend" to update Makefile.depend
3. Run "make install"

The default install prefix is ../local.  Clone vcf-split, biolibc and dependent
apps into sibling directories so that ../local represents a common path to all
of them.

To facilitate incorporation into package managers, the Makefile respects
standard make/environment variables such as CC, CFLAGS, PREFIX, etc.  

Add-on libraries required for the build, such as biolibc, should be found
under , which defaults to ../local.
The library, headers, and man pages are installed under
.  DESTDIR is empty by default and is primarily used by
package managers to stage installations.  PREFIX defaults to .

To install directly to /myprefix, assuming biolibc is installed there as well,
using a make variable:

```
make LOCALBASE=/myprefix clean depend install
```

Using an environment variable:

```
# C-shell and derivatives
setenv LOCALBASE /myprefix
make clean depend install

# Bourne shell and derivatives
LOCALBASE=/myprefix
export LOCALBASE
make clean depend install
```

View the Makefile for full details.
