# QuasiFit 0.1
David Seifert (david.seifert@bsse.ethz.ch)
Niko Beerenwinkel (niko.beerenwinkel@bsse.ethz.ch)

## Introduction
QuasiFit is an MCMC sampler that implements (relative) fitness inference for NGS data assuming a mutation-selection equilibrium. From the posterior, conclusions such as determining neutral networks and detecting epistasis can be drawn.

## Binaries
We have pre-compiled binaries for Linux and Mac users:

* Linux: `quasifit-linux-static-amd64`

  You will require a distribution having _at least_ glibc 2.3.2. Any distribution from the past 10 years should work.  
  The linux binary was built on Debian Etch 4.0r9 64-bit.

* Mac: `quasifit-mac-static-amd64`

  You will require at least Mac OS X 10.6.8.

QuasiFit was built on both platforms with GSL 1.16 and Boost 1.55 with GCC 4.8.2 on -O2 optimizations. The main code (including Eigen) was compiled with -O3 optimizations. All libraries, including the C++ runtime libraries, have been linked statically to produce a binary that has **no external dependencies**.

Download static binaries from:
http://github.com/SoapZA/QuasiFit/releases

## Prerequisites
If you wish to compile QuasiFit from source, you will require the following components (these are **not** necessary for running the statically linked binary):

1. GNU Scientific Library (GSL); _somewhat recent_ release (http://www.gnu.org/software/gsl/)

   The GNU Scientific Library is required for random number generating functions.

2. Eigen; _at least 3.1_ (http://eigen.tuxfamily.org/)

   Eigen forms the core mathematics library of QuasiFit, with all its linear algebra routines.

3. Boost; _at least 1.50_ (http://www.boost.org/)

   Boost provides the necessary abstraction for time routines and thread handling.

Furthermore, you will require a compiler that can handle **C++0x** (which includes all C++11 compilers). QuasiFit has been successfully compiled with GCC 4.4 on RHEL 6, GCC 4.8 on Gentoo and icc 12.0 on RHEL 6.

If you wish to do development, you will require parts of the extended GNU toolchain (the infamous Autotools):

1. GNU Autoconf; latest 2.69 release (http://www.gnu.org/software/autoconf/)

   Autoconf produces the ./configure script from configure.ac.

2. GNU Automake; latest 1.14 release (http://www.gnu.org/software/automake/)

   Automake produces the Makefile.in precursor, that is processed with ./configure to yield the final Makefile.

3. GNU Libtool; latest 2.4.2 release (http://www.gnu.org/software/automake/)

   Libtool is required as a dependency of boost.m4.

## Building
QuasiFit can be used as a statically compiled binary without requiring any dependencies. Should you wish to to compile QuasiFit yourself, run

```
git clone https://github.com/SoapZA/QuasiFit.git
cd QuasiFit/
./configure
make
```

The resulting executable `quasifit` can then be run.

## Options
The QuasiFit sampler has a number of options for controlling the Metropolis-Hastings algorithm and I/O. See `quasifit -h` for more information.
