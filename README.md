# QuasiFit 0.1
David Seifert (david.seifert@bsse.ethz.ch)
Niko Beerenwinkel (niko.beerenwinkel@bsse.ethz.ch)

## Introduction
QuasiFit is an MCMC sampler that implements (relative) fitness inference for NGS data assuming a mutation-selection equilibrium. From the posterior, conclusions such as determining neutral networks and detecting epistasis can be drawn.

## Binaries
We have pre-compiled 64-bit binaries for Linux and Mac users:

* Linux: `quasifit-linux-static-amd64`

  You will require a distribution having _at least_ glibc 2.3.2. Any distribution from the past 10 years should work.  
  The linux binary was built on Debian Etch 4.0r9 64-bit.

* Mac: `quasifit-mac-static-amd64`

  You will require at least Mac OS X 10.6.8.

QuasiFit was built on both platforms with GSL 1.16 and Boost 1.55 with GCC 4.8.2 on -O2 optimizations. The main code (including Eigen) was compiled with -O3 optimizations. All libraries, including the C++ runtime libraries, have been linked statically to produce a binary that has **no external dependencies**.

All static binaries are part of the main source tarball and can be found in the directory `binaries/` of the uncompressed `quasifit-0.1/` directory. Additionally, all binaries can also be downloaded from the main git tree for the most recent release.

## Prerequisites
If you wish to compile QuasiFit from source, you will require the following components (these are **not** necessary for running the statically linked binary):

1. GNU Scientific Library (GSL); _somewhat recent_ release (http://www.gnu.org/software/gsl/)

   The GNU Scientific Library is required for random number generating functions.

2. Eigen; _at least 3.1_ (http://eigen.tuxfamily.org/)

   Eigen forms the core mathematics library of QuasiFit, with all its linear algebra routines.

3. Boost; _at least 1.50_ (http://www.boost.org/)

   Boost provides the necessary abstraction for time routines and thread handling.

4. GMP (optional); _somewhat recent_ release (http://www.gmplib.org/)

   The GNU Multiple Precision Arithmetic Library (GMP) provides the basis for arbitrary precision floating-point calculations. It is only required if you wish to build an arbitrary precision sampler.

5. libquadmath (optional); _at least GCC 4.6_ (http://gcc.gnu.org/onlinedocs/libquadmath/)

   GCC's libquadmath provides the __float128 quad-precision floating point type and associated operations. This is an internal GCC library that is included with GCC since 4.6. It is only required if you wish to use a quad-precision sampler.

Furthermore, you will require a compiler that can handle **C++0x** (which includes all C++11 compilers). QuasiFit has been successfully compiled with GCC 4.4 on RHEL 6, GCC 4.8 on Gentoo/Debian Etch and icc 12.0 on RHEL 6.

If you wish to do development, you will require parts of the extended GNU toolchain (the infamous Autotools):

1. GNU Autoconf; latest 2.69 release (http://www.gnu.org/software/autoconf/)

   Autoconf produces the ./configure script from configure.ac.

2. GNU Automake; latest 1.14 release (http://www.gnu.org/software/automake/)

   Automake produces the Makefile.in precursor, that is processed with ./configure to yield the final Makefile.

3. GNU Libtool; latest 2.4.2 release (http://www.gnu.org/software/automake/)

   Libtool is required as a dependency of boost.m4.

QuasiFit is strongly intertwined with libraries and programs that heavily rely on features of UNIX-like systems, hence supporting Microsoft Windows is not a goal (in particular, building the GNU Scientific Library and using the GNU build system on Windows is a nightmare).

## Preparing
To install the aforementioned dependencies, we provide some guidance here

### Linux
Due to the large heterogeneity of the Linux distributions landscape, we will only detail the procedure of installing dependencies for Ubuntu 14.04 LTS here. The procedure should be very similiar for Debian.
```
# install basic compiler toolchain (you will be prompted to enter your password)
sudo apt-get install build-essential

# install GSL
sudo apt-get install libgsl0-dev

# install Eigen
sudo apt-get install libeigen3-dev

# install Boost
sudo apt-get install libboost-all-dev

# (optional) install GMP for arbitrary precision arithmetic
sudo apt-get install libgmp-dev
```

If you wish to work with the bleeding-edge release of QuasiFit, you will need the complete GNU Autotools toolchain. It should be reiterated here that the recommended way of building QuasiFit is by downloading the provided tarball and using either the included static binaries or compiling from source. The Git tree needs to be bootstrapped to produce the various scripts. To install the Autotools:
```
# install the GNU toolchain
sudo apt-get install autoconf automake libtool pkg-config git
```

### Mac OS X
QuasiFit should be buildable without complications on all Mac OS X versions above and including 10.6. For older versions of Mac OS X, the build process is significantly more involved due to the C++11 requirement. In this case we recommend using the provided precompiled binaries.

In any case, you will need to install the latest version of

1. Xcode for your platform (4.2 for 10.6; 4.6.3 for 10.7; 5.1.x for 10.8 & 10.9) either via the Mac App Store or by downloading the disk image from the Apple Developer Connection (http://developer.apple.com/).

2. Command Line Tools for Xcode. Since Xcode 4.3 Apple has stopped shipping command line tools with the standard Xcode package. You will need to install these via "Downloads" in the "Xcode" -> "Preferences" menu, or (preferably) by downloading the latest appropriate "Command Line Tools" package from the Apple Developer Connection.

3. MacPorts (http://www.macports.org/install.php).

Henceforth, we assume both Xcode, the Command Line Tools and Macports to be installed.

Install the remaining libraries from MacPorts by performing
```
# GCC (optional; if you wish to use quad-precision on Mac OS X,
# you will require GCC as Clang/LLVM cannot handle libquadmath)
sudo port install gcc48

# install GSL; choose one of the following
# standard variant:
sudo port install gsl
# optimized variant, compiled with -O3 and -march=native:
sudo port install gsl +optimize

# install Eigen
sudo port install eigen3

# install Boost; choose one of the following
# standard variant, pulls in a load of dependencies:
sudo port install boost
# minimalist variant,
# disables python, avoids multiple dependencies:
sudo port install boost -python27
# minimalist variant,
# also builds static libraries that can be linked into the
# final executable to reduce dynamic linking,
# making the executable more portable:
sudo port install boost -python27 -no_static

# (optional) install GMP for arbitrary precision arithmetic
sudo port install gmp
```

## Building
After having installed all of the required dependencies, you can build QuasiFit. For this, run
```
wget http://github.com/SoapZA/QuasiFit/releases/download/v0.1/quasifit-0.1.tar.bz2
tar -xjf quasifit-0.1.tar.bz2
cd quasifit-0.1/
./configure
make -j3
```

For users wishing to do development or want to stay up-to-date with the latest development, you will need to clone the git tree. This method is not recommended for users just wishing to use QuasiFit, as it requires the complete Autotools toolchain.
```
git clone https://github.com/SoapZA/QuasiFit.git
cd QuasiFit/
./autogen.sh
./configure
make -j3
```

The resulting executable `quasifit` located in `src/` can then be run. You can also install QuasiFit into a directory of your choosing if you specify the directory to the configure script with
```
./configure --prefix=<DIR>
```
and then install after the initial `make` with `make install`.

## Options
The QuasiFit sampler has a number of options for controlling the Metropolis-Hastings algorithm and I/O. See `quasifit -h` for more information.
