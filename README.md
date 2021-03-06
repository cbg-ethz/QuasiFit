# QuasiFit 0.3
David Seifert (david.seifert@bsse.ethz.ch)
Niko Beerenwinkel (niko.beerenwinkel@bsse.ethz.ch)

## Citation
If you find QuasiFit useful, please cite our paper in Genetics

   David Seifert, Francesca Di Giallonardo, Karin J. Metzner, Huldrych F. Günthard, and Niko Beerenwinkel.
   **A Framework for Inferring Fitness Landscapes of Patient-Derived Viruses Using Quasispecies Theory.**
   _Genetics_ 2015, 199(1).
   DOI: 10.1534/genetics.114.172312

## Introduction
QuasiFit is an MCMC sampler that implements (relative) fitness inference for NGS data assuming a mutation-selection equilibrium. From the posterior, neutral networks and epistasis can be determined.

## Binaries
We have pre-compiled 64-bit binaries for Linux and Mac users:

* **Linux**: `quasifit-linux-static-amd64`

  You will require a 64-bit distribution having _at least_ glibc 2.3.2. Any distribution from the past 10 years should work.  
  The linux binary was built on Debian Etch 4.0r9 64-bit.

* **Mac**: `quasifit-mac-static-amd64`

  You will require at least Mac OS X 10.4.11 running on a 64-bit Mac.

QuasiFit was built on both platforms with GSL 1.16 and Boost 1.55 with GCC 4.8.2 on -O2 optimizations. The main code (including Eigen) was compiled with -O3 optimizations. All libraries, including the C++ runtime libraries, have been linked statically to produce a binary that has **no external dependencies**, in other words, _they are directly useable_.

All static binaries can be downloaded from the main git tree for the most recent release or from the releases page.

## Prerequisites
If you wish to compile QuasiFit from source, you will require the following components (these are **not** necessary for running the statically linked binary):

1. **GSL**; _somewhat recent_ release (http://www.gnu.org/software/gsl/)

   The GNU Scientific Library is required for random number generating functions.

2. **Eigen**; _at least 3.2_ (http://eigen.tuxfamily.org/)

   Eigen forms the core mathematics library of QuasiFit, with all its linear algebra routines.

3. **Boost**; _at least 1.50_ (http://www.boost.org/)

   Boost provides the necessary abstraction for time routines and thread handling. Also abstracts the different precision types.

4. **GMP** (optional); _somewhat recent_ release (http://www.gmplib.org/)

   The GNU Multiple Precision Arithmetic Library (GMP) provides the basis (mpf_t) for arbitrary precision floating-point calculations. It is only required if you wish to build an arbitrary precision sampler.

5. **libquadmath** (optional); _at least GCC 4.6_ (http://gcc.gnu.org/onlinedocs/libquadmath/)

   GCC's libquadmath provides the __float128 quad-precision floating point type and associated operations. This is an internal GCC library that is included with GCC since 4.6. It is only required if you wish to use a quad-precision sampler. Quad precision represents a trade-off between performance and precision.

Furthermore, you will require a compiler that can handle **C++0x** (which includes all C++11 compilers). QuasiFit has been successfully compiled with GCC 4.4 on RHEL 6, GCC 4.8 on Gentoo/Debian Etch and icc 12.0 on RHEL 6. **Please note** that building an arbitrary precision sampler requires either GCC or Clang, as the Intel C++ Compiler has a known bug in handling Boost's Multiprecision library.

If you wish to do development, you will require parts of the extended GNU toolchain (the infamous Autotools):

1. **Autoconf**; latest 2.69 release (http://www.gnu.org/software/autoconf/)

   GNU Autoconf produces the ./configure script from configure.ac.

2. **Automake**; latest 1.14 release (http://www.gnu.org/software/automake/)

   GNU Automake produces the Makefile.in precursor, that is processed with ./configure to yield the final Makefile.

3. **Libtool**; latest 2.4.2 release (http://www.gnu.org/software/libtool/)

   GNU Libtool is required as a dependency of boost.m4.

QuasiFit is strongly intertwined with libraries and programs that heavily rely on features of UNIX-like systems, hence supporting Microsoft Windows is not a goal (in particular, building the GNU Scientific Library and using the GNU build system on Windows is a nightmare).

## Preparing
To install the aforementioned dependencies, follow the guides here.

### Linux
Due to the large inherent heterogeneity of the Linux landscape, we will only detail the procedure of installing dependencies for Ubuntu 14.04 LTS here. The procedure should be very similiar for Debian.
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

1. **Xcode** for your platform (4.2 for 10.6; 4.6.3 for 10.7; 5.1.x for 10.8 & 10.9) either via the Mac App Store or by downloading the disk image from the Apple Developer Connection (http://developer.apple.com/).

2. **Command Line Tools** for Xcode. Since Xcode 4.3 Apple has stopped shipping command line tools with the standard Xcode package. You will need to install these via "Downloads" in the "Xcode" -> "Preferences" menu, or (preferably) by downloading the latest appropriate "Command Line Tools" package from the Apple Developer Connection.

3. **MacPorts** (http://www.macports.org/install.php).

Henceforth, we assume Xcode, the Command Line Tools and Macports to be installed.

Install the remaining libraries from MacPorts by performing
```
# install general prerequisites
sudo port install wget pkgconfig

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

If you wish to work with the bleeding-edge release of QuasiFit, you will need the complete GNU Autotools toolchain. It should be reiterated here that the recommended way of building QuasiFit is by downloading the provided tarball and using either the included static binaries or compiling from source. The Git tree needs to be bootstrapped to produce the various scripts. To install the Autotools:
```
# install the GNU toolchain
sudo port install autoconf automake libtool
```

## Building
After having installed all of the required dependencies, you can build QuasiFit. For this, run
```
wget --no-check-certificate http://github.com/SoapZA/QuasiFit/releases/download/v0.3/quasifit-0.3.tar.bz2 -O - | tar xj
cd quasifit-0.3/
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


## Usage
### Input file
QuasiFit is versatile in what sequence file formats it can take as input:

1. **Generic FASTA file**. One drawback of the Generic FASTA input file is that it cannot include unobserved haplotypes, as every sequence in a FASTA file represents one observation. To include unobserved haplotypes into the inference procedure, use one of the two other input formats.
2. **QuasiRecomb FASTA file**. QuasiFit was designed with QuasiRecomb's output file structure as input. QuasiRecomb writes the statistics of its output into the sequence identifier field of a FASTA file:
   
   ```
   >read0_0.5026
   TAGAAGATATGGAGTTGCCAGGGAGGTGGA
   ```
   
   QuasiFit parses this expression and extracts the counts. While QuasiRecomb by itself does not include any unobserved haplotypes, you can insert these manually into the output FASTA file.
3. **QuasiFit input file**. QuasiFit can load its own kind of input file, which consists simply of comma-separated values. Every line includes one haplotype, separated by its observed count with a comma. For instance
   
   ```
   AAA,1
   AAT,1
   ```
   
   would be a two-haplotype input file for QuasiFit. This is the same input file as used in the supplemental information of the main paper in the section "Unobserved haplotypes simulations".
   
Obviously, all haplotypes have to be of the same length for the quasispecies model (in its simplest form) to be applicable.

### Output files
QuasiFit produces multiple output files:

1. `<FILE>-f.csv`: these contain the actual fitness samples from the fitness manifold. Be aware that QuasiFit will write out the full number of decimal digits for each floating-point value, hence this file can become somewhat large.
2. `<FILE>-p.csv`: these contain the estimated population distribution samples. Every row should theoretically sum to 1 (within numerical truncation errors), as every component represents the probability of a haplotype in an asymptotically infinite population.
3. `<FILE>-r.csv`: these contain the samples from the subset of the Euclidean space, which in fact is the true sampling space. Every row will include at least one 0, as the dimensionality of the euclidean space is the same as the degrees of freedom of the quasispecies distribution simplex, namely #Haplotypes - 1.
4. `<FILE>-diag.csv`: these contain 3 columns of diagnostic data. The first column represents the logarithm of the posterior density function (up to a constant shift), the second column represents the logarithm of the absolute value of the determinant of the Jacobian of h(p), and finally, the third column represents the logarithm of the multinomial likelihood (excluding the constant prefactor).

All of these files can be analysed with standard tools. We recommend using R for its sophisticated plotting capabilities. To load one such file, fire up R and use for instance
```
diagnostic_data = read.table("<FILE>-diag.csv", header=TRUE, sep=",", colClasses="numeric")
```
to have a look at the diagnostic data. QuasiFit does not automatically detect or remove the burn-in phase, as this is generally tricky and should be left for the user to determine. To determine the burn-in phase, just plot the beginning of the log Posterior and notice when the values flatten out and converge to their supposed stationary distribution. To plot the first 25'000 values of the log Posterior, proceed with
```
plot(diagnostic_data$LogPost[1:25000], xlab="MCMC iteration", ylab="Log Posterior", type="l")
```
to get something like this
<p align="center">
  <img src="misc/BurninPlot.png?raw=true"/>
</p>

Notice how the MCMC chains converge to the stationary distribution at around 10'000 iterations - this would be considered the burn-in phase. The drop in the log Posterior from the initial value is a result of starting at the MLE of the problem and the general curse of dimensionality. For more complicated post-MCMC diagnostics, try the `coda` package from CRAN (http://cran.r-project.org/web/packages/coda/index.html). For instructions on verifying convergence, see http://www.people.fas.harvard.edu/~plam/teaching/methods/convergence/convergence_print.pdf.

### Unconnected haplotypes
QuasiFit makes strong assumptions on the connectedness of haplotypes. For instance, the haplotypes **AAA** and **TTT** are separated by a mutational step that requires 3 simultaneous mutations in one replication cycle. Mainly for numerical reasons, this causes the mutation matrix **Q** to become numerically reducible and a global equilibrium distribution of the quasispecies equation is not numerically guaranteed anymore.

To circumvent this issue, we have detailed a procedure in the main paper in the section "Haplotype space and mutation probabilities" that inserts a minimal number of unobserved haplotypes such that we arrive at a network of haplotypes, where every haplotype can mutate into every other haplotype by taking only simultaneous k mutations per replication cycle (in practice we require k = 1). To do this, we have included our MATLAB script `curateSample.m` in the `scripts/` folder. Fire up a MATLAB session and run for instance

```
curateSample('quasispecies.fasta')
```

where `quasispecies.fasta` is the output file of QuasiRecomb. The curateSample script converts QuasiRecomb's output to a QuasiFit input file and includes the minimal number of unobserved haplotypes to make the haplotype graph connected with one component for given k.