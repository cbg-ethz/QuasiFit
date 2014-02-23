dnl Initialise Autoconf
AC_PREREQ([2.69])
AC_INIT([QuasiFit],
	[0.1],
	[david.seifert@bsse.ethz.ch],
	[quasifit],
	[http://www.cbg.ethz.ch/software/quasifit])
: ${CXXFLAGS="-O3"}
: ${LDFLAGS="-Wl,-O1 -Wl,--as-needed"}
AC_CONFIG_HEADERS([config.h])
AC_LANG([C++])

dnl Initialise Automake
AM_INIT_AUTOMAKE([1.14 foreign dist-bzip2])
AM_MAINTAINER_MODE([enable])
AC_PROG_LIBTOOL
AC_CONFIG_MACRO_DIR([m4])

dnl Standard checks
AC_CHECK_FUNCS([pow])
AC_CHECK_FUNCS([strtoull])
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_FUNC_STRTOD
AC_TYPE_INT64_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([], [mandatory])

dnl ./configure options
dnl Quad Precision
AC_ARG_ENABLE([quad-precision],
    AS_HELP_STRING([--enable-quad-precision], [Enable quad precision floating-point arithmetic. This greatly increases the accuracy of the estimates, as floating-point truncation errors are avoided. For the time being, this requires the Intel C++ compiler, due to performance considerations.]))

dnl Multi Precision
AC_ARG_ENABLE([multi-precision],
    AS_HELP_STRING([--enable-multi-precision], [Enable arbitrary precision floating-point arithmetic with correct rounding. This greatly increases the accuracy beyond quad precision floating-point arithmetic, incurring a large performance penalty. Currently, arbitrary precision floating-point arithmetic is implemented with libmpfr, Pavel Holoborodko's MPFR C++ and the MPreal headers in Eigen.]))

dnl Check whether both precision modes are enabled at once
AS_IF([test "x$enable_quad_precision" = "xyes" && test "x$enable_multi_precision" = "xyes"],
	[AC_MSG_ERROR([You cannot simultaneously have --enable-quad-precision and --enable-multi-precision!])])

PRECISION_CXXFLAGS=""
dnl Perform additional checks for extended precision modes
AS_IF([test "x$enable_quad_precision" = "xyes"], [
	dnl Check for required Intel C++ semantics for quad precision
	AX_CHECK_COMPILE_FLAG([-Qoption,cpp,--extended_float_type], [], [AC_MSG_ERROR([The required flag -Qoption,cpp,--extended_float_type is not recognized. This is probably because you aren't using the Intel C++ compiler or have omitted to specify if with CXX=])])
	PRECISION_CXXFLAGS="-Qoption,cpp,--extended_float_type"
	
	AC_DEFINE([HAVE_QUAD_PRECISION], [1], [Define if Intel's Quad precision machinery is available])
])

AS_IF([test "x$enable_multi_precision" = "xyes"], [
	dnl Check for required MPFR semantics
	
	AC_DEFINE([HAVE_MULTI_PRECISION], [1], [Define if MPFR's arbitrary precision machinery is available])
])
AC_SUBST([PRECISION_CXXFLAGS])

dnl Check for GNU Scientific Library
PKG_CHECK_MODULES([GSL], [gsl])

dnl Check for Eigen v3
PKG_CHECK_MODULES([Eigen], [eigen3])

dnl Check for Boost
BOOST_REQUIRE([1.50])
BOOST_SYSTEM
BOOST_CHRONO
BOOST_THREADS

AC_CONFIG_FILES([Makefile])
AC_OUTPUT