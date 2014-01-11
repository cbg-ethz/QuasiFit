AC_INIT([QuasiFit], [0.1], [david.seifert@bsse.ethz.ch], [quasifit], [http://www.cbg.ethz.ch/software/quasifit])
: ${CXXFLAGS="-O3 -fno-rtti"}
AC_PREREQ([2.69])

AM_INIT_AUTOMAKE([1.14 foreign dist-xz])
AM_MAINTAINER_MODE([enable])
AC_PROG_LIBTOOL
AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([], [mandatory])

dnl Check for GNU Scientific Library
PKG_CHECK_MODULES([GSL], [gsl])

dnl Check for Eigen v3
PKG_CHECK_MODULES([Eigen], [eigen3])

dnl Check for Boost
dnl m4_include([m4/boost.m4])
BOOST_REQUIRE([1.50])
BOOST_SYSTEM
BOOST_CHRONO
BOOST_THREADS

AC_CHECK_FUNCS([pow])
AC_CHECK_FUNCS([strtoull])
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_FUNC_STRTOD
AC_TYPE_INT64_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

AC_OUTPUT([Makefile])
