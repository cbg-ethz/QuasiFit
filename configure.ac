dnl Initialise Autoconf
AC_PREREQ([2.69])
AC_INIT(
	[QuasiFit],
	[0.1],
	[david.seifert@bsse.ethz.ch],
	[quasifit],
	[http://www.cbg.ethz.ch/software/quasifit])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([src/quasifit.cpp])
AC_CONFIG_MACRO_DIR([m4])
AC_LANG([C++])



dnl ========================================================================
dnl Check whether we want to set defaults for CXXFLAGS, CPPFLAGS and LDFLAGS
dnl ========================================================================
AC_MSG_CHECKING([whether configure should try to set CXXFLAGS])
AS_IF(
	[test "x${CXXFLAGS+set}" = "xset"],
		[enable_cxxflags_setting=no],
		[enable_cxxflags_setting=yes]
)
AC_MSG_RESULT([$enable_cxxflags_setting])

AC_MSG_CHECKING([whether configure should try to set CPPFLAGS])
AS_IF(
	[test "x${CPPFLAGS+set}" = "xset"],
		[enable_cppflags_setting=no],
		[enable_cppflags_setting=yes]
)
AC_MSG_RESULT([$enable_cppflags_setting])

AC_MSG_CHECKING([whether configure should try to set LDFLAGS])
AS_IF(
	[test "x${LDFLAGS+set}" = "xset"],
		[enable_ldflags_setting=no],
		[enable_ldflags_setting=yes]
)
AC_MSG_RESULT([$enable_ldflags_setting])

AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug], [Disables setting of NDEBUG, which affects debugging asserts in Eigen])]
)




dnl =======================================
dnl Check for standard headers and programs
dnl =======================================
AC_PROG_SED
AC_PROG_CXX
AC_CHECK_HEADERS([cstdint])



dnl ==================================
dnl Set CXXFLAGS, CPPFLAGS and LDFLAGS
dnl ==================================
AS_IF(
	[test "x${enable_cxxflags_setting}" = "xyes"],
		[
		CXXFLAGS=""
		AX_APPEND_COMPILE_FLAGS([-O3], [CXXFLAGS])
		CXXFLAGS=$( echo $CXXFLAGS | $SED -e 's/^ *//' -e 's/ *$//' )
		]
)

AS_IF(
	[test "x${enable_cppflags_setting}" = "xyes"],
		[
			CPPFLAGS=""
			AS_IF(
				[test "x$enable_debug" != "xyes"],
					[
					AX_APPEND_COMPILE_FLAGS([-DNDEBUG], [CPPFLAGS])
					]
			)
			CPPFLAGS=$( echo $CPPFLAGS | $SED -e 's/^ *//' -e 's/ *$//' )
		]
)

AS_IF(
	[test "x${enable_ldflags_setting}" = "xyes"],
		[
		LDFLAGS=""
		AX_APPEND_LINK_FLAGS([-Wl,-O1 -Wl,--as-needed], [LDFLAGS])
		LDFLAGS=$( echo $LDFLAGS | $SED -e 's/^ *//' -e 's/ *$//' )
		]
)



dnl ===============
dnl Check for C++11
dnl ===============
AX_CXX_COMPILE_STDCXX_11([], [mandatory])



dnl ===================
dnl Initialise Automake
dnl ===================
AM_INIT_AUTOMAKE([1.14 foreign dist-bzip2 no-dist-gzip])
AM_MAINTAINER_MODE([enable])
LT_INIT



dnl ========================
dnl Check for precision type
dnl ========================
have_quad=no
have_truemulti=no
have_longdouble=no
AC_ARG_WITH([multi-precision],
	[AS_HELP_STRING(
		[--with-multi-precision],
		[Enable multi-precision floating-point arithmetic. This greatly increases the accuracy beyond standard long double floating-point arithmetic, at a large performance penalty.])]dnl
	m4_newline()dnl
		m4_text_wrap([enables GCC's libquadmath and uses __float128 as the extended floating-point type. This type is about 5-10x slower than using a native floating-point type.], [                          ], [--with-multi-precision=quad], [79])dnl
	m4_newline()dnl
		m4_text_wrap([enables GMP's floating-point types. These are the fastest known arbitrary precision floating-point types. Depending on the exact precision, you should expect a drop in performance by a factor of at least 100x compared to standard long double.], [                          ], [--with-multi-precision=truemulti], [79]),
	[dnl --with-multi-precision was provided
		AS_IF(
		[test "x$withval" = "xquad"],
			[dnl Quad precision requested
				AC_LANG_PUSH([C])
				save_LIBS="$LIBS"
				LIBS=
				AC_SEARCH_LIBS([quadmath_snprintf], [quadmath], [], [AC_MSG_ERROR([*** libquadmath not found])])
				PRECISION_LIBS="$LIBS"
				AC_SUBST(PRECISION_LIBS)
				LIBS="$save_LIBS"
				AC_LANG_POP([C])
			
				AC_DEFINE([HAVE_QUAD_PRECISION], [1], [Define if quad precision machinery is available])
				have_quad=yes
			],
		[test "x$withval" = "xtruemulti"],
			[dnl GMP arbitrary precision requested
				AC_LANG_PUSH([C])
				save_LIBS="$LIBS"
				LIBS=
				AC_SEARCH_LIBS([__gmpz_init], [gmp], [], [AC_MSG_ERROR([*** GMP library not found])])
				PRECISION_LIBS="$LIBS"
				AC_SUBST(PRECISION_LIBS)
				LIBS="$save_LIBS"
				AC_LANG_POP([C])
			
				AC_DEFINE([HAVE_TRUEMULTI_PRECISION], [1], [Define if arbitrary precision machinery is available])
				have_truemulti=yes
			],
		[dnl Unrecognised option
			AC_MSG_ERROR([--with-multi-precision=$withval is not a recognized argument])
		]
		)
	],
	[dnl --with-multi-precision was NOT provided		
		AC_DEFINE([HAVE_LONG_DOUBLE_PRECISION], [1], [Define if standard long double machinery is used])
		have_longdouble=yes
	]
)



dnl ================================
dnl Check for GNU Scientific Library
dnl ================================
PKG_CHECK_MODULES([GSL], [gsl])



dnl ==================
dnl Check for Eigen v3
dnl ==================
PKG_CHECK_MODULES([EIGEN], [eigen3])



dnl ===============
dnl Check for Boost
dnl ===============
BOOST_REQUIRE([1.50])
BOOST_SYSTEM
BOOST_CHRONO
BOOST_THREADS



dnl ========
dnl Finalise
dnl ========
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT



dnl ============================
dnl Report configuration to user
dnl ============================
AC_MSG_RESULT([
	$PACKAGE_NAME $VERSION

	CXX:        $CXX
	CXXFLAGS:   $CXXFLAGS
	CPPFLAGS:   $CPPFLAGS
	LDFLAGS:    $LDFLAGS
])

AS_IF(
	[test "x$have_quad" = "xyes"],
		[AC_MSG_RESULT([	Employing quad-precision floating point arithmetic (__float128, libquadmath)])],
	[test "x$have_truemulti" = "xyes"],
		[AC_MSG_RESULT([	Employing arbitrary-precision floating point arithmetic (mpf_t, GMP)])],
	[test "x$have_longdouble" = "xyes"],
		[AC_MSG_RESULT([	Employing standard long-double floating point arithmetic (long double, 80-bit on x86)])],
	[]
)
