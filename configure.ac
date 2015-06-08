#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Autoconf 2.50 wants to disallow AM_ names.  We explicitly allow
# the ones we care about.
ifdef([m4_pattern_allow],
      [m4_pattern_allow([^AM_(CONFIG_HEADER|PROG_LEX)$])])dnl


AC_PREREQ(2.57)
AC_COPYRIGHT(
------------------------------- libplCurve.a ---------------------------------

                  Copyright (c) 2006 University of Georgia.

 Distributed under terms of the GNU General Public License version 2.1
                       or later (see COPYING.LIB-2.1).
)

# Don't forget to change the LIB_VERSION (and read the note if you are
# wondering what the numbers mean.
AC_INIT(libplCurve, 7.4.1, cantarel@math.uga.edu)
#AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE([])

#
# The version number below has a very specific meaning.  It is in the form
# a:b:c.  The first number, a, is the interface version number.  Any time we
# add or remove a function, or change the order of parameters or change the way
# that the library reads or writes files, this number needs to be incremented.
# The second number, b, is the implementation number *of this interface*.  When
# a changes, b goes back to zero.  This is the number that needs to change
# every time we put out a new library distribution.  The third number, c, is
# the interface backward-compatibility count.  Suppose that we simply add a
# needed function to the library.  The new interface still provides everything
# that the old version provided and in the same way.  So the
# backward-compatibility count gets incremented.  Of course it gets reset to 0
# when a new interface is not binary-compatible with the old one.
#
# Due to various broken ideas, we have to have two versions of this number, one
# with colons below and one with periods above.  Please keep them in sync!
#
AC_SUBST(LIB_VERSION, 7:4:1)
AC_CONFIG_SRCDIR([src/plCurve.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

#AC_GNU_SOURCE
LT_INIT

# Checks for programs.
#AC_PROG_AWK
#AC_PROG_YACC
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CPP

#AM_PROG_LEX
AC_PROG_LN_S

AC_PROG_MAKE_SET
AC_PROG_LIBTOOL

# Checks for libraries. The syntax here is
#
# AC_SEARCH_LIBS([a function in the library],[the name of the library])
#
# if you find it, it will be automatically added to paths and the link commands
# for the project.
AC_SEARCH_LIBS([sqrt],[m])
AC_SEARCH_LIBS([arg_parse],[argtable2],[have_argtable2=true],[have_argtable2=false])
AC_SEARCH_LIBS([cblas_dgemm],[gslcblas])
AC_SEARCH_LIBS([gsl_ran_gaussian],[gsl], ,AC_MSG_ERROR([plCurve 6.x.x and higher require gsl (for random number generation)]),[-lgslcblas])

if test "x${have_argtable2}" = xtrue; then
   AC_CHECK_HEADERS([argtable2.h], [], [have_argtable2=false])
fi

if test "x${have_argtable2}" = xfalse; then
   AC_MSG_WARN([
  -----------------------------------------------------------
  Unable to find argtable2 on this system. Some demo programs
  will not be installed in this build of plCurve.
  -----------------------------------------------------------])
fi

# Check for SWIG and Python if old wrappers are enabled.
AC_ARG_ENABLE([swig],
        [AS_HELP_STRING([--enable-swig],
        [build old swig Python bindings])],
        [enable_swig=$enableval],
        [enable_swig=no])

if test "x${enable_swig}" = xyes; then
   AX_PKG_SWIG(2.0.0,[have_swig=true],[AC_MSG_ERROR([SWIG interface is to be built (--enable-swig), but SWIG was not found.]); have_swig=false])
   if test "x${have_swig}" = xtrue; then
      AX_SWIG_PYTHON #I think this forces a dependency on python... TODO: make optional?
      AM_PATH_PYTHON(2.7)
      CIT_NUMPY_PYTHON_MODULE(1.8.0,[have_numpy=true],[have_numpy=false])
      if test "x${have_numpy}" = xtrue; then
         CIT_NUMPY_INCDIR
      else
        AC_MSG_WARN([An appropriate numpy installation was not found with Python. numpy functionality will not be included.])
      fi
   fi
fi

AM_CONDITIONAL([HAVE_SWIG],[test "x${have_swig}" = xtrue])

# Python extension library checks (Python, Cython, numpy)
AC_ARG_ENABLE([python],
        [AS_HELP_STRING([--disable-python],
                [do not build Python bindings])],
        [enable_python=$enableval],
        [enable_python=yes])

if test "x${enable_python}" = xyes; then
   AM_PATH_PYTHON(2.7)
   AX_PYTHON_DEVEL(>= '2.7')
   CYTHON_CHECK_VERSION([0.21.0],
        [have_cython=true],
        [AC_MSG_FAILURE([Python bindings were not disabled but Cython was not found.]);
         have_cython=false])
   CIT_NUMPY_PYTHON_MODULE(1.8.0,
        [have_numpy=true;
         CIT_NUMPY_INCDIR],
        [AC_MSG_WARN([Python bindings were not disabled but Numpy was not found.]);
         have_numpy=false]])
fi

AM_CONDITIONAL([BUILD_PYTHON],[test "x${enable_python}" = xyes])

# Checks for header files. Just add the name of the header file that you
# require inside the [ ] in AC_CHECK_HEADERS. If you have it, you'll have
# a preprocessor symbol set in the form HAVE_MYHEADER_H, so you use it with
#
# #ifdef HAVE_MYHEADER_H
#    #include <myheader.h>
# #endif
#
# in your source files. Paths are all handled for you, so you shouldn't need
# to specify more about where it is.

AC_HEADER_STDC
AC_HEADER_ASSERT
AC_CHECK_HEADERS([stdlib.h math.h stdarg.h ctype.h stdio.h assert.h \
                  string.h stdbool.h float.h wait.h sys/wait.h unistd.h \
                  time.h sys/times.h sys/types.h fcntl.h signal.h complex.h \
		  gsl/gsl_rng.h gsl/gsl_randist.h gsl/gsl_combination.h gsl/gsl_permutation.h \
		  libintl.h limits.h strings.h stdint.h ])


# Signal to AUTOMAKE if we have argtable2

AM_CONDITIONAL([HAVE_ARGTABLE2],[test "x${have_argtable2}" = xtrue])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_MODE_T
AC_TYPE_SIZE_T

# Checks for library functions.
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([strcasestr mkstemp strcasecmp time floor sqrt ceil srandomdev sranddev random srandom])
AC_CHECK_FUNCS([bzero memset mkdir regcomp strchr strcspn strerror strrchr strstr strtol getline fstat64 memcpy])

# Can we inline and/or const?
AC_C_INLINE
AC_C_CONST

# Output
AC_CONFIG_FILES([Makefile
		src/Makefile
                planarmap/Makefile
                pysrc/Makefile
                tools/Makefile
		test/Makefile])

AC_CONFIG_SUBDIRS([judy])
AC_OUTPUT

echo \
" ----------------------------------------------------------------------

  ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

  Installing to prefix: ${prefix}

  Compiler: '${CC} ${CFLAGS} ${CPPFLAGS}'

  Package Features:
      Found argtable2 to install demo programs: ${have_argtable2}
      Found SWIG to build python interface: ${have_swig}

  Now type 'make @<:@<target>@:>@'

    where the optional <target> is:

    all					- build all binaries
    install				- install everything
    check 				- run self-tests


  -----------------------------------------------------------------------"
