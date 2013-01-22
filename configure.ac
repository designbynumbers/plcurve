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
AC_INIT(libplCurve, 6.0.0, cantarel@math.uga.edu)
#AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE

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
AC_SUBST(LIB_VERSION, 6:0:0)
AC_CONFIG_SRCDIR([configure])
AM_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])

AC_GNU_SOURCE
LT_INIT 

# Checks for programs.
AC_PROG_AWK
AC_PROG_YACC
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_F77
AC_PROG_CPP
AM_PROG_LEX
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL

# Checks for libraries.
AC_SEARCH_LIBS([sqrt],[m])
AC_SEARCH_LIBS([arg_parse],[argtable2],[have_argtable2=true],[have_argtable2=false])
AC_SEARCH_LIBS([gsl_ran_gaussian],[gsl], ,AC_MSG_ERROR([plCurve 6.x.x requires gsl (for random number generation)]),[-lgslcblas])

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

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h math.h malloc.h stdarg.h ctype.h stdio.h assert.h \
                  string.h stdbool.h float.h wait.h sys/wait.h unistd.h \
                  time.h sys/times.h sys/types.h fcntl.h signal.h complex.h gsl/gsl_rng.h gsl/gsl_randist.h gsl/gsl_combination.h])

# Signal to AUTOMAKE if we have argtable2

AM_CONDITIONAL([HAVE_ARGTABLE2],[test "x${have_argtable2}" = xtrue])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_CHECK_FUNCS([strcasestr mkstemp strcasecmp time])

# Can we inline and/or const?
AC_C_INLINE
AC_C_CONST

# Output
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

echo \
" ----------------------------------------------------------------------

  ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

  Installing to prefix: ${prefix}

  Compiler: '${CC} ${CFLAGS} ${CPPFLAGS}'
  
  Package Features:
      Found argtable2 to install demo programs: ${have_argtable2}

  Now type 'make @<:@<target>@:>@'

    where the optional <target> is:

    all					- build all binaries
    install				- install everything
    check 				- run self-tests


  -----------------------------------------------------------------------"
