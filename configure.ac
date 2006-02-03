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
AC_INIT(libplCurve, 1.0, cantarel@math.uga.edu)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([config.h])
AM_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_AWK
AC_PROG_YACC
AC_PROG_CC
AC_PROG_CPP
AM_PROG_LEX
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB(m,sqrt)
AC_CHECK_LIB(argtable2,arg_parse,[builddemos=true],[builddemos=false; AC_MSG_WARN([[Demo programs require the argtable 2.0 library.]])])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h math.h malloc.h stdarg.h ctype.h stdio.h assert.h])

# Signal to AUTOMAKE if we have argtable2

AM_CONDITIONAL([BUILD_DEMOS],[test "$builddemos" = true])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_CHECK_FUNCS([floor sqrt ceil srandomdev sranddev random srandom])

# Can we inline and/or const?
AC_C_INLINE
AC_C_CONST

# Output
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
