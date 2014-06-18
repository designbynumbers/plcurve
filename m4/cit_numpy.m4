# -*- Autoconf -*-


## --------------------------- ##
## Autoconf macros for Numpy. ##
## --------------------------- ##

# CIT_NUMPY_PYTHON_MODULE
# Determine whether the numpy Python module is available.
AC_DEFUN([CIT_NUMPY_PYTHON_MODULE], [
AC_REQUIRE([AM_PATH_PYTHON])
AC_MSG_CHECKING(for numpy python module)
$PYTHON -c "import numpy; import sys; maj,min,rev = tuple(int(d) for d in numpy.version.full_version.split('.')); sys.exit(0 if (maj >= 1 and min >= 8) else 1)" 2>/dev/null
if test $? == 0; then
  AC_MSG_RESULT(found)
  m4_ifval([$2],[$2],[:])
else
  AC_MSG_WARN(not found)
  AC_DEFINE([HAVE_NUMPY], [0], [Define to 1 if numpy header files found on the system])
  m4_ifval([$3],[$3],[:])
fi
]) dnl CIT_NUMPY_PYTHON_MODULE

# NUMPY_INCDIR
# -----------------
# Determine the directory containing <numpy/arrayobject.h>
AC_DEFUN([CIT_NUMPY_INCDIR], [
AC_REQUIRE([AM_PATH_PYTHON])
AC_CACHE_CHECK([for numpy include directory],
    [_cv_numpy_incdir],
    [_cv_numpy_incdir=`$PYTHON -c "import numpy; numpypath=numpy.__path__[[0]]; print '%s/core/include' % numpypath"`])
AC_SUBST([NUMPY_INCDIR], [$_cv_numpy_incdir])
AC_DEFINE([HAVE_NUMPY], [1], [Define to 1 if numpy header files found on the system])
])dnl CIT_NUMPY_INCDIR


dnl end of file
