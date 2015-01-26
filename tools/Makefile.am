AM_CFLAGS = -g -Wall #-finline-functions -ftree-vectorize -Wall
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = test_pyembed
noinst_PROGRAMS = test_pyembed

test_pyembed_SOURCES = test_pyembed.c
test_pyembed_LDADD = -lplCurve $(PYTHON_EXTRA_LIBS)
test_pyembed_LDFLAGS = -static $(PYTHON_LDFLAGS)
test_pyembed_CFLAGS = -I${top_builddir}/src -I${top_builddir}/pysrc/libpl $(PYTHON_CPPFLAGS)
