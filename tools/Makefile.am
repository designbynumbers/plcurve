AM_CFLAGS = -g -Wall #-finline-functions -ftree-vectorize -Wall
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = processpdstor
noinst_PROGRAMS = test_pyembed

test_pyembed_SOURCES = test_pyembed.c
test_pyembed_LDADD = -lplCurve $(PYTHON_EXTRA_LIBS)
test_pyembed_LDFLAGS = -static -L${top_builddir}/src/.libs $(PYTHON_LDFLAGS)
test_pyembed_CFLAGS = -I${top_srcdir}/src -I${top_srcdir}/pysrc/libpl $(PYTHON_CPPFLAGS)

processpdstor_SOURCES = processpdstor.c
processpdstor_LDADD = -lplCurve $(PYTHON_EXTRA_LIBS)
processpdstor_LDFLAGS = -static $(PYTHON_LDFLAGS)
processpdstor_CFLAGS = -I${top_srcdir}/src -I${top_srcdir}/pysrc/libpl -L${top_builddir}/src/.libs $(PYTHON_CPPFLAGS)
