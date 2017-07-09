SUBDIRS = pdstor_classify diagram_to_img

AM_CFLAGS = -g -Wall #-finline-functions -ftree-vectorize -Wall
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = processpdstor orientedshadows shadow2diagram
noinst_PROGRAMS = test_pyembed test_random_memleak

test_pyembed_SOURCES = test_pyembed.c
test_pyembed_LDADD = -lplCurve $(PYTHON_EXTRA_LIBS)
test_pyembed_LDFLAGS = -static -L${top_builddir}/src/.libs $(PYTHON_LDFLAGS)
test_pyembed_CFLAGS = -I${top_srcdir}/src -I${top_srcdir}/pysrc/libpl $(PYTHON_CPPFLAGS)

test_random_memleak_SOURCES = test_random_memleak.c
test_random_memleak_LDADD = -lplCurve $(PYTHON_EXTRA_LIBS)
test_random_memleak_LDFLAGS = -static -L${top_builddir}/src/.libs $(PYTHON_LDFLAGS)
test_random_memleak_CFLAGS = -I${top_srcdir}/src -I${top_srcdir}/pysrc/libpl $(PYTHON_CPPFLAGS)

processpdstor_SOURCES = processpdstor.c
processpdstor_LDADD = -lplCurve $(PYTHON_EXTRA_LIBS)
processpdstor_LDFLAGS = -static $(PYTHON_LDFLAGS)
processpdstor_CFLAGS = -I${top_srcdir}/src -I${top_srcdir}/pysrc/libpl -L${top_builddir}/src/.libs $(PYTHON_CPPFLAGS)

orientedshadows_SOURCES = orientedshadows.c
orientedshadows_LDADD = -lplCurve $(PYTHON_EXTRA_LIBS)
orientedshadows_LDFLAGS = -static $(PYTHON_LDFLAGS)
orientedshadows_CFLAGS = -I${top_srcdir}/src -I${top_srcdir}/pysrc/libpl -L${top_builddir}/src/.libs $(PYTHON_CPPFLAGS)

shadow2diagram_SOURCES = shadow2diagram.c
shadow2diagram_LDADD = -lplCurve $(PYTHON_EXTRA_LIBS)
shadow2diagram_LDFLAGS = -static $(PYTHON_LDFLAGS)
shadow2diagram_CFLAGS = -I${top_srcdir}/src -I${top_srcdir}/pysrc/libpl -L${top_builddir}/src/.libs $(PYTHON_CPPFLAGS)
