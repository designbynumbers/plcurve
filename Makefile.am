AM_CFLAGS = -g -finline-functions -ftree-vectorize -Wall
# -I$(top_srcdir)/argtable2-13/src -I$(top_srcdir)/../ridgerunner/argtable2-13/src  -I$(top_srcdir)/../argtable2-13/src
# AM_LDFLAGS = -L$(top_srcdir)/argtable2-13/src -L$(top_srcdir)/../argtable2-13/src -largtable2
ACLOCAL_AMFLAGS = -I m4

include_HEADERS = plCurve.h nplCurve.h
lib_LTLIBRARIES = libplCurve.la libplCurve_safe.la libnplCurve.la libnplCurve_safe.la libplCurveStatic.la
#noinst_LTLIBRARIES = libplCurve_cov.la

lib_LIBRARIES = libplCurve.a libplCurveStatic.a

libplCurve_a_SOURCES = plCurve.c plCurve.h vector.c splines.c pllmpoly02.c plcTopology.c \
		       plwhitten.c nearest_neighbor.c symmetries.c matrix.c matrix.h plcRandomPolygon.c plc_dihedrals.c
libplCurve_a_CFLAGS = $(AM_CFLAGS)

libplCurveStatic_a_SOURCES = plCurve.c plCurve.h vector.c splines.c pllmpoly02.c plcTopology.c \
		       plwhitten.c nearest_neighbor.c symmetries.c matrix.c matrix.h plcRandomPolygon.c plc_dihedrals.c
libplCurveStatic_a_CFLAGS = $(AM_CFLAGS)

check_PROGRAMS =  nn_test symmetry_test randompolygon_test eq_primitives_test loopclosure_test

if HAVE_ARGTABLE2
bin_PROGRAMS = knot_diagram resample plcurvature plcompEric
else
bin_PROGRAMS = knot_diagram plcompEric
endif

noinst_PROGRAMS = spline_view # gcov_tests
noinst_HEADERS = homfly.h

TESTS = nn_test symmetry_test randompolygon_test eq_primitives_test loopclosure_test

dist_data_DATA = plCurve.typemap

EXTRA_DIST = bad_vects/bad_1.vect bad_vects/bad_2.vect bad_vects/bad_3.vect \
	     bad_vects/bad_4.vect bad_vects/bad_5.vect bad_vects/bad_6.vect \
	     bad_vects/bad_7.vect bad_vects/bad_8.vect bad_vects/bad_9.vect \
	     bad_vects/bad_10.vect bad_vects/bad_11.vect bad_vects/bad_12.vect \
	     bad_vects/bad_13.vect bad_vects/bad_14.vect \
	     data/circle_closure_0.vect \
             data/circle_closure_1.vect \
             data/circle_closure_2.vect \
             data/circle_closure_3.vect \
             data/circle_closure_4.vect \
             data/circle_closure_5.vect \
             data/circle_closure_6.vect \
             data/circle_closure_7.vect \
             data/circle_closure_8.vect \
             data/circle_closure_9.vect \
             data/random_closure_0.vect \
             data/random_closure_1.vect \
             data/random_closure_2.vect \
             data/random_closure_3.vect \
             data/random_closure_4.vect \
             data/random_closure_5.vect \
             data/random_closure_6.vect \
             data/random_closure_7.vect \
             data/random_closure_8.vect \
             data/random_closure_9.vect

LDADD = libplCurve.la

libplCurveStatic_la_SOURCES = plCurve.c plCurve.h vector.c splines.c pllmpoly02.c plcTopology.c \
				plwhitten.c nearest_neighbor.c symmetries.c matrix.c matrix.h plcRandomPolygon.c plc_dihedrals.c
libplCurveStatic_la_CFLAGS = -DNDEBUG  $(AM_CFLAGS)
libplCurveStatic_la_LDFLAGS = -version-info $(LIB_VERSION) -static

libplCurve_la_SOURCES = plCurve.c plCurve.h vector.c splines.c pllmpoly02.c plcTopology.c \
			plwhitten.c nearest_neighbor.c symmetries.c matrix.c matrix.h plcRandomPolygon.c plc_dihedrals.c
libplCurve_la_CFLAGS = -DNDEBUG  $(AM_CFLAGS)
libplCurve_la_LDFLAGS = -version-info $(LIB_VERSION)

libplCurve_safe_la_SOURCES = plCurve.c plCurve.h vector.c splines.c pllmpoly02.c plcTopology.c \
				plwhitten.c nearest_neighbor.c symmetries.c matrix.c matrix.h plcRandomPolygon.c plc_dihedrals.c
libplCurve_safe_la_LDFLAGS = -version-info $(LIB_VERSION)

libnplCurve_la_SOURCES = nplCurve.c nplCurve.h nvector.c
libnplCurve_la_CFLAGS = -DNDEBUG  $(AM_CFLAGS)
libnplCurve_la_LDFLAGS = -version-info $(LIB_VERSION)

libnplCurve_safe_la_SOURCES = nplCurve.c nplCurve.h nvector.c
libnplCurve_safe_la_LDFLAGS = -version-info $(LIB_VERSION)

#libplCurve_cov_la_SOURCES = plCurve.c plCurve.h vector.c splines.c pllmpoly02.c plcTopology.c \
			plwhitten.c nearest_neighbor.c symmetries.c matrix.c matrix.h
#libplCurve_cov_la_CFLAGS = -g -Wall -fprofile-arcs -ftest-coverage
#libplCurve_cov_la_LDFLAGS = -static -fprofile-arcs

#run_tests_SOURCES = run_tests.c libplCurve.a
#run_tests_LDADD = -lplCurve

nn_test_SOURCES = nn_test.c libplCurve.a
nn_test_LDADD = -lplCurve

symmetry_test_SOURCES = symmetry_test.c libplCurve.a
symmetry_test_LDADD = -lplCurve

eq_primitives_test_SOURCES = eq_primitives_test.c libplCurve.a
eq_primitives_test_LDADD = -lplCurve

loopclosure_test_SOURCES = loopclosure_test.c libplCurve.a
loopclosure_test_LDADD = -lplCurve

#gcov_tests_SOURCES = run_tests.c libplCurve_cov.la
#gcov_tests_LDADD = libplCurve_cov.la
#gcov_tests_LDFLAGS = -static -fprofile-arcs
#gcov_tests_CFLAGS = -g -Wall -fprofile-arcs -ftest-coverage

#exit_failure_tests_SOURCES = exit_failure_tests.c libplCurve.a
#exit_failure_tests_LDADD = -lplCurve
#exit_failure_tests_LDFLAGS = -static
#exit_failure_tests_CFLAGS = -g -Wall

knot_diagram_SOURCES = knot_diagram.c libplCurve.a
knot_diagram_LDADD = -lplCurve
if HAVE_ARGTABLE2
knot_diagram_LDFLAGS = -static
else
knot_diagram_LDFLAGS = -static
endif

resample_SOURCES = resample.c libplCurve.a
resample_LDADD = -lplCurve
resample_LDFLAGS = -static

plcurvature_SOURCES = plcurvature.c libplCurve.a
plcurvature_LDADD = -lplCurve
plcurvature_LDFLAGS = -static

if HAVE_F77
AM_FFLAGS = -g
plcompEric_SOURCES = compEric.f
endif

if HAVE_SWIG
BUILT_SOURCES = $(srcdir)/plCurve_wrap.c
SWIG_SOURCES = plCurve.i

pkgpython_PYTHON = plcurve.py __init__.py
pkgpyexec_LTLIBRARIES = _plcurve.la
_plcurve_la_SOURCES = $(srcdir)/plCurve_wrap.c $(SWIG_SOURCES)
_plcurve_la_CFLAGS = $(AX_SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)
_plcurve_la_LDFLAGS = -module
_plcurve_la_LIBADD = libplCurve.la

$(srcdir)/plCurve_wrap.c : $(SWIG_SOURCES)
	$(SWIG) $(AX_SWIG_PYTHON_OPT) -builtin -I$(top_srcdir) -o $@ $<
endif

#clean-local:
#	rm *.gcno *.gcda
