AM_CFLAGS = -g -finline-functions -Wall

include_HEADERS = plCurve.h
lib_LTLIBRARIES = libplCurve.la libplCurve_safe.la
noinst_LTLIBRARIES = libplCurve_cov.la

check_PROGRAMS = run_tests exit_failure_tests
bin_PROGRAMS = knot_diagram
noinst_PROGRAMS = spline_view gcov_tests

TESTS = run_tests

EXTRA_DIST = bad_vects

LDADD = libplCurve.la

libplCurve_la_SOURCES = plCurve.c plCurve.h vector.c splines.c
libplCurve_la_CFLAGS = -DNDEBUG -O2 $(AM_CFLAGS)
libplCurve_la_LDFLAGS = -version-info $(LIB_VERSION)

libplCurve_safe_la_SOURCES = plCurve.c plCurve.h vector.c splines.c
libplCurve_safe_la_LDFLAGS = -version-info $(LIB_VERSION)

libplCurve_cov_la_SOURCES = plCurve.c plCurve.h vector.c splines.c
libplCurve_cov_la_CFLAGS = -g -Wall -fprofile-arcs -ftest-coverage
libplCurve_cov_la_LDFLAGS = -static -fprofile-arcs

run_tests_SOURCES = run_tests.c libplCurve_safe.la
run_tests_LDADD = libplCurve_safe.la

gcov_tests_SOURCES = run_tests.c libplCurve_cov.la
gcov_tests_LDADD = libplCurve_cov.la
gcov_tests_LDFLAGS = -static -fprofile-arcs
gcov_tests_CFLAGS = -g -Wall -fprofile-arcs -ftest-coverage

exit_failure_tests_SOURCES = exit_failure_tests.c libplCurve_cov.la
exit_failure_tests_LDADD = libplCurve_cov.la
exit_failure_tests_LDFLAGS = -static -fprofile-arcs
exit_failure_tests_CFLAGS = -g -Wall -fprofile-arcs -ftest-coverage

knot_diagram_SOURCES = knot_diagram.c libplCurve_safe.la
knot_diagram_LDADD = libplCurve_safe.la
if HAVE_ARGTABLE2
knot_diagram_LDFLAGS = -static -largtable2
else
knot_diagram_LDFLAGS = -static
endif
