AM_CFLAGS = -g -O2 -finline-functions -Wall

lib_LTLIBRARIES = libplCurve.la libplCurve_safe.la
noinst_LTLIBRARIES = libplCurve_cov.la

check_PROGRAMS = run_tests exit_failure_tests
noinst_PROGRAMS = knot_diagram

TESTS = run_tests

EXTRA_DIST = bad_vects

libplCurve_la_SOURCES = plCurve.c plCurve.h vector.c splines.c
libplCurve_la_CFLAGS = -DNDEBUG $(AM_CFLAGS)
libplCurve_la_LDFLAGS = -version-info $(LIB_VERSION)

libplCurve_safe_la_SOURCES = plCurve.c plCurve.h vector.c splines.c
libplCurve_safe_la_LDFLAGS = -version-info $(LIB_VERSION)

libplCurve_cov_la_SOURCES = plCurve.c plCurve.h vector.c splines.c
libplCurve_cov_la_CFLAGS = -g -Wall -fprofile-arcs -ftest-coverage
libplCurve_cov_la_LDFLAGS = -static -version-info -fprofile-arcs $(LIB_VERSION)

run_tests_SOURCE = run_tests.c libplCurve_cov.la
run_tests_LDADD = libplCurve_cov.la
run_tests_LDFLAGS = -static -fprofile-arcs
run_tests_CFLAGS = -g -Wall -fprofile-arcs -ftest-coverage

exit_failure_tests_SOURCE = exit_failure_tests libplCurve_cov.la
exit_failure_tests_LDADD = libplCurve_cov.la
exit_failure_tests_LDFLAGS = -static -fprofile-arcs
exit_failure_tests_CFLAGS = -g -Wall -fprofile-arcs -ftest-coverage

knot_diagram_SOURCE = knot_diagram.c libplCurve.la
knot_diagram_LDADD = libplCurve.la
