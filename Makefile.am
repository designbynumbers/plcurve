AM_CFLAGS = -g -O2 -finline-functions -Wall

lib_LTLIBRARIES = libplCurve.la libplCurve_safe.la
noinst_LTLIBRARIES = libplCurve_cov.la
noinst_PROGRAMS = test_readwrite test_cst test_rand run_tests

libplCurve_la_SOURCES = plCurve.c plCurve.h vector.c splines.c 
libplCurve_la_CFLAGS = -DNDEBUG $(AM_CFLAGS)
libplCurve_la_LDFLAGS = -version-info $(LIB_VERSION)

libplCurve_safe_la_SOURCES = plCurve.c plCurve.h vector.c splines.c
libplCurve_safe_la_LDFLAGS = -version-info $(LIB_VERSION)

libplCurve_cov_la_SOURCES = plCurve.c plCurve.h vector.c splines.c
libplCurve_cov_la_CFLAGS = -g -Wall -fprofile-arcs -ftest-coverage
libplCurve_cov_la_LDFLAGS = -version-info $(LIB_VERSION)

test_readwrite_SOURCE = test_readwrite.c libplCurve_safe.la
test_readwrite_LDADD = libplCurve_safe.la

test_cst_SOURCE = test_cst.c libplCurve_safe.la
test_cst_LDADD = libplCurve_safe.la

test_rand_SOURCE = test_rand.c libplCurve_safe.la
test_rand_LDADD = libplCurve_safe.la
test_rand_CFLAGS = -pg $(AM_CFLAGS)
test_rand_LDFLAGS = -pg $(AM_LDFLAGS)

run_tests_SOURCE = run_tests.c libplCurve_cov.la
run_tests_LDADD = libplCurve_cov.la
run_tests_CFLAGS = -g -Wall -fprofile-arcs -ftest-coverage
