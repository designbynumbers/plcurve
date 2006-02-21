AM_CFLAGS = -g -Wall -fprofile-arcs -ftest-coverage

noinst_LIBRARIES = libplCurve.a libplCurve_safe.a
noinst_PROGRAMS = test_readwrite test_cst test_rand

libplCurve_a_SOURCES = plCurve.c plCurve.h vector.c splines.c 
libplCurve_a_CFLAGS = -DNDEBUG -O2 $(AM_CFLAGS)

libplCurve_safe_a_SOURCES = plCurve.c plCurve.h vector.c splines.c

test_readwrite_SOURCE = test_readwrite.c libplCurve.a
test_readwrite_LDADD = libplCurve_safe.a

test_cst_SOURCE = test_cst.c libplCurve.a
test_cst_LDADD = libplCurve_safe.a

test_rand_SOURCE = test_rand.c libplCurve.a
test_rand_LDADD = libplCurve_safe.a
test_rand_CFLAGS = -pg $(AM_CFLAGS)
test_rand_LDFLAGS = -pg $(AM_LDFLAGS)
