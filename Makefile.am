AM_CFLAGS = -g -Wall -fprofile-arcs -ftest-coverage

noinst_LIBRARIES = libplCurve.a
noinst_PROGRAMS = test_readwrite test_cst

libplCurve_a_SOURCES = plCurve.c plCurve.h vector.c \
                       spline_links.h spline_links.c 

test_readwrite_SOURCE = test_readwrite.c
test_readwrite_LDADD = libplCurve.a
test_cst_SOURCE = test_cst.c
test_cst_LDADD = libplCurve.a
