AM_CFLAGS = -Wall

noinst_LIBRARIES = libplCurve.a
noinst_PROGRAMS = test_readwrite

libplCurve_a_SOURCES = plCurve.c plCurve.h vector.c \
                       spline_links.h spline_links.c 

test_readwrite_SOURCE = test_readwrite.c
test_readwrite_LDADD = libplCurve.a
