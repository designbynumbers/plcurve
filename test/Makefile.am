LDADD = libplCurve.la

check_PROGRAMS =  nn_test symmetry_test randompolygon_test eq_primitives_test loopclosure_test
TESTS = nn_test symmetry_test randompolygon_test eq_primitives_test loopclosure_test

nn_test_SOURCES = nn_test.c
nn_test_LDADD = -L$(top_builddir)/src -lplCurve
nn_test_CFLAGS = -I$(top_builddir)/src

symmetry_test_SOURCES = symmetry_test.c
symmetry_test_LDADD = -L$(top_builddir)/src -lplCurve
symmetry_test_CFLAGS = -I$(top_builddir)/src

eq_primitives_test_SOURCES = eq_primitives_test.c 
eq_primitives_test_LDADD = -L$(top_builddir)/src -lplCurve
eq_primitives_test_CFLAGS = -I$(top_builddir)/src

loopclosure_test_SOURCES = loopclosure_test.c
loopclosure_test_LDADD = -L$(top_builddir)/src -lplCurve
loopclosure_test_CFLAGS = -I$(top_builddir)/src

randompolygon_test_SOURCES = randompolygon_test.c 
randompolygon_test_LDADD = -L$(top_builddir)/src -lplCurve
randompolygon_test_CFLAGS = -I$(top_builddir)/src
