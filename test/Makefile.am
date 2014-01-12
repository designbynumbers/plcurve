LDADD = libplCurve.la

check_PROGRAMS =  nn_test symmetry_test randompolygon_test eq_primitives_test loopclosure_test run_tests exit_failure_tests
TESTS = nn_test symmetry_test eq_primitives_test loopclosure_test run_tests

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

exit_failure_tests_SOURCES = exit_failure_tests.c
exit_failure_tests_LDADD = -L$(top_builddir)/src -lplCurve
exit_failure_tests_CFLAGS = -I$(top_builddir)/src

run_tests_SOURCES = run_tests.c
run_tests_LDADD = -L$(top_builddir)/src -lplCurve
run_tests_CFLAGS = -I$(top_builddir)/src
