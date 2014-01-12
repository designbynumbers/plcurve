LDADD = libplCurve.la

EXTRA_DIST = bad_vects/bad_1.vect \
	bad_vects/bad_2.vect \
	bad_vects/bad_3.vect \
	bad_vects/bad_4.vect \
	bad_vects/bad_5.vect \
	bad_vects/bad_6.vect \
	bad_vects/bad_7.vect \
	bad_vects/bad_8.vect \
	bad_vects/bad_9.vect \
	bad_vects/bad_10.vect \
	bad_vects/bad_11.vect \
	bad_vects/bad_12.vect \
	bad_vects/bad_13.vect \
	bad_vects/bad_14.vect 

check_PROGRAMS =  nn_test symmetry_test randompolygon_test eq_primitives_test loopclosure_test run_tests exit_failure_tests
TESTS = nn_test symmetry_test eq_primitives_test loopclosure_test runtest_wrapper.sh

nn_test_SOURCES = nn_test.c
nn_test_LDADD = -L$(top_builddir)/src -lplCurve
nn_test_CFLAGS = -I$(top_srcdir)/src

symmetry_test_SOURCES = symmetry_test.c
symmetry_test_LDADD = -L$(top_builddir)/src -lplCurve
symmetry_test_CFLAGS = -I$(top_srcdir)/src

eq_primitives_test_SOURCES = eq_primitives_test.c 
eq_primitives_test_LDADD = -L$(top_builddir)/src -lplCurve
eq_primitives_test_CFLAGS = -I$(top_srcdir)/src

loopclosure_test_SOURCES = loopclosure_test.c
loopclosure_test_LDADD = -L$(top_builddir)/src -lplCurve
loopclosure_test_CFLAGS = -I$(top_srcdir)/src

randompolygon_test_SOURCES = randompolygon_test.c 
randompolygon_test_LDADD = -L$(top_builddir)/src -lplCurve
randompolygon_test_CFLAGS = -I$(top_srcdir)/src

exit_failure_tests_SOURCES = exit_failure_tests.c
exit_failure_tests_LDADD = -L$(top_builddir)/src -lplCurve
exit_failure_tests_CFLAGS = -I$(top_srcdir)/src

run_tests_SOURCES = run_tests.c
run_tests_LDADD = -L$(top_builddir)/src -lplCurve
run_tests_CFLAGS = -I$(top_srcdir)/src

runtest_wrapper.sh:
	echo './run_tests $(top_srcdir)/test/' > runtest_wrapper.sh
	chmod +x runtest_wrapper.sh

CLEANFILES = runtest_wrapper.sh
