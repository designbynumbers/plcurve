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
	bad_vects/bad_14.vect \
	make_hopf

check_PROGRAMS = tsmcmc_random_test \
		 nn_test \
	         symmetry_test \
	         randompolygon_test \
	         eq_primitives_test \
	         loopclosure_test \
	         exit_failure_tests  \
	         test_struts \
                 test_minrad \
                 test_findpoca \
                 call_check \
                 hopf_test \
		 test_pdcode \
		 test_isomorphism \
		 test_pdcontainer \
	         test_multidx \
		 test_dihedral \
	         test_perm \
		 ccode_test

TESTS = tsmcmc_random_test \
	nn_test \
	symmetry_test \
	eq_primitives_test \
	loopclosure_test \
	test_findpoca \
        test_struts \
        call_check \
        hopf_test \
	test_pdcode \
	test_isomorphism \
	test_pdcontainer \
	test_multidx \
	test_dihedral \
	test_perm \
	ccode_test

ccode_test_SOURCES = ccode_test.c 
ccode_test_LDADD = -L$(top_builddir)/src -lplCurveStatic
ccode_test_LDFLAGS = -static
ccode_test_CFLAGS = -I$(top_srcdir)/src

tsmcmc_random_test_SOURCES = tsmcmc_random_test.c
tsmcmc_random_test_LDFLAGS = -static
tsmcmc_random_test_LDADD = -L$(top_builddir)/src -lplCurve
tsmcmc_random_test_CFLAGS = -I$(top_srcdir)/src 

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

CLEANFILES = runtest_wrapper.sh hopftest.vect treftest.vect

noinst_PROGRAMS = test_levels octrope_timing

test_levels_SOURCES = test_levels.c
test_levels_LDADD = -L$(top_builddir)/src -lplCurve
test_levels_CFLAGS = -I$(top_srcdir)/src

test_struts_SOURCES = test_struts.c
test_struts_LDADD = -L$(top_builddir)/src -lplCurve
test_struts_CFLAGS = -I$(top_srcdir)/src

test_minrad_SOURCES = test_minrad.c
test_minrad_LDADD = -L$(top_builddir)/src -lplCurve
test_minrad_CFLAGS = -I$(top_srcdir)/src

test_findpoca_SOURCES = test_findpoca.c
test_findpoca_LDADD = -L$(top_builddir)/src -lplCurve
test_findpoca_CFLAGS = -I$(top_srcdir)/src

call_check_SOURCES = call_check.c
call_check_LDADD = -L$(top_builddir)/src -lplCurve
call_check_CFLAGS = -I$(top_srcdir)/src

hopf_test_SOURCES = hopf_test.c
hopf_test_LDADD = -L$(top_builddir)/src -lplCurve
hopf_test_CFLAGS = -I$(top_srcdir)/src

octrope_timing_SOURCES = octrope_timing.c
octrope_timing_LDADD = -L$(top_builddir)/src -lplCurve
octrope_timing_CFLAGS = -I$(top_srcdir)/src

test_isomorphism_SOURCES = test_isomorphism.c
test_isomorphism_LDADD = -L$(top_builddir)/src -lplCurve
test_isomorphism_CFLAGS = -I$(top_srcdir)/src

test_pdcode_SOURCES = test_pdcode.c 
test_pdcode_LDADD = -L$(top_builddir)/src -lplCurve
test_pdcode_CFLAGS = -I$(top_srcdir)/src

test_pdcontainer_SOURCES = test_pdcontainer.c 
test_pdcontainer_LDADD = -L$(top_builddir)/src -lplCurve
test_pdcontainer_CFLAGS = -I$(top_srcdir)/src

test_multidx_SOURCES = test_multidx.c
test_multidx_LDADD = -L$(top_builddir)/src -lplCurve
test_multidx_CFLAGS = -I$(top_srcdir)/src

test_dihedral_SOURCES = test_dihedral.c
test_dihedral_LDADD = -L$(top_builddir)/src -lplCurve
test_dihedral_CFLAGS = -I$(top_srcdir)/src

test_perm_SOURCES = test_perm.c
test_perm_LDADD = -L$(top_builddir)/src -lplCurve
test_perm_CFLAGS = -I$(top_srcdir)/src

