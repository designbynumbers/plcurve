#LDADD = libplCurve.la

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
	rolfsentable.txt \
	thistlethwaitetable.txt
#	make_hopf

if HAVE_ARGTABLE2
check_PROGRAMS = csu_random_polygon_test \
		 test_hypersimplex_sampling \
		 test_k_ascent_permutations \
		 test_tangle_slide_ops \
		 test_splitdiagram \
		 test_readkt \
		 test_tangle_regenerate \
		 test_tangle_slide_input \
		 test_bad_homfly \
		 test_diagram_isotopic \
		 test_2crossing \
		 test_crossingops \
	         test_operations \
	         test_orientation \
		 test_homfly \
		 tsmcmc_random_test \
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

TESTS = csu_random_polygon_test \
	test_hypersimplex_sampling \
	test_k_ascent_permutations \
	test_tangle_slide_ops \
	test_splitdiagram \
	test_readkt \
	test_tangle_regenerate \
	test_tangle_slide_input \
	test_bad_homfly \
	test_diagram_isotopic \
	test_2crossing \
	test_crossingops \
	test_operations \
	test_orientation \
	test_homfly  \
	tsmcmc_random_test \
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
else
check_PROGRAMS = csu_random_polygon_test \
		 test_hypersimplex_sampling \
		 test_k_ascent_permutations \
		 test_tangle_slide_ops \
		 test_splitdiagram \
	         test_readkt \
	         test_tangle_regenerate \
		 test_tangle_slide_input \
		 test_bad_homfly \
		 test_diagram_isotopic \
		 test_2crossing \
		 test_crossingops \
	         test_operations \
		 test_homfly \
	         tsmcmc_random_test \
		 nn_test \
	         symmetry_test \
	         randompolygon_test \
	         eq_primitives_test \
	         loopclosure_test \
	         exit_failure_tests  \
	         test_struts \
                 test_minrad \
                 test_findpoca \
		 test_dihedral \
                 call_check \
                 hopf_test \
		 ccode_test

TESTS = csu_random_polygon_test \
	test_hypersimplex_sampling \
	test_k_ascent_permutations \
	test_tangle_slide_ops \
	test_splitdiagram \
	test_readkt \
	test_tangle_regenerate \
	test_tangle_slide_input \
	test_bad_homfly \
	test_diagram_isotopic \
	test_2crossing \
	test_crossingops \
	test_operations \
	test_homfly \
	tsmcmc_random_test \
	nn_test \
	symmetry_test \
	eq_primitives_test \
	loopclosure_test \
	test_findpoca \
        test_struts \
        call_check \
        hopf_test \
	test_dihedral \
	ccode_test
endif

test_homfly_SOURCES = test_homfly.c 
test_homfly_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_homfly_LDFLAGS = -static
test_homfly_CFLAGS = -I$(top_srcdir)/src

csu_random_polygon_test_SOURCES = csu_random_polygon_test.c 
csu_random_polygon_test_LDADD = -L$(top_builddir)/src -lplCurveStatic
csu_random_polygon_test_LDFLAGS = -static
csu_random_polygon_test_CFLAGS = -I$(top_srcdir)/src

test_k_ascent_permutations_SOURCES = test_k_ascent_permutations.c 
test_k_ascent_permutations_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_k_ascent_permutations_LDFLAGS = -static
test_k_ascent_permutations_CFLAGS = -I$(top_srcdir)/src

test_hypersimplex_sampling_SOURCES = test_hypersimplex_sampling.c 
test_hypersimplex_sampling_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_hypersimplex_sampling_LDFLAGS = -static
test_hypersimplex_sampling_CFLAGS = -I$(top_srcdir)/src

test_tangle_slide_ops_SOURCES = test_tangle_slide_ops.c 
test_tangle_slide_ops_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_tangle_slide_ops_LDFLAGS = -static
test_tangle_slide_ops_CFLAGS = -I$(top_srcdir)/src

test_splitdiagram_SOURCES = test_splitdiagram.c 
test_splitdiagram_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_splitdiagram_LDFLAGS = -static
test_splitdiagram_CFLAGS = -I$(top_srcdir)/src

test_readkt_SOURCES = test_readkt.c 
test_readkt_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_readkt_LDFLAGS = -static
test_readkt_CFLAGS = -I$(top_srcdir)/src

test_tangle_regenerate_SOURCES = test_tangle_regenerate.c 
test_tangle_regenerate_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_tangle_regenerate_LDFLAGS = -static
test_tangle_regenerate_CFLAGS = -I$(top_srcdir)/src

test_tangle_slide_input_SOURCES = test_tangle_slide_input.c 
test_tangle_slide_input_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_tangle_slide_input_LDFLAGS = -static
test_tangle_slide_input_CFLAGS = -I$(top_srcdir)/src

test_bad_homfly_SOURCES = test_bad_homfly.c 
test_bad_homfly_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_bad_homfly_LDFLAGS = -static
test_bad_homfly_CFLAGS = -I$(top_srcdir)/src

test_diagram_isotopic_SOURCES = test_diagram_isotopic.c 
test_diagram_isotopic_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_diagram_isotopic_LDFLAGS = -static
test_diagram_isotopic_CFLAGS = -I$(top_srcdir)/src

test_2crossing_SOURCES = test_2crossing.c 
test_2crossing_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_2crossing_LDFLAGS = -static
test_2crossing_CFLAGS = -I$(top_srcdir)/src

test_operations_SOURCES = test_operations.c 
test_operations_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_operations_LDFLAGS = -static
test_operations_CFLAGS = -I$(top_srcdir)/src

test_crossingops_SOURCES = test_crossingops.c 
test_crossingops_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_crossingops_LDFLAGS = -static
test_crossingops_CFLAGS = -I$(top_srcdir)/src

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
test_isomorphism_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_isomorphism_CFLAGS = -I$(top_srcdir)/src
test_isomorphism_LDFLAGS = -static

test_pdcode_SOURCES = test_pdcode.c 
test_pdcode_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_pdcode_CFLAGS = -I$(top_srcdir)/src
test_pdcode_LDFLAGS = -static

test_pdcontainer_SOURCES = test_pdcontainer.c 
test_pdcontainer_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_pdcontainer_CFLAGS = -I$(top_srcdir)/src
test_pdcontainer_LDFLAGS = -static

test_multidx_SOURCES = test_multidx.c
test_multidx_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_multidx_CFLAGS = -I$(top_srcdir)/src
test_multidx_LDFLAGS = -static

test_dihedral_SOURCES = test_dihedral.c
test_dihedral_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_dihedral_CFLAGS = -I$(top_srcdir)/src
test_dihedral_LDFLAGS = -static

test_perm_SOURCES = test_perm.c
test_perm_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_perm_CFLAGS = -I$(top_srcdir)/src
test_perm_LDFLAGS = -static

test_orientation_SOURCES = test_orientation.c
test_orientation_LDADD = -L$(top_builddir)/src -lplCurveStatic
test_orientation_CFLAGS = -I$(top_srcdir)/src
test_orientation_LDFLAGS = -static
