LDADD = libplCurve.la

check_PROGRAMS =  nn_test symmetry_test randompolygon_test eq_primitives_test loopclosure_test
TESTS = nn_test symmetry_test randompolygon_test eq_primitives_test loopclosure_test

nn_test_SOURCES = nn_test.c libplCurve.a
nn_test_LDADD = -lplCurve

symmetry_test_SOURCES = symmetry_test.c libplCurve.a
symmetry_test_LDADD = -lplCurve

eq_primitives_test_SOURCES = eq_primitives_test.c libplCurve.a
eq_primitives_test_LDADD = -lplCurve

loopclosure_test_SOURCES = loopclosure_test.c libplCurve.a
loopclosure_test_LDADD = -lplCurve
