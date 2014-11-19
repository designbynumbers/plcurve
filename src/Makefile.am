AM_CFLAGS = -g -Wall #-finline-functions -ftree-vectorize -Wall
ACLOCAL_AMFLAGS = -I m4

include_HEADERS = plCurve.h nplCurve.h octrope.h plcTopology.h plcRandomPolygon.h
lib_LTLIBRARIES = libplCurve.la libnplCurve.la libplCurveStatic.la

if HAVE_ARGTABLE2
bin_PROGRAMS = knot_diagram resample plcurvature ropelength struts roundout_rl mrlocs plcrwalk momentpolytope carbonizepd
else
bin_PROGRAMS = knot_diagram
endif

noinst_PROGRAMS = spline_view # gcov_tests
noinst_HEADERS = homfly.h

libplCurveStatic_la_SOURCES = plCurve.c \
	plCurve.h \
	vector.c \
	splines.c \
	pllmpoly02.c \
	plcTopology.c \
	plwhitten.c \
	nearest_neighbor.c \
	symmetries.c \
	matrix.c \
	matrix.h \
	plcRandomPolygon.c \
	plc_dihedrals.c \
	octrope.h \
	octrope_internal.h \
        octrope.c \
        octmem.h \
        octmem.c \
        minrad.c \
        octcnv2.h \
        octcnv4.h \
	ordie.h \
	ordie.c \
	pdcode.c \
	pd_isomorphisms.h \
	pd_isomorphisms.c \
	pd_container.h \
	pd_container.c \
	pd_multidx.h \
	pd_multidx.c \
	pdcode_testcodes.c \
	pd_dihedral.h \
	pd_dihedral.c \
	pd_cyclic.h \
	pd_cyclic.c \
	pd_orientation.h \
	pd_orientation.c \
	pd_permdata.h \
	pd_perm.h \
	pd_perm.c \
	plcToricSymplectic.c \
	plcToricSymplecticTriangulation.c \
	tsmcmc.h \
	geyer_estimators.h \
	geyer_estimators.c \
	plcRandomPolygon.h \
	polynomials.h \
	polynomials.c \
	pd_operations.c \
	pd_crossingmoves.c \
	pd_sortedbuf.h \
	pd_sortedbuf.c \
	pd_tangles.c \
	pd_splitdiagram.h \
	pd_splitdiagram.c \
	pd_deletions.h \
	pd_deletions.c

libplCurveStatic_la_CFLAGS = $(AM_CFLAGS)
libplCurveStatic_la_LDFLAGS = -version-info $(LIB_VERSION) -static

libplCurve_la_SOURCES = plCurve.c \
	plCurve.h \
	vector.c \
	splines.c \
	pllmpoly02.c \
	plcTopology.c \
	plwhitten.c \
	nearest_neighbor.c \
	symmetries.c \
	matrix.c \
	matrix.h \
	plcRandomPolygon.c \
	plc_dihedrals.c \
	octrope.h \
	octrope_internal.h \
        octrope.c \
        octmem.h \
        octmem.c \
        minrad.c \
        octcnv2.h \
        octcnv4.h \
	ordie.h \
	ordie.c \
	pdcode.c \
	pd_isomorphisms.h \
	pd_isomorphisms.c \
	pd_container.h \
	pd_container.c \
	pd_multidx.h \
	pd_multidx.c \
	pdcode_testcodes.c \
	pd_dihedral.h \
	pd_dihedral.c \
	pd_cyclic.h \
	pd_cyclic.c \
	pd_orientation.h \
	pd_orientation.c \
	pd_permdata.h \
	pd_perm.h \
	pd_perm.c \
	plcToricSymplectic.c \
	plcToricSymplecticTriangulation.c \
	tsmcmc.h \
	geyer_estimators.h \
	geyer_estimators.c \
	plcRandomPolygon.h \
	polynomials.h \
	polynomials.c \
	pd_operations.c \
	pd_crossingmoves.c \
	pd_sortedbuf.h \
	pd_sortedbuf.c \
	pd_tangles.c \
	pd_splitdiagram.h \
	pd_splitdiagram.c \
	pd_deletions.h \
	pd_deletions.c

libplCurve_la_CFLAGS = $(AM_CFLAGS)
libplCurve_la_LDFLAGS = -version-info $(LIB_VERSION)

libnplCurve_la_SOURCES = nplCurve.c nplCurve.h nvector.c
libnplCurve_la_CFLAGS = $(AM_CFLAGS)
libnplCurve_la_LDFLAGS = -version-info $(LIB_VERSION)

knot_diagram_SOURCES = knot_diagram.c
knot_diagram_LDADD = -lplCurve
if HAVE_ARGTABLE2
knot_diagram_LDFLAGS = -static
else
knot_diagram_LDFLAGS = -static
endif

#ccode_test_SOURCES = ccode_test.c
#ccode_test_LDADD = -lplCurve
#ccode_test_LDFLAGS = -static

momentpolytope_SOURCES = momentpolytope.c
momentpolytope_LDADD = -lplCurve

carbonizepd_SOURCES = carbonizepd.c
carbonizepd_LDADD = -lplCurve 
carbonizepd_LDFLAGS = -static
carbonizepd_CFLAGS = -I${top_builddir}/src

spline_view_SOURCES = spline_view.c
spline_view_LDADD = -lplCurve
spline_view_LDFLAGS = -static

resample_SOURCES = resample.c
resample_LDADD = -lplCurve
resample_LDFLAGS = -static

plcurvature_SOURCES = plcurvature.c
plcurvature_LDADD = -lplCurve
plcurvature_LDFLAGS = -static

plcrwalk_SOURCES = plcrwalk.c
plcrwalk_LDADD = -lplCurve

ropelength_SOURCES = ropelength.c
ropelength_LDADD = -lplCurve

struts_SOURCES = struts.c
struts_LDADD = -lplCurve

roundout_rl_SOURCES = roundout_rl.c
roundout_rl_LDADD = -lplCurve

mrlocs_SOURCES = mrlocs.c
mrlocs_LDADD = -lplCurve

if HAVE_SWIG
BUILT_SOURCES = plCurve_wrap.c \
		plcRandomPolygon_wrap.c
SWIG_SOURCES = plCurve.i plcRandomPolygon.i
PYTHON_GENERATED = plcurve.py tsmcmc.py

pkgpython_PYTHON = $(PYTHON_GENERATED) __init__.py
pkgpyexec_LTLIBRARIES = _plcurve.la _tsmcmc.la
_plcurve_la_SOURCES = plCurve_wrap.c plCurve.i
_plcurve_la_CFLAGS = $(AX_SWIG_PYTHON_CPPFLAGS) -I"$(top_builddir)" -I$(NUMPY_INCDIR)
_plcurve_la_LDFLAGS = -module
_plcurve_la_LIBADD = libplCurve.la
_tsmcmc_la_SOURCES = plcRandomPolygon_wrap.c plcRandomPolygon.i
_tsmcmc_la_CFLAGS = $(AX_SWIG_PYTHON_CPPFLAGS) -I"$(top_builddir)" -I$(NUMPY_INCDIR)
_tsmcmc_la_LDFLAGS = -module
_tsmcmc_la_LIBADD = libplCurve.la

plCurve_wrap.c : plCurve.i
	$(SWIG) $(AX_SWIG_PYTHON_OPT) -builtin -I"$(top_builddir)" -o $@ $<
plcRandomPolygon_wrap.c: plcRandomPolygon.i
	$(SWIG) $(AX_SWIG_PYTHON_OPT) -builtin -I"$(top_builddir)" -o $@ $<

clean-local:
	-rm $(BUILT_SOURCES) $(PYTHON_GENERATED)
endif
