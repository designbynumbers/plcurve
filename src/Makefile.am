AM_CFLAGS = -g -finline-functions -ftree-vectorize -Wall
ACLOCAL_AMFLAGS = -I m4

include_HEADERS = plCurve.h nplCurve.h
lib_LTLIBRARIES = libplCurve.la libnplCurve.la libplCurveStatic.la

if HAVE_ARGTABLE2
bin_PROGRAMS = knot_diagram resample plcurvature ropelength struts roundout_rl mrlocs 
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
        octcnv4.h 

libplCurveStatic_la_CFLAGS = -DNDEBUG  $(AM_CFLAGS)
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
        octcnv4.h 

libplCurve_la_CFLAGS = -DNDEBUG  $(AM_CFLAGS)
libplCurve_la_LDFLAGS = -version-info $(LIB_VERSION)

libnplCurve_la_SOURCES = nplCurve.c nplCurve.h nvector.c
libnplCurve_la_CFLAGS = -DNDEBUG  $(AM_CFLAGS)
libnplCurve_la_LDFLAGS = -version-info $(LIB_VERSION)

knot_diagram_SOURCES = knot_diagram.c 
knot_diagram_LDADD = -lplCurve
if HAVE_ARGTABLE2
knot_diagram_LDFLAGS = -static
else
knot_diagram_LDFLAGS = -static
endif

spline_view_SOURCES = spline_view.c 
spline_view_LDADD = -lplCurve
spline_view_LDFLAGS = -static

resample_SOURCES = resample.c 
resample_LDADD = -lplCurve
resample_LDFLAGS = -static

plcurvature_SOURCES = plcurvature.c 
plcurvature_LDADD = -lplCurve
plcurvature_LDFLAGS = -static

ropelength_SOURCES = ropelength.c
ropelength_LDADD = -lplCurve

struts_SOURCES = struts.c 
struts_LDADD = -lplCurve

roundout_rl_SOURCES = roundout_rl.c 
roundout_rl_LDADD = -lplCurve

mrlocs_SOURCES = mrlocs.c
mrlocs_LDADD = -lplCurve

# There's a sort of aborted project to provide real Python wrappers for this
# entire library using SWIG (eventually, future extensions to the library 
# will be build natively in Python, but we're not there yet).
# At the moment, all this is commented out.

#if HAVE_SWIG
#BUILT_SOURCES = plCurve_wrap.c
#SWIG_SOURCES = plCurve.i

#pkgpython_PYTHON = plcurve.py __init__.py
#pkgpyexec_LTLIBRARIES = _plcurve.la
#_plcurve_la_SOURCES = plCurve_wrap.c $(SWIG_SOURCES)
#_plcurve_la_CFLAGS = $(AX_SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)
#_plcurve_la_LDFLAGS = -module
#_plcurve_la_LIBADD = libplCurve.la

#plCurve_wrap.c : $(SWIG_SOURCES)
#	$(SWIG) $(AX_SWIG_PYTHON_OPT) -builtin -I$(top_srcdir) -o $@ $<
#endif
