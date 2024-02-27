AM_CFLAGS = -g -Wall -finline-functions -ftree-vectorize -Wall
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = processpdstor orientedshadows shadow2diagram knottype

processpdstor_SOURCES = processpdstor.c
processpdstor_LDADD = ${top_builddir}/src/libplCurve.la
processpdstor_LDFLAGS = -static
processpdstor_CFLAGS = -I${top_srcdir}/src

orientedshadows_SOURCES = orientedshadows.c
orientedshadows_LDADD = ${top_builddir}/src/libplCurve.la
orientedshadows_LDFLAGS = -static
orientedshadows_CFLAGS = -I${top_srcdir}/src

shadow2diagram_SOURCES = shadow2diagram.c
shadow2diagram_LDADD = ${top_builddir}/src/libplCurve.la
shadow2diagram_LDFLAGS = -static
shadow2diagram_CFLAGS = -I${top_srcdir}/src

knottype_SOURCES = knottype.c polynomials.c polynomials.h ordie.h ordie.c mangle.c mangle.h
knottype_LDADD = ${top_builddir}/src/libplCurve.la
knottype_LDFLAGS = -static
knottype_CFLAGS = -I${top_srcdir}/src
