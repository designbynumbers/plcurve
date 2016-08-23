lib_LTLIBRARIES = libplanarmap.la

bin_PROGRAMS = planarmap

include_HEADERS = PMdef.h PMconjugation.h PMenlight.h PMextract.h PMplanmap.h PMdisplay.h stats.h
libplanarmap_la_SOURCES = PMdef.c \
	PMconjugation.c \
	PMenlight.c \
	PMextract.c \
	PMdisplay.c \
	PMplanmap.c \
	interface.c \
	stats.c \
	PMdef.h \
	PMconjugation.h \
	PMenlight.h \
	PMextract.h \
	PMdisplay.h \
	PMplanmap.h \
	interface.h \
	stats.h

libplanarmap_la_CFLAGS = $(AM_CFLAGS)
libplanarmap_la_LDFLAGS = -version-info 1:2
libplanarmap_la_LIBADD = -lm

planarmap_SOURCES = PMdef.c \
	PMconjugation.c \
	PMenlight.c \
	PMextract.c \
	PMdisplay.c \
	PMplanmap.c \
	interface.c \
	stats.c \
	PMdef.h \
	PMconjugation.h \
	PMenlight.h \
	PMextract.h \
	PMdisplay.h \
	PMplanmap.h \
	interface.h \
	stats.h \
	planarmap.c
planarmap_CFLAGS = $(AM_CFLAGS) -W -Wall
planarmap_LDADD = -lm
