SUBDIRS = JudyCommon JudyL Judy1 JudySL JudyHS . # Reversed order of visiting current dir.

# Now generate the main library here -- JC.

noinst_LTLIBRARIES = libJudy.la
libJudy_la_SOURCES =                # We have to be explicit that there are no sources.
noinst_HEADERS = Judy.h
libJudy_la_LIBADD = ./JudyCommon/*.lo ./JudyL/*.lo ./Judy1/*.lo ./JudyHS/*.lo ./JudySL/*.lo 

DISTCLEANFILES = .deps Makefile
