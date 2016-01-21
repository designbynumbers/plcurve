INCLUDES =  -I$(srcdir)/. -I$(srcdir)/..
AM_CFLAGS = @CFLAGS@ @WARN_CFLAGS@ 

include_HEADERS = JudyPrivate.h \
	JudyPrivate1L.h \
	JudyPrivateBranch.h \
	JudyByCount.c \
	JudyDecascade.c \
	JudyGet.c \
	JudyPrevNext.c \
	JudyCascade.c \
	JudyDel.c \
	JudyIns.c \
	JudyMallocIF.c \
	JudyPrevNextEmpty.c \
	JudyCount.c \
	JudyFirst.c \
	JudyInsArray.c \
	JudyMemActive.c \
	JudyPrintJP.c \
	JudyTables.c \
	JudyCreateBranch.c \
	JudyFreeArray.c \
	JudyInsertBranch.c \
	JudyMemUsed.c

noinst_LTLIBRARIES = libJudyMalloc.la

libJudyMalloc_la_SOURCES = JudyMalloc.c

DISTCLEANFILES = .deps Makefile
