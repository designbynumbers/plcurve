INCLUDES =  -I$(srcdir)/. -I$(srcdir)/.. -I$(srcdir)/../JudyCommon/ 
AM_CFLAGS = -DJUDYL @WARN_CFLAGS@ 

include_HEADERS = JudyL.h
noinst_LTLIBRARIES = libJudyL.la libnext.la libprev.la libcount.la libinline.la

libJudyL_la_SOURCES = JudyLCascade.c JudyLTables.c JudyLCount.c JudyLCreateBranch.c JudyLDecascade.c JudyLDel.c JudyLFirst.c JudyLFreeArray.c JudyLGet.c JudyLInsArray.c JudyLIns.c JudyLInsertBranch.c JudyLMallocIF.c JudyLMemActive.c JudyLMemUsed.c 

libnext_la_SOURCES = JudyLNext.c JudyLNextEmpty.c
libnext_la_CFLAGS = $(AM_CFLAGS) -DJUDYNEXT

libprev_la_SOURCES = JudyLPrev.c JudyLPrevEmpty.c
libprev_la_CFLAGS = $(AM_CFLAGS) -DJUDYPREV

libcount_la_SOURCES = JudyLByCount.c
libcount_la_CFLAGS = $(AM_CFLAGS) -DNOSMARTJBB -DNOSMARTJBU -DNOSMARTJLB

libinline_la_SOURCES = j__udyLGet.c
libinline_la_CFLAGS = $(AM_CFLAGS) -DJUDYGETINLINE

JudyLTables.c: JudyLTablesGen.c
	$(CC) $(INCLUDES) $(AM_CFLAGS) @CFLAGS@ -o JudyLTablesGen JudyLTablesGen.c; ./JudyLTablesGen 

JudyLByCount.c:   $(srcdir)/../JudyCommon/JudyByCount.c
	cp -f $(srcdir)/../JudyCommon/JudyByCount.c      		JudyLByCount.c   

JudyLCascade.c:$(srcdir)/../JudyCommon/JudyCascade.c
	cp -f $(srcdir)/../JudyCommon/JudyCascade.c       	JudyLCascade.c

JudyLCount.c:$(srcdir)/../JudyCommon/JudyCount.c
	cp -f $(srcdir)/../JudyCommon/JudyCount.c        		JudyLCount.c

JudyLCreateBranch.c:$(srcdir)/../JudyCommon/JudyCreateBranch.c 
	cp -f $(srcdir)/../JudyCommon/JudyCreateBranch.c 		JudyLCreateBranch.c

JudyLDecascade.c:$(srcdir)/../JudyCommon/JudyDecascade.c    
	cp -f $(srcdir)/../JudyCommon/JudyDecascade.c    		JudyLDecascade.c

JudyLDel.c:$(srcdir)/../JudyCommon/JudyDel.c       
	cp -f $(srcdir)/../JudyCommon/JudyDel.c          		JudyLDel.c

JudyLFirst.c:$(srcdir)/../JudyCommon/JudyFirst.c
	cp -f $(srcdir)/../JudyCommon/JudyFirst.c        		JudyLFirst.c

JudyLFreeArray.c:$(srcdir)/../JudyCommon/JudyFreeArray.c    
	cp -f $(srcdir)/../JudyCommon/JudyFreeArray.c    		JudyLFreeArray.c

JudyLGet.c:$(srcdir)/../JudyCommon/JudyGet.c       
	cp -f $(srcdir)/../JudyCommon/JudyGet.c          		JudyLGet.c

j__udyLGet.c:$(srcdir)/../JudyCommon/JudyGet.c 
	cp -f $(srcdir)/../JudyCommon/JudyGet.c          		j__udyLGet.c

JudyLInsArray.c:$(srcdir)/../JudyCommon/JudyInsArray.c     
	cp -f $(srcdir)/../JudyCommon/JudyInsArray.c     		JudyLInsArray.c

JudyLIns.c:$(srcdir)/../JudyCommon/JudyIns.c      
	cp -f $(srcdir)/../JudyCommon/JudyIns.c          		JudyLIns.c

JudyLInsertBranch.c:$(srcdir)/../JudyCommon/JudyInsertBranch.c
	cp -f $(srcdir)/../JudyCommon/JudyInsertBranch.c 		JudyLInsertBranch.c

JudyLMallocIF.c:$(srcdir)/../JudyCommon/JudyMallocIF.c     
	cp -f $(srcdir)/../JudyCommon/JudyMallocIF.c     		JudyLMallocIF.c

JudyLMemActive.c:$(srcdir)/../JudyCommon/JudyMemActive.c
	cp -f $(srcdir)/../JudyCommon/JudyMemActive.c    		JudyLMemActive.c

JudyLMemUsed.c:$(srcdir)/../JudyCommon/JudyMemUsed.c   
	cp -f $(srcdir)/../JudyCommon/JudyMemUsed.c      		JudyLMemUsed.c

JudyLNext.c:$(srcdir)/../JudyCommon/JudyPrevNext.c
	cp -f $(srcdir)/../JudyCommon/JudyPrevNext.c     		JudyLNext.c

JudyLPrev.c:$(srcdir)/../JudyCommon/JudyPrevNext.c
	cp -f $(srcdir)/../JudyCommon/JudyPrevNext.c     		JudyLPrev.c

JudyLNextEmpty.c:$(srcdir)/../JudyCommon/JudyPrevNextEmpty.c
	cp -f $(srcdir)/../JudyCommon/JudyPrevNextEmpty.c		JudyLNextEmpty.c

JudyLPrevEmpty.c:$(srcdir)/../JudyCommon/JudyPrevNextEmpty.c
	cp -f $(srcdir)/../JudyCommon/JudyPrevNextEmpty.c		JudyLPrevEmpty.c

JudyLTablesGen.c:$(srcdir)/../JudyCommon/JudyTables.c
	cp -f $(srcdir)/../JudyCommon/JudyTables.c       	 	JudyLTablesGen.c

DISTCLEANFILES = .deps Makefile 

CLEANFILES = JudyLByCount.c  \
	     JudyLCascade.c \
	     JudyLCount.c \
	     JudyLCreateBranch.c \
	     JudyLDecascade.c \
	     JudyLDel.c \
	     JudyLFirst.c \
	     JudyLFreeArray.c \
	     JudyLGet.c \
	     j__udyLGet.c \
	     JudyLInsArray.c \
	     JudyLIns.c \
	     JudyLInsertBranch.c \
	     JudyLMallocIF.c \
	     JudyLMemActive.c \
	     JudyLMemUsed.c \
	     JudyLNext.c \
	     JudyLPrev.c \
	     JudyLNextEmpty.c \
	     JudyLPrevEmpty.c \
	     JudyLTablesGen.c \
	     JudyLTables.c \
	     JudyLTablesGen \
	     .libs \
	     *.o \
	     *.lo \
	     *.la 
