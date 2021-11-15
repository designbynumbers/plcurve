# was INCLUDES --JC
AM_CPPFLAGS =  -I$(srcdir)/. -I$(srcdir)/.. -I$(srcdir)/../JudyCommon/ 
AM_CFLAGS = -DJUDY1 @WARN_CFLAGS@ 

include_HEADERS = Judy1.h
noinst_LTLIBRARIES = libJudy1.la libnext.la libprev.la libcount.la libinline.la

libJudy1_la_SOURCES = Judy1Test.c Judy1Tables.c Judy1Set.c Judy1SetArray.c Judy1Unset.c Judy1Cascade.c Judy1Count.c Judy1CreateBranch.c Judy1Decascade.c Judy1First.c Judy1FreeArray.c Judy1InsertBranch.c Judy1MallocIF.c Judy1MemActive.c Judy1MemUsed.c 

libnext_la_SOURCES = Judy1Next.c Judy1NextEmpty.c 
libnext_la_CFLAGS = $(AM_CFLAGS) -DJUDYNEXT

libprev_la_SOURCES = Judy1Prev.c Judy1PrevEmpty.c 
libprev_la_CFLAGS = $(AM_CFLAGS) -DJUDYPREV

libcount_la_SOURCES = Judy1ByCount.c
libcount_la_CFLAGS = $(AM_CFLAGS) -DNOSMARTJBB -DNOSMARTJBU -DNOSMARTJLB

libinline_la_SOURCES = j__udy1Test.c
libinline_la_CFLAGS = $(AM_CFLAGS) -DJUDYGETINLINE

Judy1Tables.c: Judy1TablesGen.c
	$(CC) $(INCLUDES) $(AM_CFLAGS) @CFLAGS@ -o Judy1TablesGen Judy1TablesGen.c; ./Judy1TablesGen

Judy1ByCount.c:$(srcdir)/../JudyCommon/JudyByCount.c
	cp -f $(srcdir)/../JudyCommon/JudyByCount.c      		Judy1ByCount.c   

Judy1Cascade.c:$(srcdir)/../JudyCommon/JudyCascade.c
	cp -f $(srcdir)/../JudyCommon/JudyCascade.c      		Judy1Cascade.c

Judy1Count.c:$(srcdir)/../JudyCommon/JudyCount.c
	cp -f $(srcdir)/../JudyCommon/JudyCount.c        		Judy1Count.c

Judy1CreateBranch.c:$(srcdir)/../JudyCommon/JudyCreateBranch.c
	cp -f $(srcdir)/../JudyCommon/JudyCreateBranch.c 		Judy1CreateBranch.c

Judy1Decascade.c:$(srcdir)/../JudyCommon/JudyDecascade.c
	cp -f $(srcdir)/../JudyCommon/JudyDecascade.c    		Judy1Decascade.c

Judy1Unset.c:$(srcdir)/../JudyCommon/JudyDel.c
	cp -f $(srcdir)/../JudyCommon/JudyDel.c          		Judy1Unset.c

Judy1First.c:$(srcdir)/../JudyCommon/JudyFirst.c
	cp -f $(srcdir)/../JudyCommon/JudyFirst.c        		Judy1First.c

Judy1FreeArray.c:$(srcdir)/../JudyCommon/JudyFreeArray.c
	cp -f $(srcdir)/../JudyCommon/JudyFreeArray.c    		Judy1FreeArray.c

Judy1Test.c:$(srcdir)/../JudyCommon/JudyGet.c
	cp -f $(srcdir)/../JudyCommon/JudyGet.c          		Judy1Test.c

j__udy1Test.c:$(srcdir)/../JudyCommon/JudyGet.c
	cp -f $(srcdir)/../JudyCommon/JudyGet.c          		j__udy1Test.c

Judy1SetArray.c:$(srcdir)/../JudyCommon/JudyInsArray.c
	cp -f $(srcdir)/../JudyCommon/JudyInsArray.c     		Judy1SetArray.c

Judy1Set.c:$(srcdir)/../JudyCommon/JudyIns.c
	cp -f $(srcdir)/../JudyCommon/JudyIns.c          		Judy1Set.c

Judy1InsertBranch.c:$(srcdir)/../JudyCommon/JudyInsertBranch.c
	cp -f $(srcdir)/../JudyCommon/JudyInsertBranch.c 		Judy1InsertBranch.c

Judy1MallocIF.c:$(srcdir)/../JudyCommon/JudyMallocIF.c
	cp -f $(srcdir)/../JudyCommon/JudyMallocIF.c     		Judy1MallocIF.c

Judy1MemActive.c:$(srcdir)/../JudyCommon/JudyMemActive.c
	cp -f $(srcdir)/../JudyCommon/JudyMemActive.c    		Judy1MemActive.c

Judy1MemUsed.c:$(srcdir)/../JudyCommon/JudyMemUsed.c
	cp -f $(srcdir)/../JudyCommon/JudyMemUsed.c      		Judy1MemUsed.c

Judy1Next.c:$(srcdir)/../JudyCommon/JudyPrevNext.c
	cp -f $(srcdir)/../JudyCommon/JudyPrevNext.c     		Judy1Next.c

Judy1Prev.c:$(srcdir)/../JudyCommon/JudyPrevNext.c
	cp -f $(srcdir)/../JudyCommon/JudyPrevNext.c     		Judy1Prev.c

Judy1NextEmpty.c:$(srcdir)/../JudyCommon/JudyPrevNextEmpty.c
	cp -f $(srcdir)/../JudyCommon/JudyPrevNextEmpty.c		Judy1NextEmpty.c

Judy1PrevEmpty.c:$(srcdir)/../JudyCommon/JudyPrevNextEmpty.c
	cp -f $(srcdir)/../JudyCommon/JudyPrevNextEmpty.c		Judy1PrevEmpty.c

Judy1TablesGen.c:$(srcdir)/../JudyCommon/JudyTables.c
	cp -f $(srcdir)/../JudyCommon/JudyTables.c	        Judy1TablesGen.c

DISTCLEANFILES = .deps Makefile 

CLEANFILES = Judy1ByCount.c \
	     Judy1Cascade.c \
	     Judy1Count.c \
	     Judy1CreateBranch.c \
	     Judy1Decascade.c \
	     Judy1Unset.c \
	     Judy1First.c \
	     Judy1FreeArray.c \
	     Judy1Test.c \
	     j__udy1Test.c \
	     Judy1SetArray.c \
	     Judy1Set.c \
	     Judy1InsertBranch.c \
	     Judy1MallocIF.c \
	     Judy1MemActive.c \
	     Judy1MemUsed.c \
	     Judy1Next.c \
	     Judy1Prev.c \
	     Judy1NextEmpty.c \
	     Judy1PrevEmpty.c \
	     Judy1TablesGen.c \
	     Judy1Tables.c \
	     .libs \
	     Judy1TablesGen \
	     *.o \
	     *.lo \
	     *.la
