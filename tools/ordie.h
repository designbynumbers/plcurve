/* 

  ordie.h : Versions of common library functions which terminate if they fail.

*/

#ifndef ORDIE_H__
#define ORDIE_H__

#include<config.h>

#ifdef HAVE_STDIO_H
#include<stdio.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#include<sys/types.h>
#endif

#define fopen_or_die(A,B) fopen_or_die_internal(A,B,__FILE__,__LINE__,stderr)

FILE *fopen_or_die_internal(const char *filename,const char *opentype,
			    const char *file,const int line,FILE *logstream);

#define mkdir_or_die(A,B) mkdir_or_die_internal(A,B,__FILE__,__LINE__,stderr)

void mkdir_or_die_internal(const char *filename,mode_t mode,
			   const char *file,const int line,FILE *logstream);

#define chmod_or_die(A,B) chmod_or_die_internal(A,B,__FILE__,__LINE__,stderr)

void chmod_or_die_internal(const char *filename,mode_t mode,
			   const char *file,const int line,FILE *logstream);


#define system_or_die(A) system_or_die_internal(A,__FILE__,__LINE__,stderr)

void system_or_die_internal(const char *cmdstring,
			    const char *file,const int line,FILE *logstream);

#define chdir_or_die(A) chdir_or_die_internal(A,__FILE__,__LINE__,stderr)

void chdir_or_die_internal(const char *dirname,
			   const char *file,const int line,FILE *logstream);

#define calloc_or_die(A,B) calloc_or_die_internal(A,B,__FILE__,__LINE__,stderr)

void *calloc_or_die_internal(size_t count, size_t eltsize,
			     const char *file,const int line,FILE *logstream);


#endif
