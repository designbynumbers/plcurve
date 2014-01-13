/*

   ordie.c : Versions of common library functions which terminate if they fail.

*/

#include"config.h"

#ifdef HAVE_STDIO_H
  #include<stdio.h>
#endif

#ifdef HAVE_STDLIB_H
  #include<stdlib.h>
#endif

#ifdef HAVE_SYS_STAT_H
  #include<sys/stat.h>
#endif

#ifdef HAVE_ERRNO_H
  #include<errno.h>
#endif

#ifdef HAVE_STRING_H
  #include<string.h>
#endif

#ifdef HAVE_UNISTD_H
  #include<unistd.h>
#endif

#include"ordie.h"

FILE *fopen_or_die_internal(const char *filename,const char *opentype,
			    const char *file,const int line,FILE *logstream)

{

  FILE *retfile;

  retfile = fopen(filename,opentype);

  if (retfile == NULL) {

    fprintf(stderr,"%s (%d): Failed to open file %s.\n",file,line,filename);
    exit(1);

  }

  return retfile;


}

void mkdir_or_die_internal(const char *filename,mode_t mode,
			   const char *file,const int line,FILE *logstream)

{
  int errnum;
  
  errnum = mkdir(filename,mode);

  if (errnum != 0) {

    fprintf(logstream,
	    "%s (%d): Failed to make directory %s.\n"
	    "         %s.\n",file,line,filename,strerror(errnum));

    exit(1);
  
  }

}

void chmod_or_die_internal(const char *filename,mode_t mode,
			   const char *file,const int line,FILE *logstream)

{
  int errnum;
  
  errnum = chmod(filename,mode);

  if (errnum != 0) {

    fprintf(logstream,
	    "%s (%d): Failed to change permissions of %s to %d.\n"
	    "         %s.\n",file,line,filename,mode,strerror(errnum));

    exit(1);
  
  }

}

void system_or_die_internal(const char *cmdstring,
			    const char *file,const int line,FILE *logstream)
{
  int errnum;
  
  errnum = system(cmdstring);

  if (errnum != 0) {

    fprintf(logstream,
	    "%s (%d): Failed to execute system call %s.\n"
	    ,file,line,cmdstring);

    exit(1);
  
  }
  
}

void chdir_or_die_internal(const char *dirname,
			   const char *file,const int line,FILE *logstream)

{
  int errnum;
  
  errnum = chdir(dirname);

  if (errnum != 0) {

    fprintf(logstream,
	    "%s (%d): Failed to chdir to directory %s.\n"
	    "         %s.\n",file,line,dirname,strerror(errnum));

    exit(1);
  
  }

}

void *malloc_or_die_internal(size_t size,const char *file,const int line,FILE *logstream)
{
  void *ptr;

  ptr = malloc(size);

  if (ptr == NULL) {

    fprintf(logstream,
	    "%s (%d): Failed to malloc a buffer of size %d.\n"
	    ,file,line,(int)(size));

    exit(1);
  
  }

  return ptr;

}

void *calloc_or_die_internal(size_t count, size_t eltsize,
			    const char *file,
			    const int line,FILE *logstream)

{
  void *ptr;
  
  ptr= calloc(count,eltsize);

  if (ptr == NULL) {

    fprintf(logstream,
	    "%s (%d): Failed to calloc %d elements of size %d.\n"
	    ,file,line,(int)(count),(int)(eltsize));

    exit(1);
  
  }

  return ptr;

}

void *realloc_or_die_internal(void *ptr, size_t newsize,
			     const char *file,const int line,FILE *logstream)

{
  void *newptr;
  
  newptr= realloc(ptr,newsize);

  if (newptr == NULL) {

    fprintf(logstream,
	    "%s (%d): Failed to realloc buffer to size %d.\n"
	    ,file,line,(int)(newsize));

    exit(1);
  
  }

  return ptr;

}
