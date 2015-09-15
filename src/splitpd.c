/* 

   splitpd.c : reads a pdstor and divides into sequentially numbered pieces.

*/

#ifdef HAVE_CONFIG_H
#include"config.h"
#endif

#ifdef HAVE_STDIO_H
#include<stdio.h>
#endif

#ifdef HAVE_STDBOOL_H
#include<stdbool.h>
#endif

#ifdef HAVE_STRING_H
#include<string.h>
#endif

#ifdef HAVE_STDINT_H
#include<stdint.h>
#endif

#ifdef HAVE_STDLIB_H
#include<stdlib.h>
#endif

#ifdef HAVE_MATH_H
#include<math.h>
#endif

#ifdef HAVE_TIME_H
#include<time.h>
#endif

#ifdef HAVE_ASSERT_H
#include<assert.h>
#endif

#ifdef HAVE_LIMITS_H
#include<limits.h>
#endif

#ifdef HAVE_INTTYPES_H
#include<inttypes.h>
#endif

#ifdef HAVE_TIME_H
#include<time.h>
#endif


#include"plcTopology.h"
#include"pd_multidx.h"
#include"pd_perm.h"
#include"pd_isomorphisms.h"
#include"pd_storage.h"
#include"pd_orientation.h"

#include"ordie.h"
#include"argtable2.h"
/* We use a local, compiled in copy of argtable to avoid having to
   depend on system utilities. */

struct arg_lit  *verbose;
struct arg_file *file;
struct arg_int  *pieces;
struct arg_int  *size;
struct arg_lit  *help;
struct arg_end  *end;
struct arg_end  *helpend;

int PD_VERBOSE = 0;
int VERBOSE;

int main(int argc,char *argv[])
{
  int            nerrors;
   
  void *argtable[] = 
    {
     verbose  = arg_lit0("v","verbose","print debugging information"),
     file     = arg_file1(NULL,"file","<file>","pdstor to split into pieces"),
     pieces   = arg_int0("n","pieces","<num>","number of pieces"),
     size     = arg_int0("s","size","<num>","number of pdcodes in each file"),
     help     = arg_lit0(NULL,"help","display help message"),
     end      = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};

  printf("splitpd (%s)\n",PACKAGE_STRING);

 /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("splitpd: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"splitpd");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("splitpd divides a pdstor into smaller pdstors for parallel processing\n"
	     "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  if (verbose->count > 0) {

    VERBOSE = 50;

  } else { 

    VERBOSE = 0;

  }

  if (pieces->count > 0 && size->count > 0) {
    
    printf("splitpd divides a pdstor into smaller pdstors for parallel processing\n"
	   "can specify either the number of pieces to divide into _or_ the size\n"
	   "of each piece (but not both).\n\n"
	   "usage: \n\n");
    arg_print_glossary(stdout, argtable," %-25s %s\n");
    exit(0);

  }
    
  clock_t start, end;
  double cpu_time_used = 0.0;
  start = clock();

  /* First, we try to open the file. */
  
  FILE *infile;

  infile = fopen(file->filename[0],"r");

  if (infile == NULL) {

    printf("splitpd loads a pdstor and splits it into pieces for parallel processing");

     printf("Couldn't open \n"
	    "\t  file = %s, \n"
	    "usage: \n\n",file->filename[0]);
     arg_print_glossary(stdout, argtable," %-25s %s\n");
     exit(0);

  }
  
  pd_stor_t *pdstor = pd_read_pdstor(infile,NONE); /* Read without checking */

  if (pdstor == NULL) {

    printf("splitpd: Couldn't read pdstor from afile = %s",file->filename[0]);
    exit(1);

  }

  fclose(infile);

  /* Now we have parsed the arguments and are ready to work. */
  /* First, we parse the header lines from the afile... */

  unsigned int elts;
  elts = pd_stor_nelts(pdstor);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("splitpd: %4.3f seconds spent loading file of %d pd codes\n",cpu_time_used,elts);

  unsigned int subpdsize;
  unsigned int subpdcount;
  
  if (pieces->count > 0) {

    subpdsize = ceil((float)(elts)/(float)(pieces->ival[0]));
    subpdcount = ceil((float)(elts)/(float)(subpdsize));
    
  } else if (size->count > 0) {

    subpdsize = size->ival[0];
    subpdcount = ceil((float)(elts)/(float)(subpdsize));

  } else {

     printf("splitpd divides a pdstor into smaller pdstors for parallel processing\n"
	   "must specify either the number of pieces to divide into _or_ the size\n"
	   "of each piece (but not both).\n\n"
	   "usage: \n\n");
    arg_print_glossary(stdout, argtable," %-25s %s\n");
    exit(0);
  }

  printf("splitpd: dividing %d elt pdstor into ~%d pieces of size ~%d.\n",
	 elts,subpdcount,subpdsize);

  start = clock();
    
  pd_code_t *pd;
  unsigned int processed = 1;
  int percent_complete = 0;
  unsigned int subpdnum = 0;
  unsigned int subpdelt = 0;
  pd_stor_t *subpd;

  /* Now we need to delete the trailing extension from the input filename. */

  char *basename = strdup(file->basename[0]);
  char *extension_start = strstr(basename,".pdstor");
  *extension_start = 0;
 
  for(pd = pd_stor_firstelt(pdstor),
	subpd = pd_new_pdstor();
      pd != NULL;pd = pd_stor_nextelt(pdstor),
      processed++) {

    pd_addto_pdstor(subpd,pd,NONE);
    pd_code_free(&pd);
    subpdelt++;

    if (subpdelt == subpdsize) {

      char fname[256];
      sprintf(fname,"%s-%05d.pdstor",basename,subpdnum);
      FILE *outfile = fopen(fname,"w");
      pd_write_pdstor(outfile,subpd);
      pd_free_pdstor(&subpd);
      fclose(outfile);

      subpdnum++;
      subpdelt = 0;
      subpd = pd_new_pdstor();

      /* We now report progress. */
      
      if (VERBOSE < 10) { 
	
	if ((int)(floor(100.0*(double)(processed)/(double)(elts)))
	    != percent_complete) { 
	  
	  percent_complete = (int)(floor(100.0*(double)(processed)/(double)(elts)));
	  printf("splitpd: %d/%d (%d%%) of total pdcodes, written %d pieces\r",
		 processed,elts,percent_complete,subpdnum);
	  fflush(stdout);
	  
	}
	
      }
    }
  }

  if (pd_stor_nelts(subpd) != 0) {  /* Are there leftovers? */
    
    char fname[256];
    sprintf(fname,"%s-%05d.pdstor",basename,subpdnum);
    FILE *outfile = fopen(fname,"w");
    pd_write_pdstor(outfile,subpd);
    pd_free_pdstor(&subpd);
    fclose(outfile);
    subpdnum++;
    
  }
    
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

  printf("                                                                 \n"
	 "splitpd: %d element pdstor split into %d pdstors \n"
	 "         in %-4.5f seconds.\n",
	 elts,subpdnum,cpu_time_used);
  fflush(stdout);
      
  pd_free_pdstor(&pdstor);
    
  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
  free(helpend);

  exit(0);

}
