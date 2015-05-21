/* 

   filterpd.c : reads a pdstor, applies various filters, and writes a new one.

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
struct arg_int  *ncomps;
struct arg_file *outfile;
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
     verbose = arg_lit0("v","verbose","print debugging information"),
     file    = arg_file1(NULL,"file","<file>","pdstor to filter"),
     ncomps  = arg_int0(NULL,"ncomps","<num>","number of components"),
     outfile = arg_file1("o","output","<file>","output pdstor (will be overwritten)"),
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};

  printf("filterpd (%s)\n",PACKAGE_STRING);

 /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("filterpd: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"filterpd");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("filterpd take a pdstor file and applies various filters\n"
	     "to the entries\n"
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

  clock_t start, end;
  double cpu_time_used = 0.0;
  start = clock();

  /* First, we try to open the file. */
  
  FILE *infile;

  infile = fopen(file->filename[0],"r");

  if (infile == NULL) {

    printf("filterpd loads a pdstor, extracts some of the contents \n"
	   "and writes the pdstors which match various criteria to a file\n");

     printf("Couldn't open \n"
	    "\t  file = %s, \n"
	    "usage: \n\n",file->filename[0]);
     arg_print_glossary(stdout, argtable," %-25s %s\n");
     exit(0);

  }
  
  pd_stor_t *outstor = pd_new_pdstor();
  pd_stor_t *pdstor = pd_read_pdstor(infile,NONE); /* Read without checking */

  if (pdstor == NULL) {

    printf("filterpd: Couldn't read pdstor from afile = %s",file->filename[0]);
    exit(1);

  }

  fclose(infile);

  /* Now we have parsed the arguments and are ready to work. */
  /* First, we parse the header lines from the afile... */

  unsigned int elts;

  elts = pd_stor_nelts(pdstor);
  
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("filterpd: %4.3f seconds spent loading file of %d pd codes\n",cpu_time_used,elts);
  start = clock();
    
  pd_code_t *pd;
  unsigned int processed = 1;
  int percent_complete;
  bool keep;
 
  for(pd = pd_stor_firstelt(pdstor);pd != NULL;pd = pd_stor_nextelt(pdstor),
      processed++) {

    keep = true;
    
    if (ncomps->count > 0) {

      if (pd->ncomps != ncomps->ival[0]) { keep = false; }

    }

    if (keep) {

      pd_addto_pdstor(outstor,pd,NONE);
      
    }
	      
    free(pd);
	      
    /* We now report progress. */

    if (VERBOSE < 10) { 

      if ((int)(floor(100.0*(double)(processed)/(double)(elts)))
	  != percent_complete) { 
	
	percent_complete = (int)(floor(100.0*(double)(processed)/(double)(elts)));
	printf("filterpd: %d/%d (%d%%) of A pdcodes/out contains %d elts\r",
	       processed,elts,percent_complete,pd_stor_nelts(outstor));
	fflush(stdout);
	
      }
	
    }

  }
    
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

  printf("                                                                 \n"
	 "filterpd: %d pdcodes filtered to %d pdcodes \n"
	 "          computed in %-4.5f seconds.\n",
	 elts,pd_stor_nelts(outstor),cpu_time_used);
  fflush(stdout);

  FILE *out;

  printf("filterpd: writing filtered data to %s...",outfile->filename[0]);
  out = fopen(outfile->filename[0],"w");
  assert(out != NULL);
  
  pd_write_pdstor(out,outstor);
  fclose(out);
  printf("done\n");
      
  pd_free_pdstor(&pdstor);
  pd_free_pdstor(&outstor);
    
  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
  free(helpend);

  exit(0);

}
