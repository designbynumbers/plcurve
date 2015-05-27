/* 

   computepd.c : reads a pdstor, computes some function, and writes results to file.

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
#include"pd_invariants.h"

#include"ordie.h"
#include"argtable2.h"
/* We use a local, compiled in copy of argtable to avoid having to
   depend on system utilities. */

struct arg_lit  *verbose;
struct arg_file *file;
struct arg_str  *tocompute;

struct arg_rem  *ncomps;
struct arg_rem  *uid;
struct arg_rem  *ncross;
struct arg_rem  *hash;
struct arg_rem  *interlacements;

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
     file    = arg_file1(NULL,"file","<file>","pdstor to read"),
     tocompute = arg_strn("c","compute","<quantity>",1,20,"thing to compute"),
     hash = arg_rem(NULL, "  hash : 32 character hash for pdcode"),
     uid = arg_rem(NULL,"  uid : unsigned int unique id (among codes with same hash)"),
     ncomps = arg_rem(NULL,"  ncomps : # of components"),
     ncross = arg_rem(NULL,"  ncross : # of crossings"),
     interlacements = arg_rem(NULL, "  interlacements : # of interlaced crossing pairs (by component)"),
     outfile = arg_file1("o","output","<file>","output csv file"),
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};

  printf("computepd (%s)\n",PACKAGE_STRING);

 /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("computepd: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"computepd");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("computepd loads a pdstor, computes something for each pd, \n"
	     "and writes the output to a file\n"
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

    printf("computepd loads a pdstor, computes something for each pd, \n"
	   "and writes the output to a file\n");

     printf("Couldn't open \n"
	    "\t  file = %s, \n"
	    "usage: \n\n",file->filename[0]);
     arg_print_glossary(stdout, argtable," %-25s %s\n");
     exit(0);

  }
  
  pd_stor_t *pdstor = pd_read_pdstor(infile,NONE); /* Read without checking */

  if (pdstor == NULL) {

    printf("computepd: Couldn't read pdstor from file = %s",file->filename[0]);
    exit(1);

  }

  fclose(infile);

  FILE *of = fopen(outfile->filename[0],"w");

  if (of == NULL) {

    printf("computepd: Couldn't open outfile = %s",outfile->filename[0]);
    exit(1);

  }
  
  /* Now we have parsed the arguments and are ready to work. */
  /* First, we parse the header lines from the afile... */

  unsigned int elts;

  elts = pd_stor_nelts(pdstor);
  
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("computepd: %4.3f seconds spent loading file of %d pd codes\n",cpu_time_used,elts);
  start = clock();
    
  pd_code_t *pd;
  unsigned int processed = 0;
  int percent_complete;
  
  for(pd = pd_stor_firstelt(pdstor);pd != NULL;pd = pd_stor_nextelt(pdstor),
      processed++) {

    /* Now we try to parse the tocompute options */

    int ntc = 0;

    fprintf(of,"%d",processed);
	    
    for(ntc=0;ntc<tocompute->count;ntc++) {

      if (!strcmp(tocompute->sval[ntc],"hash")) {

	fprintf(of,",%s ",pd->hash);

      } else if (!strcmp(tocompute->sval[ntc],"uid")) {

	fprintf(of,",%lu ",pd->uid);

      } else if (!strcmp(tocompute->sval[ntc],"ncomps")) {

	fprintf(of,",%d ",(int)(pd->ncomps));

      } else if (!strcmp(tocompute->sval[ntc],"ncross")) {

	fprintf(of,",%d ",(int)(pd->ncross));

      } else if (!strcmp(tocompute->sval[ntc],"interlacements")) {

	int *interlacements = pd_interlaced_crossings(pd);

	int i;
	for(i=0;i<pd->ncomps;i++) { 
	
	  fprintf(of,",%d ",interlacements[i]);

	}
	
	free(interlacements);

      } else if (!strcmp(tocompute->sval[ntc],"unsigned-interlacements")) {

	unsigned int *interlacements = pd_interlaced_crossings_unsigned(pd);

	int i;
	for(i=0;i<pd->ncomps;i++) { 
	
	  fprintf(of,",%u ",interlacements[i]);

	}
	
	free(interlacements);

      }

    }

    fprintf(of,"\n");

    /* We now report progress. */

    if (VERBOSE < 10) { 

      if ((int)(floor(100.0*(double)(processed)/(double)(elts)))
	  != percent_complete) { 
	
	percent_complete = (int)(floor(100.0*(double)(processed)/(double)(elts)));
	printf("computepd: %d/%d (%d%%) of pdcodes\r",
	       processed,elts,percent_complete);
	fflush(stdout);
	
      }
	
    }

  }
    
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

  printf("                                                                 \n"
	 "computepd: %d pdcodes processed \n"
	 "           in %-4.5f seconds.\n",
	 elts,cpu_time_used);
  fflush(stdout);

  fclose(of);
  pd_free_pdstor(&pdstor);
    
  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
  free(helpend);

  exit(0);

}
