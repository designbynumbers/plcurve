/* 

   joinpd.c : adds elements of a pdstor to another, replacing 
              the second pdstor with the new information.

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

#include"ordie.h"
#include"argtable2.h"
/* We use a local, compiled in copy of argtable to avoid having to
   depend on system utilities. */

struct arg_lit  *verbose;
struct arg_file *files;
struct arg_file *outfile;
struct arg_lit  *help;
struct arg_end  *end;
struct arg_end  *helpend;

int PD_VERBOSE;
int VERBOSE;

int main(int argc,char *argv[])
{
  int            nerrors;
   
  void *argtable[] = 
    {
     verbose = arg_lit0("v","verbose","print debugging information"),
     files   = arg_filen(NULL,NULL,"<files>",1,5000,"pdstor files to combine"),
     outfile = arg_file1("o","output","<output filename>","output pdstor name"),
     nochecking = arg_lit0(NULL,"no-checking","no isomorphism checking"),
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};

  printf("joinpd (%s)\n",PACKAGE_STRING);

 /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("joinpd: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"joinpd");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("joinpd adds contents of one pdstor to another\n"
	     "with isomorphism checking\n"
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

  /* First, we try to open the files. */

  FILE *firstfile, *bfile;

  firstfile = fopen(files->filename[0],"r");

  if (firstfile == NULL) {
    
    printf("joinpd adds contents of one pdstor to another\n"
	   "with isomorphism checking\n"
	   "\n"
	   "Couldn't open first file %s \n"	   
	   "usage: \n\n",files->filename[0]);
    arg_print_glossary(stdout, argtable," %-25s %s\n");
    exit(0);
    
  }

  printf("joinpd: Loading %s...",files->filename[0]);
  fflush(stdout);

  pd_stor_t *astor = pd_read_pdstor(firstfile,ISOMORPHISM);

  if (astor == NULL) {

    printf("joinpd: Couldn't read pdstor from first file = %s",files->filename[0]);
    exit(1);

  }

  fclose(firstfile);

  unsigned int aelts;
  aelts = pd_stor_nelts(astor);

  printf("done (%d elts)\n",aelts);

  int i;
  for(i=1;i<files->count;i++) {

    bfile = fopen(files->filename[i],"r");
    
    if (bfile == NULL) {
      
      printf("joinpd adds contents of one pdstor to another\n"
	     "with isomorphism checking\n"
	     "\n"
	     "Couldn't open file %d with name %s\n"
	     "usage: \n\n",i,files->filename[i]);
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
    pd_stor_t *bstor = pd_read_pdstor(bfile,NONE);
    
    if (bstor == NULL) {
      
      printf("joinpd: Couldn't read pdstor from file = %s",files->filename[i]);
      exit(1);
      
    }
    
    fclose(bfile);
    
    /* Now we have parsed the arguments and are ready to work. */
    /* First, we parse the header lines from the afile... */
    
    unsigned int  belts;
    belts = pd_stor_nelts(bstor);
    
    printf("joinpd: preparing to add %s (%d elts) to main pdstor (%d elts)\n",
	   files->filename[i],belts,pd_stor_nelts(astor));
    clock_t start, end;
    double cpu_time_used = 0.0;
    start = clock(); 
    
    pd_code_t *pdA;
    int percent_complete;
    int processed = 0;
    
    for(pdA = pd_stor_firstelt(bstor);pdA != NULL;pdA = pd_stor_nextelt(bstor),
	  processed++) {
      
      pd_addto_pdstor(astor,pdA,ISOMORPHISM);
      
      /* We now report progress. */
      
      if (VERBOSE < 10) { 
	
	if ((int)(floor(100.0*(double)(processed)/(double)(aelts)))
	    != percent_complete) { 
	  
	  percent_complete = (int)(floor(100.0*(double)(processed)/(double)(aelts)));
	  printf("joinpd: %d/%d (%d%%) of pdcodes/main contains %d elts\r",
		 processed,belts,percent_complete,pd_stor_nelts(astor)); 
	  
	}
	
      }
      
    }
    
    printf("joinpd: %d/%d (%d%%) of A pdcodes/B contains %d elts\r",
	   processed,belts,100,pd_stor_nelts(astor));    
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    printf("\njoinpd: %d pd codes added (%d total now) in %-4.5f seconds.\n",
	   pd_stor_nelts(bstor),pd_stor_nelts(astor),cpu_time_used);
    fflush(stdout);

    pd_free_pdstor(&bstor);
    
  }
  
  FILE *out;
  
  printf("joinpd: writing new pdstor to %s...",outfile->filename[0]);
  out = fopen(outfile->filename[0],"w");
  pd_write_pdstor(out,astor);
  fclose(out);

  printf("done\n");
  
  pd_free_pdstor(&astor);
      
  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
  free(helpend);

  exit(0);

}
