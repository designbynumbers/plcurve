/* 

   joinpd.c : adds elements of a pdstor to another, replacing 
              the second pdstor with the new information.

*/

#ifdef HAVE_CONFIG_H
#include"config.h"
#endif

#include<stdio.h>
#include<stdbool.h>
#include<string.h>
#include<stdint.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<assert.h>

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
struct arg_lit  *nochecking;
struct arg_lit  *isomorphism;
struct arg_lit  *isotopy;
struct arg_rem  *mustdo;
struct arg_lit  *help;
struct arg_end  *end;
struct arg_end  *helpEnd;

int VERBOSE;

int main(int argc,char *argv[])
{
  int            nerrors;
   
  void *argtable[] = 
    {
     verbose = arg_lit0("v","verbose","print debugging information"),
     files   = arg_filen(NULL,NULL,"<files>",1,100000,"pdstor files to combine"),
     outfile = arg_file1("o","output","<output filename>","output pdstor name"),
     help = arg_lit0(NULL,"help","display help message"),
     mustdo = arg_rem("\nUser must specify one:\n"," "),
     nochecking = arg_lit0(NULL,"no-checking","no isomorphism checking"),
     isomorphism = arg_lit0(NULL,"isomorphism-checking","structural isomorphism, ignores crossings"),
     isotopy = arg_lit0(NULL,"isotopy-checking","diagram isotopy, uses crossing data"),
     end = arg_end(20)};
  
  void *helptable[] = {help, helpEnd = arg_end(20)};

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
	     "with isomorphism/isotopy/no checking\n"
	     "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  int isotypes = nochecking->count + isomorphism->count + isotopy->count;
  
  if (isotypes != 1) {

      printf("joinpd adds contents of one pdstor to another\n"
	     "with isomorphism/isotopy/no checking. Must specify\n"
	     "type of checking\n"
	     "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

  }

  pd_equivalence_t relation;
  char relstring[256];
  char relshort[256];

  if (nochecking->count == 1) {

    relation = NONE;
    sprintf(relstring,"no checking");
    sprintf(relshort,"(NONE)");

  } else if (isomorphism->count == 1) {

    relation = ISOMORPHISM;
    sprintf(relstring,"isomorphism (ignoring crossing signs)");
    sprintf(relshort,"(MORPH)");
    
  } else if (isotopy->count == 1) {

    relation = DIAGRAM_ISOTOPY;
    sprintf(relstring,"isotopy (using crossing signs)");
    sprintf(relshort,"(TOPY)");
    
  } else {

     printf("joinpd adds contents of one pdstor to another\n"
	     "with isomorphism/isotopy/no checking. Must specify\n"
	     "type of checking\n"
	     "usage: \n\n");
     arg_print_glossary(stdout, argtable," %-25s %s\n");
     exit(0);

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
	   "with checking set to %s\n"
	   "\n"
	   "Couldn't open first file %s \n"	   
	   "usage: \n\n",relstring,files->filename[0]);
    arg_print_glossary(stdout, argtable," %-25s %s\n");
    exit(0);
    
  }

  printf("joinpd: Checking set to %s\n",relstring);
  printf("joinpd: Loading %s...",files->filename[0]);
  fflush(stdout);

  pd_stor_t *astor = pd_read_pdstor(firstfile,relation);

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
      
      pd_addto_pdstor(astor,pdA,relation);
      pd_code_free(&pdA);
      
      /* We now report progress. */
      
      if (VERBOSE < 10) { 
	
	if ((int)(floor(100.0*(double)(processed)/(double)(aelts)))
	    != percent_complete) { 
	  
	  percent_complete = (int)(floor(100.0*(double)(processed)/(double)(aelts)));
	  printf("joinpd %s: %d/%d (%d%%) of pdcodes/main contains %d elts\r",
		 relshort,processed,belts,percent_complete,pd_stor_nelts(astor)); 
	  
	}
	
      }
      
    }
    
    printf("joinpd %s: %d/%d (%d%%) of A pdcodes/B contains %d elts\r",
	   relstring,processed,belts,100,pd_stor_nelts(astor));    
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    printf("\njoinpd %s: %d pd codes added (%d total now) in %-4.5f seconds.\n",
	   relstring,pd_stor_nelts(bstor),pd_stor_nelts(astor),cpu_time_used);
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
  free(helpEnd);

  exit(0);

}
