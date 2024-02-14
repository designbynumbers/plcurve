/* 

   connectsum.c : uses the pdcode functions to connect sum all
                  elements of a pdstor with all elements of another
                  pdstor (in all ways), adding the results to a third
                  pdstor. This is a "round" in the connect sum
                  generation of all diagrams.

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
struct arg_file *Afile;
struct arg_file *Bfile;
struct arg_file *outfile;
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
     Afile   = arg_file1("a","afile","<file>","pdstor of A summands"),
     Bfile   = arg_file1("b","bfile","<file>","pdstor of B summands"),
     outfile = arg_file1("o","output","<file>","output filename"),
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help, helpEnd = arg_end(20)};

  printf("connectsum (%s)\n",PACKAGE_STRING);

 /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("connectsum: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"connectsum");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("connectsum takes all possible connect sums of pdcodes in \n"
	     "two pdstors in files A and B, and enters resulting pdcodes \n"
	     "(up to isomorphism) in file -o <file>. \n"
	     "All three arguments are required. \n"
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

  /* First, we try to open the three files. We're going to have to
     read the pdstor in C into memory to work with it, but the pdstors
     in A and B can be left largely on disk to save memory. We will have
     a problem if the Afile and bfile are the same, but we're going to 
     deal with that problem at the level of scripting. */

  FILE *afile, *bfile;

  afile = fopen(Afile->filename[0],"r");
  bfile = fopen(Bfile->filename[0],"r");

  if (afile == NULL || bfile == NULL) {

    printf("connectsum takes all possible connect sums of pdcodes in \n"
	   "two pdstors in files A and B, and enters resulting pdcodes \n"
	   "(up to isomorphism) in file -o <file>. \n"
	   "All three arguments are required. \n");

     printf("Couldn't open one of\n"
	    "\t afile = %s, \n"
	    "\t bfile = %s, \n"
	    "usage: \n\n",Afile->filename[0],Bfile->filename[0]);
     arg_print_glossary(stdout, argtable," %-25s %s\n");
     exit(0);

  }
  
  pd_stor_t *cstor = pd_new_pdstor();
  
  pd_stor_t *astor = pd_read_pdstor(afile,NONE); /* Read without checking */

  if (astor == NULL) {

    printf("connectsum: Couldn't read pdstor from afile = %s",Afile->filename[0]);
    exit(1);

  }

  fclose(afile);

  pd_stor_t *bstor = pd_read_pdstor(bfile,NONE);  /* Read without checking */

  if (bstor == NULL) {

    printf("connectsum: Couldn't read pdstor from bfile = %s",Bfile->filename[0]);
    exit(1);

  }

  fclose(bfile);

  /* Now we have parsed the arguments and are ready to work. */
  /* First, we parse the header lines from the afile... */

  unsigned int aelts, belts;

  aelts = pd_stor_nelts(astor);
  belts = pd_stor_nelts(bstor);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("connectsum: %4.3f seconds spent preparing for connect sum run\n",cpu_time_used);
  start = clock();
    
  pd_code_t *pdA, *pdB;
  unsigned int aprocessed = 1;
  int percent_complete;
  unsigned int total_cs = 0;

  for(pdA = pd_stor_firstelt(astor);pdA != NULL;pdA = pd_stor_nextelt(astor),
      aprocessed++) {

    for(pdB = pd_stor_firstelt(bstor);pdB != NULL;pdB = pd_stor_nextelt(bstor)) {

      /* Since we are only operating up to isomorphism anyway, we null the crossing
	 signs (if any) on pdA and pdB right now. */

      pd_idx_t cr;
      for(cr=0;cr<pdA->ncross;cr++) { pdA->cross[cr].sign = PD_UNSET_ORIENTATION; }
      for(cr=0;cr<pdB->ncross;cr++) { pdB->cross[cr].sign = PD_UNSET_ORIENTATION; }

      /* Iterate over all orientations of components of pdA and pdB */
      /* To do this, we need to make a multidx with orientation_ops */

      pd_multidx_t *a_orientation_iter = pd_new_multidx(pdA->ncomps,NULL,orientation_ops);
      unsigned int amax = pd_multidx_nvals(a_orientation_iter);
      unsigned int acount;

      /* Now we're going to check up to isomorphism, which allows a single, global
	 orientation reversal (on all components). This means that we can assume 
         w.l.o.g. that the first component of a has POSITIVE orientation. */ 
      
      for(acount=0;acount<amax;acount++,pd_increment_multidx(a_orientation_iter)) {

	if (((pd_orientation_t *)(a_orientation_iter->obj[0]))->orient == PD_POS_ORIENTATION) {
	  
	  pd_multidx_t *b_orientation_iter = pd_new_multidx(pdB->ncomps,NULL,orientation_ops);
	  unsigned int bmax = pd_multidx_nvals(b_orientation_iter);
	  unsigned int bcount;
	  
	  pd_code_t *workingA = pd_copy(pdA);
	  pd_idx_t   compA;
	  
	  for(compA=0;compA<pdA->ncomps;compA++) {
	    
	    pd_reorient_component(workingA,compA,
				  ((pd_orientation_t *)(a_orientation_iter->obj[compA]))->orient) ;
	    
	  }
	  
	  for(bcount=0;bcount<bmax;bcount++,pd_increment_multidx(b_orientation_iter)) {
	    
	    pd_code_t *workingB = pd_copy(pdB);
	    pd_idx_t   compB;
	    
	    for(compB=0;compB<pdB->ncomps;compB++) {
	      
	      pd_reorient_component(workingB,compB,
				    ((pd_orientation_t *)(b_orientation_iter->obj[compB]))->orient) ;
	      
	    }
	    
	    pd_idx_t edgeA, edgeB;
	    
	    /* Iterate over all edge pairs. This requires 
	       a little infrastructure, I think.*/
	    
	    for(edgeA = 0; edgeA < workingA->nedges; edgeA++) {
	      
	      for(edgeB = 0; edgeB < workingB->nedges; edgeB++) {
		
		pd_code_t *pdC = pd_connect_sum(workingA,edgeA,workingB,edgeB);
		pd_addto_pdstor(cstor,pdC,ISOMORPHISM);
	      total_cs++;
	      
	      free(pdC);
	      
	      }
	      
	    }

	    /* Now do housekeeping */

	    pd_code_free(&workingB);

	  }

	  pd_free_multidx(&b_orientation_iter);
	  pd_code_free(&workingA);

	}

      }

      pd_free_multidx(&a_orientation_iter);

    }

    /* We now report progress. */

    if (VERBOSE < 10) { 

      if ((int)(floor(100.0*(double)(aprocessed)/(double)(aelts)))
	  != percent_complete) { 
	
	percent_complete = (int)(floor(100.0*(double)(aprocessed)/(double)(aelts)));
	printf("connectsum: %d/%d (%d%%) of A pdcodes/C contains %d elts, %u cs\r",
	       aprocessed,aelts,percent_complete,pd_stor_nelts(cstor),total_cs);
	fflush(stdout);
	
      }
	
    }

  }
    
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

  

  printf("\r"
	 "connectsum: %d connect sums (%d isomorphism classes) \n"
	 "            computed in %-4.5f seconds (%g sec/cs).\n",
	 total_cs,pd_stor_nelts(cstor),cpu_time_used,cpu_time_used/(double)(total_cs));
  fflush(stdout);

  FILE *out;

  printf("connectsum: writing connect sums to %s...",outfile->filename[0]);
  out = fopen(outfile->filename[0],"w");
  assert(out != NULL);
  
  pd_write_pdstor(out,cstor);
  fclose(out);
  printf("done\n");
      
  pd_free_pdstor(&astor);
  pd_free_pdstor(&bstor);
  pd_free_pdstor(&cstor);
    
  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
  free(helpEnd);

  exit(0);

}
