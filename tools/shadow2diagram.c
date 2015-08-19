/*

  This program spits out the pdstor files in Ken Millett's ccode format
  in order to allow checking against Eric R's knot identification code.

*/

#include <plCurve.c>
#include <plcTopology.h>
#include <pd_multidx.h>
#include <pd_cyclic.h>
#include <pd_perm.h>
#include <pd_isomorphisms.h>
#include <pd_storage.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include <argtable2.h>
#include <assert.h>
#include <string.h>
#include <sys/stat.h>

#include <Python.h>
#include <pdcode/diagram_api.h>

// Turn asserts ON.
#define DEBUG 1

/* Global variables live here. */

struct arg_file *infile;  //
struct arg_file *outfile; // optional outfile override

struct arg_lit  *ccodes;
struct arg_lit  *verbose;
struct arg_lit  *help;
struct arg_lit  *quiet;
struct arg_lit  *KnotTheory;
struct arg_lit  *countonly;
struct arg_lit  *simplify;
struct arg_end  *end;
struct arg_end  *helpend;

char *mangle(const char *filename,const char *oldextension,const char *newextension);
void  nmangle(char *newname,int nnsize,
	      const char *filename,const char *oldextension,const char *newextension);


FILE *fmangle(const char *filename,const char *oldextension,const char *newextension)
{
  FILE *outfile;
  char *name;

  name = mangle(filename,oldextension,newextension);

  outfile = fopen(name,"w");

  if (outfile == NULL) {

    fprintf(stderr,"fmangle: Could not open filename %s.\n",
	    name);

    exit(1);

  }

  free(name);
  return outfile;
}

char *mangle(const char *filename,const char *oldextension,const char *newextension)
/* Allocate space, then perform an "nmangle". */
{
  char *newname;
  int nnsize;

  nnsize = strlen(filename) + strlen(newextension) + 10;
  newname = calloc(nnsize,sizeof(char));

  if (newname == NULL) {

    fprintf(stderr,"mangle: Couldn't allocate string of size %d.\n",nnsize);
    exit(1);

  }

  nmangle(newname,nnsize,filename,oldextension,newextension);

  return newname;
}

void  nmangle(char *newname,int nnsize,
	      const char *filename,const char *oldextension,const char *newextension)

/* Replace the (terminating) string "oldextension" with "newextension" if present in
"filename". Return the results in "newname". */
{

  if (nnsize < strlen(filename) + strlen(newextension) + 2) {

    fprintf(stderr,"nmangle: Buffer newname is not long enough to add extension %s to filename %s.\n",
	    newextension,filename);
    exit(1);

  }

  /* We have already checked that we have enough space, but we're cautious. */

  strncpy(newname,filename,nnsize);

  if (strstr(newname,oldextension) != NULL) {

    strcpy(strstr(newname,oldextension),newextension);

  } else {

    strncat(newname,newextension,nnsize);

  }

}

char *pdcode_to_ccode(pd_code_t *pd);

int main(int argc,char *argv[]) {

  int nerrors;

  void *argtable[] =
    {

      infile = arg_file1(NULL,NULL,"<filename>","text format (pdstor) pd_code file"),
      outfile = arg_filen("o","outfile","<filename>",0,1,"filename for (single) output file"),      
      verbose = arg_lit0(NULL,"verbose","print debugging information"),
      quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"),
      KnotTheory = arg_lit0("K","KnotTheory","print pd codes in the style of knottheory"),
      countonly = arg_lit0(NULL,"count-only","count diagrams, but don't write them to disk"),
      ccodes = arg_lit0("C","ccodes","print pd codes in Millett/Ewing ccode format"),
      help = arg_lit0(NULL,"help","display help message"),
      end = arg_end(20)
    };

  void *helptable[] = {help,helpend = arg_end(20)};
  void *helpendtable[] = {helpend};

  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("shadow2diagram: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */

    nerrors = arg_parse(argc,argv,helptable);

    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      fprintf(stderr,"shadow2diagram compiled " __DATE__ " " __TIME__ "\n");
      arg_print_errors(stdout,end,"shadow2diagram");

      printf("usage\n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(1);

    } else {  /* The help table matched, which means we
		 asked for help or gave nothing */

      fprintf(stderr,"shadow2diagram compiled " __DATE__ " " __TIME__ "\n");
      printf("shadow2diagram converts pdcodes without crossing information into\n"
	     "diagrams with all possible crossing signs and orientations. The\n"
	     "diagrams are identified up to diagram-isotopy.\n"
	     "\n"
	     "By default, the output is a new pdstor file.\n"
	     "\n"
	     "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }

  }

  /******************** Now we deal with filenames *******************/

  char *ofname = NULL;

  if (outfile->count != 0) { /* We'll have to make something up */

    ofname = calloc(4096,sizeof(char));
    strncpy(ofname,outfile->filename[0],4096);

  }

  /* Now we actually do the work */

  printf("shadow2diagram crossing sign/orientation generator for pdstor files\n");
  printf("input is %d pdstor files\n",infile->count);

  FILE *out;

  if (outfile->count != 0) {

    printf("opening output file %s...",ofname);
    out = fopen(ofname,"w");

    if (out == NULL) {

      printf("fail\n");
      exit(1);

    }

    printf("done\n");

  }

  /* Otherwise, we're going to have a separate outfile for each infile. */

  int i;
  for(i=0;i<infile->count;i++) {

    if (outfile->count == 0) { /* We don't have a file open yet. */

      if (ofname != NULL) { free(ofname); }

      if (KnotTheory->count > 0) { 
	ofname = mangle(infile->basename[i],".pdstor","-diagrams.KnotTheory");
      } else if (ccodes->count > 0) {
	ofname = mangle(infile->basename[i],".pdstor","-diagrams.ccodes");
      }	else {
	ofname = mangle(infile->basename[i],".pdstor","-diagrams.pdstor");
      }	
	
      printf("opening output file %s...",ofname);
      out = fopen(ofname,"rw");

      if (out == NULL) {

	printf("fail\n");
	exit(1);

      }

      printf("done\n");

    }

    printf("opening file %s...",infile->filename[i]);
    FILE *in;
    in = fopen(infile->filename[i],"r");
    if (in == NULL) {
      printf("fail.\n");
      exit(1);
    }
    printf("done\n");

    printf("parsing header...");
    int nelts_claimed,nelts_actual,nhashes;
    if (fscanf(in,
	       "pdstor \n"
	       "nelts %d/%d (claimed/actual) nhashes %d\n\n",
	       &nelts_claimed,&nelts_actual,&nhashes) != 3) {

      printf("fail (couldn't read header)\n");
      exit(1);

    }

    if (nelts_claimed != nelts_actual) {

      printf("fail (nelts claimed %d and actual %d don't match)\n",
	     nelts_claimed,nelts_actual);

      exit(1);

    }

    printf("done (%d pd codes, %d hashes)\n",nelts_claimed,nhashes);

    printf("expanding by assigning crossing/orientation information...");
    int j;
    int codes_added = 0;

    pd_start_incremental_pdstor(out);
    unsigned int total_nhashes = 0, total_nelts = 0;
    
    for(j=0;!feof(in);j++) {

      pd_stor_t *diagrams = pd_new_pdstor();
      pd_code_t *inpd;
      inpd = pd_read(in);

      assert(inpd != NULL);

      if (!pd_ok(inpd)) {

	  pd_printf("pdcode read from file does not pass pd_ok\n",inpd);
	  exit(1);

      }
   
      pd_multidx_t *orientation_idx;
      pd_idx_t*    *comp_orientations;
      pd_idx_t      two = 2;
      comp_orientations = calloc(inpd->ncomps,sizeof(pd_idx_t *));
      
      for(i=0;i<inpd->ncomps;i++) { comp_orientations[i] = &two; }
	   
      orientation_idx = pd_new_multidx(inpd->ncomps,
				       (void **)(comp_orientations),
				       cyclic_ops);
      
      unsigned int norientations = pd_multidx_nvals(orientation_idx);
      unsigned int orcount;

      /* 
	 DEBUGGING CODE! We're checking for diagrams that have an
	 automorphism group size of one. For these, EVERY assignment
	 of crossings and orientations ought to result in a unique
	 diagram added to the list. 
      */

      bool no_orientation_preserving_automorphisms = false;
      unsigned int nisos;
      pd_iso_t **isos = pd_build_isos(inpd,inpd,&nisos);
      if (nisos == 1) {
	no_orientation_preserving_automorphisms = true;
      }
      pd_free_isos(&nisos,&isos);

      bool no_orientation_reversing_automorphisms = false;
      pd_code_t *reversed_pd = pd_copy(inpd);
      for(i=0;i<inpd->ncomps;i++) {
	pd_reorient_component(reversed_pd,i,PD_NEG_ORIENTATION);
      }			      
      isos = pd_build_isos(inpd,reversed_pd,&nisos);
      if (nisos == 0) {
	no_orientation_reversing_automorphisms = true;
      }
      pd_free_isos(&nisos,&isos);
      pd_code_free(&reversed_pd);

      unsigned int nhashes_in, nelts_in;
      pd_stor_stats(diagrams,&nhashes_in,&nelts_in);
      
      for(orcount=0;
	  orcount < norientations;
	  orcount++,pd_increment_multidx(orientation_idx)) {

	pd_multidx_t *crossing_idx;
	pd_idx_t*    *comp_crossings;
	
	comp_crossings = calloc(inpd->ncross,sizeof(pd_idx_t *));
	for(i=0;i<inpd->ncross;i++) { comp_crossings[i] = &two; }

	crossing_idx = pd_new_multidx(inpd->ncross,
				      (void **)(comp_crossings),
				      cyclic_ops);
	
	unsigned int ncrossings = pd_multidx_nvals(crossing_idx);
	unsigned int crosscount;
	
	for(crosscount=0;
	    crosscount < ncrossings;
	    crosscount++,pd_increment_multidx(crossing_idx)) {
	  
	  pd_code_t *working_pd = pd_copy(inpd);
	  
	  for(i=0;i<working_pd->ncomps;i++) {
	    
	    if (orientation_idx->i[i] == 1) {
	      pd_reorient_component(working_pd,i,PD_NEG_ORIENTATION);
	    } else {
	      /* Don't reorient the component */
	    }
	    
	  }
	  
	  for(i=0;i<working_pd->ncross;i++) {
	    
	    if (crossing_idx->i[i] == 1) {
	      working_pd->cross[i].sign = PD_NEG_ORIENTATION;
	    } else {
	      working_pd->cross[i].sign = PD_POS_ORIENTATION;
	    }
	    
	  }

	  pd_addto_pdstor(diagrams,working_pd,DIAGRAM_ISOTOPY);
	  pd_code_free(&working_pd);

	  codes_added++;

	}
	
	pd_free_multidx(&crossing_idx);
	free(comp_crossings);
	
      }
      
      pd_free_multidx(&orientation_idx);
      free(comp_orientations);

      unsigned int nhashes_out, nelts_out;
      pd_stor_stats(diagrams,&nhashes_out,&nelts_out);

      if (no_orientation_preserving_automorphisms &&
	  no_orientation_reversing_automorphisms) {

	int expected_elts = (int)(pow(2.0,inpd->ncross+1));
	if (expected_elts != nelts_out - nelts_in) {

	  pd_printf("%d-crossing pd code \n"
		    "%PD\n"
		    "has no automorphisms (orientation \n"
		    "reversing OR preserving), but running \n"
		    "through crossing/orientation loop \n"
		    "added %d new diagram-isotopy types\n"
		    "instead of 2^(%d + 1) = %d.\n",inpd,
		    inpd->ncross,nelts_out-nelts_in,inpd->ncross,expected_elts);
	  exit(1);

	}

      }

      /* if (no_orientation_preserving_automorphisms && */
      /* 	  !no_orientation_reversing_automorphisms) { */

      /* 	int expected_elts = (int)(pow(2.0,inpd->ncross)); */
      /* 	if (expected_elts != nelts_out - nelts_in) { */

      /* 	  pd_printf("%d-crossing pd code \n" */
      /* 		    "%PD\n" */
      /* 		    "has only the orientation-reversing \n" */
      /* 		    "automorphism, but running \n" */
      /* 		    "through crossing/orientation loop \n" */
      /* 		    "added %d new diagram-isotopy types\n" */
      /* 		    "instead of 2^(%d) = %d.\n",inpd, */
      /* 		    inpd->ncross,nelts_out-nelts_in,inpd->ncross,expected_elts); */
      /* 	  exit(1); */

      /* 	} */

      /*}      */

      pd_code_free(&inpd);

      /* How we update depends on the flags... */

      if (KnotTheory->count > 0 || ccodes->count > 0) {
	
	pd_code_t *popcode;
	for(popcode = pd_stor_firstelt(diagrams);
	    popcode != NULL;
	    popcode = pd_stor_nextelt(diagrams)) {
	  
	  if (KnotTheory->count > 0) {
	    
	    pd_write_KnotTheory(out,popcode);
	    
	  } else if (ccodes->count > 0) {
	    
	    fprintf(out,"%s",pdcode_to_ccode(popcode));
	    
	  }
	
	  pd_code_free(&popcode);
	
	}

      } else {
      
	pd_addto_incremental_pdstor(out,diagrams,
				    &total_nhashes,
				    &total_nelts);

      }
      
      pd_free_pdstor(&diagrams);
     
      
    }
    
    printf("done\n");
    
    printf("wrote %d diagram pdcodes/%d diagram-isotopy types.\n",
	   codes_added,total_nelts);

    pd_finish_incremental_pdstor(out,total_nhashes,total_nelts);
    fclose(out);    

  }

  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
  arg_freetable(helpendtable,sizeof(helpendtable)/sizeof(helpendtable[0]));

  exit(0);
}
