/*

  This program generates orientations (but not crossing information) for 
  shadows in a pdstor file, and outputs the results (with counts) to a 
  new pdstor file.

*/

#include <config.h>

#include <plCurve.c>
#include <plcTopology.h>
#include <pd_multidx.h>
#include <pd_cyclic.h>
#include <pd_isomorphisms.c>
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

#ifdef HAVE_PYTHON
  #include <Python.h>
  #include <pdcode/diagram_api.h>
#endif

// Turn asserts ON.
#define DEBUG 1


/* Global variables live here. */

struct arg_file *infile;  //
struct arg_file *outfile; // optional outfile override

struct arg_lit  *verbose;
struct arg_lit  *help;
struct arg_lit  *quiet;
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

      infile = arg_filen(NULL,NULL,"<filename>",1,1000,"text format (pdstor) pd_code file"),
      outfile = arg_filen("o","outfile","<filename>",0,1,"filename for (single) output file in ccode format"),
      verbose = arg_lit0(NULL,"verbose","print debugging information"),
      quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"),
      help = arg_lit0(NULL,"help","display help message"),
      end = arg_end(20)
    };

  void *helptable[] = {help,helpend = arg_end(20)};
  void *helpendtable[] = {helpend};

  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("orientedshadows: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */

    nerrors = arg_parse(argc,argv,helptable);

    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      fprintf(stderr,"orientedshadows compiled " __DATE__ " " __TIME__ "\n");
      arg_print_errors(stdout,end,"orientedshadows");

      printf("usage\n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(1);

    } else {  /* The help table matched, which means we
		 asked for help or gave nothing */

      fprintf(stderr,"orientedshadows compiled " __DATE__ " " __TIME__ "\n");
      printf("orientedshadows converts a file of pd_codes in the pdstor text file format\n"
	     "to pd codes stored _with_ orientation data (per component), but \n"
	     "_without_ crossing data. The output is written to a new pdstor file.\n"
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

  printf("orientedshadows operating on %d pdstor files\n",infile->count);

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
      ofname = mangle(infile->basename[i],".pdstor","-oriented.pdstor");
      printf("opening output file %s...",ofname);
      out = fopen(ofname,"w");

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

    printf("expanding and generating new pdstor...\n");
    int j;
    int codes_added = 0;
    pd_stor_t *newpdstor = pd_new_pdstor();
    int by_components[10] = {0,0,0,0,0,0,0,0,0,0};

    for(j=0;!feof(in);j++) {

      pd_code_t *inpd;
      inpd = pd_read(in);

      assert(inpd != NULL);

      if (!pd_ok(inpd)) {

	  pd_printf("pdcode read from file does not pass pd_ok\n",inpd);
	  exit(1);

      }

      by_components[inpd->ncomps]++;

      int i;

      for(i=0;i<inpd->ncross;i++) {

	if (inpd->cross[i].sign != PD_UNSET_ORIENTATION) {

	  if (verbose->count > 0) {

	    pd_printf("warn: input pd %PD HAS crossing information (will null)",inpd);
	  }
	  inpd->cross[i].sign = PD_UNSET_ORIENTATION;

	}

      }
      
      pd_or_t *component_orientations;
      component_orientations = calloc(inpd->ncomps,sizeof(pd_or_t));
      for(i=0;i<inpd->ncomps;i++) {
	component_orientations[i] = PD_POS_ORIENTATION;
      }

      pd_multidx_t *orientation_idx;
      pd_idx_t*    *comp_orientations;
      pd_idx_t      /*one = 1,*/ two = 2;
      comp_orientations = calloc(inpd->ncomps,sizeof(pd_idx_t *));
      
      for(i=0;i<inpd->ncomps;i++) { comp_orientations[i] = &two; }
	
      orientation_idx = pd_new_multidx(inpd->ncomps,
				       (void **)(comp_orientations),
				       cyclic_ops);
      
      unsigned int norientations = pd_multidx_nvals(orientation_idx);
      unsigned int orcount;

      for(orcount=0;
	  orcount < norientations;
	  orcount++,pd_increment_multidx(orientation_idx)) {
	
	pd_code_t *working_pd = pd_copy(inpd);

	for(i=0;i<working_pd->ncomps;i++) {
	  if (orientation_idx->i[i] == 1) {
	    pd_reorient_component(working_pd,i,PD_NEG_ORIENTATION);
	    component_orientations[i] = PD_NEG_ORIENTATION;
	  } else {
	    component_orientations[i] = PD_POS_ORIENTATION;
	  }
	}
	
	pd_addto_pdstor(newpdstor,working_pd,ISOMORPHISM);
	pd_code_free(&working_pd);
	
	codes_added++;
	
      }
      
      pd_free_multidx(&orientation_idx);
      free(comp_orientations);
      free(component_orientations);
      
      pd_code_free(&inpd);

    }

    unsigned int pds_nhashes, pds_nelts;
    pd_stor_stats(newpdstor,&pds_nhashes,&pds_nelts);
    
    printf("done (added %d oriented shadows/%d unique).\n",
	   codes_added,pds_nelts);

    printf("writing %s...",ofname);
    pd_write_pdstor(out,newpdstor);
    printf("done\n");

    printf("closing output file %s...",ofname);
    fclose(out);

    pd_free_pdstor(&newpdstor);

    double maxtotal = 0;
    
    printf("\n"
	   "Report (ISOMORPHISM)          \n"
	   "--------------------------------\n");
    for(i=1;i<10;i++) {
      printf("%5d %2d-component link shadows\n",by_components[i],i);
      maxtotal += pow(2.0,i)*by_components[i];
    }
    printf("-----------------------------------\n\n");
    printf("total number of orientations: %g\n",maxtotal);
    
  }

  printf("done\n");

  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
  arg_freetable(helpendtable,sizeof(helpendtable)/sizeof(helpendtable[0]));

  printf("done\n");

  exit(0);
}
