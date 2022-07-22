/*

  This program generates random polygonal knots, filtering by given knot type.

*/

#ifdef HAVE_CONFIG_H
#include"config.h"
#endif

#include <plCurve.h>
#include <plcTopology.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include <argtable2.h>
#include <assert.h>
#include <string.h>
#include <sys/stat.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//#include <Python.h>
//#include <pdcode/diagram_api.h>

// Turn asserts ON.
#define DEBUG 1
int PD_VERBOSE=0;

/* Global variables live here. */

struct arg_file *outfile; 
struct arg_int  *crossings;
struct arg_int  *knotindex;
struct arg_int  *nverts;
struct arg_int  *nsamples;
struct arg_int  *maxtrials;
struct arg_int  *seed;

struct arg_lit  *verbose;
struct arg_lit  *help;
struct arg_lit  *quiet;

struct arg_end  *end;
struct arg_end  *helpend;

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void print_progress (double percentage,int found,int target)
{
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r%3d%% (%d/%d) [%.*s%*s]", val, found, target, lpad, PBSTR, rpad, "");
    fflush (stdout);
}


int main(int argc,char *argv[]) {

  int nerrors;

  void *argtable[] =
    {

      outfile = arg_filen("o","outfile","<filename>",0,1,"output VECT files will be named <filename>-1.vect, ..., <filename>-<s>.vect"),
      crossings = arg_int0("c","crossings","<crossing-number>","number of crossings in knot type"),
      knotindex = arg_int0("i","rolfsen-index","<index>","index in Rolfsen table for knot type"),
      nverts = arg_int1("n","number-of-vertices","<n>","number of vertices in polygon"),
      nsamples = arg_int1("s","number-of-samples","<s>","number of polygons of given knot typeto generate"),
      maxtrials = arg_int0("m","maximum-trials","<m>","maximum number of polygons to check for correct knot type"),
      seed = arg_int0(NULL,"seed","<s>","seed for random number generator"),
      verbose = arg_lit0(NULL,"verbose","print debugging information"),
      quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"),
      help = arg_lit0(NULL,"help","display help message"),
      end = arg_end(20)
    };

  void *helptable[] = {help,helpend = arg_end(20)};
  void *helpendtable[] = {helpend};

  const gsl_rng_type * T;
  gsl_rng *rng;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);

  int seedi;

  if (seed->count > 0) { seedi = seed->ival[0]; }
  else { seedi = time(0); }
  gsl_rng_set(rng,seedi);

  printf("randomknot compiled " __DATE__ " " __TIME__ "\n");
  printf("with %s random number gen, seeded with %d.\n",gsl_rng_name(rng),seedi);

  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("randomknot: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */

    nerrors = arg_parse(argc,argv,helptable);

    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      fprintf(stderr,"randomknot compiled " __DATE__ " " __TIME__ "\n");
      arg_print_errors(stdout,end,"randomknot");

      printf("usage\n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(1);

    } else {  /* The help table matched, which means we
		 asked for help or gave nothing */

      printf("\n"
	     "randomknot generates random equilateral polygons using the moment polytope\n"
	     "algorithm of Cantarella, Duplantier, Shonkwiler, and Uehara. If a knot type\n"
	     "is specified, the program calculates knot type (using the HOMFLY polynomial)\n"
	     "and saves only polygons of the given knot type.\n"
	     "\n"
	     "The number of vertices and samples must always be specified. If a knot type\n"
	     "is (crossing number and index) is specified, the user must also specify a\n"
	     "maximum number of trials.\n"
	     "\n"
	     "In cases where the same HOMFLY corresponds to more than one knot type,\n"
	     "all samples with the correct HOMFLY are saved, regardless of actual knot type.\n"
	     "\n"
	     "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }

  }

  if (crossings->count != knotindex->count) {

    fprintf(stderr,
	    "randomknot: Must specify both crossing number and index to specify knot type\n"
	    "            or must specify neither to obtain random polygons of any knot type\n");

    printf("usage\n\n");
    arg_print_glossary(stdout, argtable," %-25s %s\n");
    exit(1);

  }

  if (crossings->count > 0 && maxtrials->count == 0) {

    fprintf(stderr,
	    "randomknot: If knot type is specified, must also specify --maximum-trials\n"
	    "            to provide upper bound on number of polygons to check for correct\n"
	    "            knot type.\n");

    printf("usage\n\n");
    arg_print_glossary(stdout, argtable," %-25s %s\n");
    exit(1);

  }

  /******************** Now we deal with filenames *******************/

  char *ofname = NULL;

  if (outfile->count != 0) { 

    ofname = calloc(4096,sizeof(char));
    strncpy(ofname,outfile->filename[0],4096);

  } else if (crossings->count > 0) {

    ofname = calloc(256,sizeof(char));
    sprintf(ofname,"%d-%d-%d-vert-example",crossings->ival[0],knotindex->ival[0],nverts->ival[0]);

  } else {

    ofname = calloc(256,sizeof(char));
    sprintf(ofname,"%d-vert-example",nverts->ival[0]);

  }

  /* Now we actually do the work */

  printf("plcurve version %s\n",PACKAGE_VERSION);

  char svntag[1024];
  sprintf(svntag,"%s",SVNVERSION);
  if (!strstr("exported",svntag)) {  /* We were built from svn */
      printf("svn version %s\n",SVNVERSION);
  }

  if (crossings->count > 0) {
    
    printf("generating %d %d-vertex examples of knot %d.%d (%d max trials)\n",
	   nsamples->ival[0],nverts->ival[0],crossings->ival[0],knotindex->ival[0],
	   maxtrials->ival[0]);
    fflush(stdout);

  } else {

     printf("generating %d %d-vertex random polygons\n",
	    nsamples->ival[0],nverts->ival[0]);
     fflush(stdout);

  }

  int found=0,trial=0;

  for(;trial<maxtrials->ival[0] && found<nsamples->ival[0];trial++) {

    plCurve *sample;
    plc_knottype *knottype;
    int nposs;

    sample = plc_random_equilateral_closed_polygon(rng,nverts->ival[0]);

    if (crossings->count > 0) {
      
      knottype = plc_classify(rng,sample,&nposs);

      if (knottype != NULL) { /* The classify call can fail. */
	
	int poss;

	for(poss=0;poss<nposs;poss++) {

	  if (knottype->nf == 1 && knottype->cr[0] == crossings->ival[0] &&
	      knottype->ind[0] == knotindex->ival[0]) {
	    
	    found++;
	    FILE *outfile;
	    char name[512];
	    
	    sprintf(name,"%s-%05d.vect",ofname,found);
	    outfile = fopen(name,"w");
	    if (outfile == NULL) {
	      fprintf(stderr,"randomknot: Couldn't open %s for output.",name);
	      exit(1);
	    }
	    
	    plc_write(outfile,sample);
	    fclose(outfile);

	    if (quiet->count == 0) { 
	      
	      print_progress((double)(found)/(double)(nsamples->ival[0]),found,nsamples->ival[0]);
	      
	    }
	    
	  }
	}
	
	free(knottype);
	
      }
      
    } else { /* knot type not specified; there's no reason to check it */

        found++;
	FILE *outfile;
	char name[512];

	sprintf(name,"%s-%05d.vect",ofname,found);
	outfile = fopen(name,"w");
	if (outfile == NULL) {
	  fprintf(stderr,"randomknot: Couldn't open %s for output.",name);
	  exit(1);
	}

	plc_write(outfile,sample);
	fclose(outfile);

	if (quiet->count == 0) {
	  
	  print_progress((double)(found)/(double)(nsamples->ival[0]),found,nsamples->ival[0]);

	}
	
    }

    plc_free(sample);

  }

  printf("done\n");
  free(ofname);
  
  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
  arg_freetable(helpendtable,sizeof(helpendtable)/sizeof(helpendtable[0]));

  exit(0);
}
