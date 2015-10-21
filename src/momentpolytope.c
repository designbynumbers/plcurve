/*

  This file generates closed and confined random walks using tsmcmc. 

*/

#include"plCurve.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <argtable2.h>
#include <assert.h>
#include "tsmcmc.h"
#include <string.h>
#include <sys/stat.h>

// Turn asserts ON.
#define DEBUG 1 
#define DIBUF_SIZE 10000

/* Global variables live here. */

struct arg_int  *eqn;          // equilateral polygon with n edges of length 1

struct arg_dbl  *ftc;          // Failure to close polygon with n-1 edges of length 1 and one edge of length ftc.
struct arg_int  *ftcn;         // Number of edges in failure-to-close polygon

struct arg_dbl  *edgelengths_arg;   // arbitrary edgelengths

struct arg_str  *triangulation; 

struct arg_int  *seed;         // specify a seed for the random number generator (for debugging)
struct arg_file *outfilename;  // name for outfile (optional)

struct arg_lit  *integral_arg;     // Are edgelengths all integers?
struct arg_lit  *mathematica;  // Also output triangulation in Mathematica form.

struct arg_lit  *verbose;
struct arg_lit  *help;
struct arg_lit  *quiet;
struct arg_end  *end;
struct arg_end  *helpend;

struct arg_rem *tab1;
struct arg_rem *tab2;
struct arg_rem *tab3;

gsl_rng *rng; /* The global random number generator */

int main(int argc,char *argv[]) {

  int nerrors;

  void *argtable[] = 
    {
      tab2 = arg_rem("Type of polygon","(one of these must be given)"),
      tab3 = arg_rem("-----","------"),
      
      eqn = arg_int0("n","equilateral-n","<int>","closed walk with <int> edges of length 1"),

      ftc = arg_dbl0(NULL,"ftc-dist","<float>","distance between first and last vertex (optional)"),
      ftcn = arg_int0(NULL,"ftc-n","<n>","number of edges in failure to close polytope"),

      edgelengths_arg = arg_dbln("e","edgelength","<float or int>",0,1000,"arbitrary edgelengths in order"),

      tab1 = arg_rem("-----","------"),

      triangulation = arg_str0("t","triangulation","<fan|spiral|teeth|random>","triangulation type (default is fan)"),
      integral_arg = arg_lit0("i","integral","interpret edgelengths as integers"),

      seed = arg_int0(NULL,"seed","<int>","seed for random number generator (optional)"),
      outfilename = arg_file0("o","outfile","<filename>","filename for output file"),

      //mathematica = arg_lit0(NULL,"Mathematica","also output triangulation in Mathematica form"),

      verbose = arg_lit0(NULL,"verbose","print debugging information"),
      quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"), 
      help = arg_lit0(NULL,"help","display help message"),
      end = arg_end(20)
    };
  
  void *helptable[] = {help,helpend = arg_end(20)};
  void *helpendtable[] = {helpend};

  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("momentpolytope: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      fprintf(stderr,"momentpolytope compiled " __DATE__ " " __TIME__ "\n");
      arg_print_errors(stdout,end,"momentpolytope");

      printf("usage\n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      fprintf(stderr,"momentpolytope compiled " __DATE__ " " __TIME__ "\n");
      printf("momentpolytope generates polymake-format polytopes for the \n"
	     "moment polytope of three kinds of polygons:  \n"
	     "\n"
	     "equilateral polygons with n length 1 edges\n"
	     "   (--equilateral-n option)\n"
	     "\n"
	     "'failure-to-close' polygons (n-1 length 1 edges) and one\n"
	     "   edge of arbitrary length (--ftc, and --ftcn options)\n"
	     "\n"
	     "polygons with arbitrary edgelengths\n"
	     "   (-e or --edgelength options)\n"
	     "\n"
	     "Optionally, you can also output a description of the\n"
	     "triangulation in Mathematica format.\n\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  if (eqn->count == 0 && ftc->count == 0 && ftcn->count == 0 && edgelengths_arg->count == 0) { 

    fprintf(stderr,"momentpolytope compiled " __DATE__ " " __TIME__ "\n");
    printf("momentpolytope generates polymake-format polytopes for the \n"
	   "moment polytope of three kinds of polygons:  \n"
	   "\n"
	   "equilateral polygons with n length 1 edges\n"
	   "   (--equilateral-n option)\n"
	   "\n"
	   "'failure-to-close' polygons (n-1 length 1 edges) and one\n"
	   "   edge of arbitrary length (--ftc, and --ftcn options)\n"
	   "\n"
	   "polygons with arbitrary edgelengths\n"
	   "   (-e or --edgelength options)\n"
	   "\n"
	   "Optionally, you can also output a description of the\n"
	   "triangulation in Mathematica format.\n\n"
	   "usage: \n\n");
    arg_print_glossary(stdout, argtable," %-25s %s\n");
    exit(0);

  }

  /* If we got this far, it's time to try to parse the arguments */

  const gsl_rng_type * rng_T;
     
  gsl_rng_env_setup();  
  rng_T = gsl_rng_default;
  rng = gsl_rng_alloc(rng_T);
  
  int seedi;
  
  if (seed->count > 0) { seedi = seed->ival[0]; }
  else { seedi = time(0); }
  gsl_rng_set(rng,seedi);

  tsmcmc_triangulation_t T;
  bool integral = false;
  FILE *outfile;
  double *edgelengths;
  int nedges;

  char edgelength_description[256];
  char triangulation_description[256];

  /* We need to EITHER have equilateral mode, FTC mode, or arbitrary edgelength mode. */

  if (eqn->count > 0 && (ftc->count > 0 || ftcn->count > 0 || edgelengths_arg->count > 0) ) { 

    printf("momentpolytope: Can't use (-n or --equilateral-n) options and --failure-to-close or -e. (Pick one or the other.)\n");
    exit(1);

  } 

  if ((ftc->count > 0 || ftcn->count > 0) && edgelengths_arg->count > 0) {

    printf("momentpolytope: Can\'t use both --failure-to-close and --edgelength options. (Pick one or the other.)\n");
    exit(1);

  }

  /* Now we need to create the vector of edgelengths */

  int i;

  if (eqn->count > 0) { 

    nedges = eqn->ival[0];
    assert(nedges > 3 && nedges < 1000000);
    edgelengths = calloc(nedges,sizeof(double));
    
    for(i=0;i<nedges;i++) { 
      
      edgelengths[i] = 1.0;

    }

    integral = true;
    sprintf(edgelength_description,"equilateral");


  } else if ((ftc->count > 0) || (ftcn->count > 0)) { 

    if (ftc->count == 0 || ftcn->count == 0) { 

      printf("momentpolytope: Can't specify --failure-to-close-dist without --failure-to-close-n.\n");
      exit(1);

    }

    if (ftc->dval[0] < 0) { 

      printf("momentpolytope: Can't have a negative failure to close distance.\n");
      exit(1);

    }

    if (ftcn->ival[0]-1 < ftc->dval[0]) { 

      printf("momentpolytope: Can't have a number of unit length edges (%d) less than failure to close distance of %g.\n",ftcn->ival[0]-1,ftc->dval[0]);
      exit(1);

    }

    nedges = ftcn->ival[0];
    assert(nedges > 3 && nedges < 1000000);

    edgelengths = calloc(nedges,sizeof(double));
    assert(edgelengths != NULL);

    edgelengths[0] = ftc->dval[0]; 
    for(i=1;i<nedges;i++) { 
     
      edgelengths[i] = 1.0;
    
    }

    if (fabs(edgelengths[0] - round(edgelengths[0])) < 1e-8 && integral_arg->count == 0) { 

      printf("momentpolytope: Warning! Failure to close distance is an integer (%g), but --integral option not specified.\n"
	     "                Will continue the computation in floating point mode.\n",edgelengths[0]);

    }

    sprintf(edgelength_description,"failure-to-close-%03.3g",edgelengths[0]);

  } else { 

    if (edgelengths_arg->count < 4) { 

      printf("momentpolytope: In arbitrary edgelength mode, at least 4 edgelengths must be supplied (read %d edgelengths).\n",edgelengths_arg->count);
      exit(1);

    }

    nedges = edgelengths_arg->count;
    assert(nedges < 1000000);

    edgelengths = calloc(nedges,sizeof(double));
    assert(edgelengths != NULL);

    for(i=0;i<nedges;i++) { 

      edgelengths[i] = edgelengths_arg->dval[i];

    }

    sprintf(edgelength_description,"arbitrary-edgelength");
    
  }

  /* We now check for integral edgelengths. */

  if (integral_arg->count > 0) { 

    for(i=0;i<nedges;i++) { 

      edgelengths[i] = round(edgelengths[i]);

    }

    integral = true;

  } 

  /* Now we parse the triangulation type and build the triangulation itself. */

  if (triangulation->count == 0 || strstr(triangulation->sval[0],"fan") != NULL) { 

    T = tsmcmc_fan_triangulation(nedges);
    sprintf(triangulation_description,"fan");

  } else if (strstr(triangulation->sval[0],"spiral") != NULL) { 

    T = tsmcmc_spiral_triangulation(nedges);
    sprintf(triangulation_description,"spiral");

  } else if (strstr(triangulation->sval[0],"teeth") != NULL) {

    T = tsmcmc_teeth_triangulation(nedges);
    sprintf(triangulation_description,"teeth");

  } else if (strstr(triangulation->sval[0],"random") != NULL) { 

    T = tsmcmc_random_triangulation(rng,nedges);
    sprintf(triangulation_description,"random");

  } else {

    printf("momentpolytope: The type of triangulation (%s) wasn't recognized.\n"
	   "                Valid options are fan, teeth, spiral, and random.\n",triangulation->sval[0]);
    exit(1);

  }

  /******************** Now we deal with filenames *******************/

  char filename[4096];

  if (outfilename->count == 0) { /* We'll have to make something up */
    
    sprintf(filename,"%s-%s-%d-edge-polytope.txt",edgelength_description,triangulation_description,nedges);
    
  } else {

    strncpy(filename,outfilename->filename[0],sizeof(filename));

  }

  /* Now we actually do the work */
     
  printf("momentpolytope polymake-format polytope generator for polygons\n");   
  printf("with %s random number gen, seeded with %d.\n",gsl_rng_name(rng),seedi);

  printf("generating moment polytope for %d edge polygon of type %s\n"
	 "with %s triangulation and edgelengths\n\n\t",nedges,edgelength_description,triangulation_description);

  for(i=0;i<nedges;i++) {

    printf("%03.3g ",edgelengths[i]);
    if ((i+1) % 5 == 0) { printf("\n\t"); }

  } 

  printf("\n\n");

  printf("opening output file %s...",filename);
  outfile = fopen(filename,"w");

  if (outfile == NULL) { 

    printf("fail\n"); 
    exit(1);

  }

  printf("done\n");

  printf("writing polytope...");
  tsmcmc_triangulation_polymake(outfile,T,edgelengths,integral);
  printf("done\n");

  printf("clearing memory...");

  fclose(outfile);
  free(edgelengths);
  tsmcmc_triangulation_free(T);

  gsl_rng_free(rng);
  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
  arg_freetable(helpendtable,sizeof(helpendtable)/sizeof(helpendtable[0]));

  printf("done\n");

  exit(0);
}
