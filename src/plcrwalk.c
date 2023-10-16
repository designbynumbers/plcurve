/*

  This file generates closed and confined random walks using tsmcmc
  or the moment polytope algorithm. It can generate chordlength data
  as well; this is primarily used to make sure that the algorithm
  is operating correctly. 

*/

#include"plCurve.h"
#include"plcTopology.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <string.h>
#include <dirent.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include"argtable2.h"
#include <assert.h>
#include "tsmcmc.h"
#include <string.h>
#include <sys/stat.h>

// Turn asserts ON.
#define DEBUG 1 
#define DIBUF_SIZE 10000

/* Global variables live here. */

struct arg_int  *samples;      // number of samples to generate
struct arg_int  *n;            // number of edges
struct arg_int  *skip;         // number of walks to skip between samples
struct arg_int  *seed;         // specify a seed for the random number generator (for debugging)
struct arg_int  *mpr;          // number of times to repeat each moment polytope step (optional)
struct arg_dbl  *beta;         // fraction of hit-and-run vs dihedral angle steps (optional)
struct arg_dbl  *delta;        // fraction of permutation steps vs hit-and-run and dihedral steps (optional)
struct arg_dbl  *radius;       // radius of confining sphere (if any)
struct arg_int  *burnin;       // number of steps to discard in burn-in period (optional)
struct arg_int  *interval;     // sampling interval (optional)
struct arg_file *outfile;      // name for outfile (optional)
struct arg_str  *format;       // Vect, Mathematica, ChordLength (optional)
struct arg_str  *knottype;     // knot type of samples to look for
struct arg_dbl  *ftc;          // Failure to close
struct arg_int  *chordskip;    // If the output format is ChordLength, which chord length to record

struct arg_lit  *verbose;
struct arg_lit  *help;
struct arg_lit  *quiet;

struct arg_end *end;
struct arg_end *helpend;

gsl_rng *rng; /* The global random number generator */

int PD_VERBOSE = 0;

int main(int argc,char *argv[]) {

  int nerrors;

  void *argtable[] = 
    {
      n = arg_int1("n","nedges","<int>","number of edges in walk"),
      samples = arg_int1("s","samples","<int>","number of samples to generate"),
      knottype = arg_str0("k","knot-type","<3_1>","knot type of samples to generate"),

      radius = arg_dbl0("r","radius","<float>","radius of confinement sphere (optional)"),
      ftc = arg_dbl0(NULL,"failure-to-close","<float>","distance between first and last vertex (optional, don't use for closed polygons"),
      skip = arg_int0(NULL,"skip","<int>","number of walks to skip between samples"),
      burnin = arg_int0("b","burnin","<int < nsamples>","number of samples to discard"),
      mpr = arg_int0(NULL,"moment-polytope-repeat","<10>","number of times to repeat moment polytope step (optional)"),
      seed = arg_int0(NULL,"seed","<int>","seed for random number generator (optional)"),
      beta = arg_dbl0(NULL,"beta","<float> in [0,1]","fraction of non-permutation steps which are moment polytope steps"),
      delta = arg_dbl0(NULL,"delta","<float> in [0,1]","fraction of all steps which are permutation steps"),
     
      outfile = arg_file0("o","outfile","<filename>","filename for output tar file"),
      format = arg_str0("f","format","<Mathematica|VECT|ChordLength>","output format (optional)"),
      chordskip = arg_int0(NULL,"chordskip","<int> in [0,n]","number of edges to skip when measuring chord length"),

      verbose = arg_lit0(NULL,"verbose","print debugging information"),
      quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"), 
      help = arg_lit0(NULL,"help","display help message"),
      end = arg_end(20)
    };
  
  void *helptable[] = {help,helpend = arg_end(20)};
  void *helpendtable[] = {helpend};

  plc_knottype *searchtype = NULL;

  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("randompolygon: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      fprintf(stderr,"randompolygon compiled " __DATE__ " " __TIME__ "\n");
      arg_print_errors(stdout,end,"randompolygon");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      fprintf(stderr,"randompolygon compiled " __DATE__ " " __TIME__ "\n");
      printf("randompolygon generates closed (or fixed failure-to-close) \n"
	     "and possibly confined equilateral random walks using\n"
	     "the toric-symplectic moment polytope algorithm or\n"
	     "the (improved) moment polytope rejection sampling algorithm.\n"
	     "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  /* If we got this far, it's time to try to parse the arguments */

  clock_t bigstart,bigend;
  double cpu_time_used;

  const gsl_rng_type * rng_T;
     
  gsl_rng_env_setup();  
  rng_T = gsl_rng_default;
  rng = gsl_rng_alloc(rng_T);
  
  int seedi;
  
  if (seed->count > 0) { seedi = seed->ival[0]; }
  else { seedi = time(0); }
  gsl_rng_set(rng,seedi);

  tsmcmc_triangulation_t T;
  double runbeta, rundelta;
  int runmpr = 10,runburnin = 1;
  int runn, runskip;
  int chordskipk;
  double runradius;
  double runftc = 0.0;
  FILE *summaryfile = NULL;

  runn = n->ival[0];

  if (ftc->count > 0) { 

    runftc=ftc->dval[0];
    if (runftc <= 0.0 || runftc >= (double)(runn)) { fprintf(stderr,"randompolygon: failure-to-close must be in (0,n) for any polygons to exist"); exit(1);}

  }
  
  if (burnin->count > 0) { 

    runburnin = burnin->ival[0];
    if (runburnin < 0) { fprintf(stderr,"randompolygon: number of samples to discard during burnin = %d must be > 0",runburnin); exit(1);}

  }

  if (skip->count > 0) { 

    runskip = skip->ival[0];
    if (runskip < 0) { fprintf(stderr,"randompolygon: number of walks to skip between samples = %d must be >= 0",runskip); exit(1);}

  } else {

    runskip = 2*runn;
    
  }

  if (radius->count > 0) { 

    if (ftc->count > 0) { 

      fprintf(stderr,"randompolygon: can't generate confined walks with nonzero failure-to-close (yet)\n");
      exit(1);

    }

    T = tsmcmc_fan_triangulation(n->ival[0]);
    runbeta = 0.5; rundelta = 0;
    runradius = radius->dval[0];

  } else { 

    T = tsmcmc_spiral_triangulation(n->ival[0]);
    runbeta = 0.5; rundelta = 0.9;

  }


  if (beta->count > 0) { 

    runbeta = beta->dval[0]; 
    if (runbeta < 0 || runbeta > 1.0) { fprintf(stderr,"randompolygon: beta = %g must be in [0,1].",runbeta); exit(1); }

  } 

  if (delta->count > 0) { 

    rundelta = delta->dval[0];
    if (rundelta < 0 || rundelta > 1.0) { fprintf(stderr,"randompolygon: delta = %g must be in [0,1].",rundelta); exit(1); }

  }

  if (mpr->count > 0) { 

    runmpr = mpr->ival[0];
    if (runmpr < 0) { fprintf(stderr,"randompolygon: number of times to repeat moment polytope steps = %d must be > 0",runmpr); exit(1);}

  }

  typedef enum {Mathematica,VECT,ChordLength} outputformat_t;
  
  outputformat_t runof = VECT;

  if (format->count > 0) { 

    if (!strcmp(format->sval[0],"Mathematica") || 
	!strcmp(format->sval[0],"mathematica") ||
	!strcmp(format->sval[0],"MATHEMATICA") ||
	!strcmp(format->sval[0],"CSV") ||
	!strcmp(format->sval[0],"csv")) { 

      runof = Mathematica;

    } else if (!strcmp(format->sval[0],"VECT") || !strcmp(format->sval[0],"vect")) {

      runof = VECT;

    } if (!strcmp(format->sval[0],"ChordLength") || !strcmp(format->sval[0],"chordlength")) {

      runof = ChordLength;

    } else { 

      fprintf(stderr,"randompolygon: Can't parse output format %s (should be Mathematica, CSV or VECT)\n",format->sval[0]); 
      exit(1);

    }

  }

  if (runof == ChordLength) {

    if (chordskip->count > 0) {

      chordskipk = chordskip->ival[0];
      if (chordskipk < 0 || chordskipk > runn) {

	fprintf(stderr,"randompolygon: Can't collect chordlength data on skip %d chords in a %d-gon.\n",
		chordskipk,runn);
	exit(1);
      }

    } else {

      chordskipk = floor(runn/2.0);

    }

  }

  char outfile_name[512];
  
  if (outfile->count > 0) { 

    strncpy(outfile_name,outfile->filename[0],512);

  } else {

    if (runof == Mathematica || runof == VECT) {

      if (radius->count > 0) { 
	
	sprintf(outfile_name,"rws-%03d-edge-%.2f-rad",runn,runradius);
	
      } else {

	sprintf(outfile_name,"rws-%03d-edge",runn);
	
      }

    } else if (runof == ChordLength) {

      sprintf(outfile_name,"rws-chordlength-%03d-edge-%03d-skip.csv",runn,chordskipk);

    }
  }

  if (runof == ChordLength) {
    
    summaryfile = fopen(outfile_name,"w");
    if (summaryfile == NULL) {
      fprintf(stderr,"randompolygon: Couldn't open summary file %s.\n",outfile_name);
      exit(1);
    }
    
  }

  if (ftc->count == 0 && radius->count == 0) { runskip = 1; }
  /* There's no need to skip samples with the direct sampler */

  if (knottype->count > 0) {
    searchtype = plc_read_knottype(knottype->sval[0]);
    if (searchtype == NULL) {
      fprintf(stderr,"randompolygon: Can't parse desired knot type %s.\n",knottype->sval[0]);
      exit(1);
    }
  }
  
  printf("randompolygon closed/confined equilateral random walk generator\n");   
  printf("with %s random number gen, seeded with %d.\n",gsl_rng_name(rng),seedi);
  
  if (ftc->count == 0) { 
    printf("generating %d %d-edge closed random walks",samples->ival[0],runn);
  } else {
    printf("generating %d %d-edge equilateral random walks with failure-to-close %g\n",samples->ival[0],runn,runftc);
  }
  if (radius->count > 0) { printf(" confined in sphere of radius %4f around vertex 0",runradius); }
  printf("\n");

  if (runof == VECT) {
    printf("as polygons stored as VECT files in %s",outfile_name);
  } else if (runof == Mathematica) {
    printf("as polygons stored as Mathematica (csv) files in %s",outfile_name);
  } else if (runof == ChordLength) {
    printf("as skip-%d chordlengths only in %s",chordskipk,outfile_name);
  }
  printf("\n\n");

  if (verbose->count > 0) { 

    if (radius->count > 0 || ftc->count > 0) {
      
      printf("Algorithm Parameters\n"
	     "-------------------------------------------------\n"
	     "# of walks to skip between recorded samples: %04d\n"
	     "# of walks to discard during burnin period : %04d\n"
	     "# of times to repeat moment polytope steps : %04d\n"
	     "fraction of permutation steps (delta)      : %04g\n"
	     "fraction of moment polytope steps \n"
	     "  among non-permutation steps (beta)       : %04g\n"
	     "failure to close                           : %04g\n"
	     "-------------------------------------------------\n",
	     runskip,runburnin,runmpr,rundelta,runbeta,runftc);

    } else {

      printf("Algorithm Parameters\n"
	     "--------------------------------------------------\n"
	     "Using CSS/CDSU direct sampling algorithm \n");
    }
  }

  fflush(stdout);
 
  /* We start by creating a subdirectory for the randompolygon samples */

  if (runof == Mathematica || runof == VECT) {
    
    printf("opening sample directory (%s)...",outfile_name);
    
    DIR *sampledir;
    sampledir = opendir(outfile_name);
    
    if (sampledir != NULL) { /* Directory exists; kill it */
    
      struct dirent *thisfile;
      for(thisfile=readdir(sampledir);thisfile!= NULL;thisfile=readdir(sampledir)) { 
	
	if (strcmp(thisfile->d_name,".") && strcmp(thisfile->d_name,"..")) { 
	  
	  char fullname[2048];
	  sprintf(fullname,"./%s/%s",outfile_name,thisfile->d_name);
	  
	  if(remove(fullname) != 0) { 
	    
	    fprintf(stderr,"randompolygon: Couldn't remove file %s to make room for samples",fullname);
	    exit(1);
	    
	  }

	}
    
      }
    
      closedir(sampledir);
      if (remove(outfile_name) != 0) {
	
	fprintf(stderr,"randompolygon: Couldn't remove existing directory ./%s to make room for samples",outfile_name);
	exit(1);

      }
      
    }
   
    if (mkdir(outfile_name,S_IRWXU) != 0) { 

      fprintf(stderr,"randompolygon: Couldn't create directory ./%s to store sample files\n",outfile_name);
      exit(1);

    }

    printf("done\n");

  }

  /* Now initialize the sampler. */

  double *diagonal_lengths = NULL,
    *edge_lengths = NULL,
    *dihedral_angles = NULL;

  if (radius->count > 0) { 

    tsmcmc_confined_equilateral_ngon(rng,T,runradius,
				     &edge_lengths,&diagonal_lengths,&dihedral_angles);

  } else if (ftc->count > 0) {

    tsmcmc_failure_to_close_ngon(rng,T,runftc,&edge_lengths,&diagonal_lengths,&dihedral_angles);

  } else {

    /* The direct sampler requires no initialization */

  }

  /* Start the main loop */

  printf("generating samples...");
  bigstart = clock();   

  int i,samps;

  for(i=0,samps=0; 
      samps < samples->ival[0] && i < 10000000;
      i++) {

    if (radius->count > 0) { 

      tsmcmc_confined_dihedral_diagonal_step(rng,T,runradius,edge_lengths,
					     diagonal_lengths,dihedral_angles,runbeta,runmpr);

    } else if (ftc->count > 0) { 

      tsmcmc_dihedral_diagonal_permute_step(rng,T,edge_lengths,
					    diagonal_lengths,
					    dihedral_angles,
					    runbeta,rundelta,runmpr);

    }

    if (i > runburnin && (i - runburnin) % (runskip+1) == 0) { 

      plCurve *L;

      if (radius->count > 0 || ftc->count > 0) {
	
	L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedral_angles);

      } else {

	L = plc_random_equilateral_closed_polygon(rng,runn);

      }

      int passes_knottype = true;
      int nposs = 0;
      plc_knottype *thistype;
      int check;

      if (searchtype != NULL) {

	thistype = plc_classify(rng,L,&nposs);

	if (thistype == NULL) { /* We couldn't classify the knot: fail */

	  passes_knottype = false;

	} else {

	  passes_knottype = false;

	  for(check=0;check < nposs;check++) {

	    if (thistype->cr[check] == searchtype->cr[check] && thistype->ind[check] == searchtype->ind[check]) {

	      passes_knottype = true;

	    }

	  }

	}

      }

      if (passes_knottype) {

	if (runof == Mathematica || runof == VECT) {
     
	  char filename[2048];
	
	  if (runof == VECT) { 
	    
	    sprintf(filename,"./%s/walk-%05d.vect",outfile_name,samps);
	  
	  } else if (runof == Mathematica) {

	    sprintf(filename,"./%s/walk-%05d.csv",outfile_name,samps);

	  }

	  FILE *sample_outfile = fopen(filename,"w");
	  if (sample_outfile == NULL) {
	    fprintf(stderr,"randompolygon: Couldn't open output filename %s.\n",filename);
	    exit(1);
	  }
	  
	  if (runof == VECT) { 
	  
	    plc_write(sample_outfile,L);
	  
	  } else { 

	    int vt;
	    for(vt=0;vt<L->cp[0].nv;vt++) { 
	  
	      fprintf(sample_outfile,"%12g, %12g, %12g \n",plc_M_clist(L->cp[0].vt[vt]));

	    }
	
	  }

	  fclose(sample_outfile);

	} else if (runof == ChordLength) {

	  double cl = plc_distance(L->cp[0].vt[0],L->cp[0].vt[chordskipk]);
	  fprintf(summaryfile,"%.17g \n",cl);

	}
	
	samps++;

      }

      plc_free(L);

    }

  }
      
  bigend = clock();
  cpu_time_used = ((double)(bigend - bigstart))/CLOCKS_PER_SEC;
  printf("done (%-5.3g sec/%-5.3g sec/polygon/%-5.3f polygons/sec)\n\n",cpu_time_used,cpu_time_used/samps,samps/cpu_time_used);

  if (summaryfile != NULL) {
    fclose(summaryfile);
    summaryfile = NULL;
  }
 
  /***************************************/

  printf("clearing memory...");

  tsmcmc_triangulation_free(T);  
  gsl_rng_free(rng);
  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
  arg_freetable(helpendtable,sizeof(helpendtable)/sizeof(helpendtable[0]));

  if (edge_lengths != NULL) { free(edge_lengths); edge_lengths = NULL; }
  if (diagonal_lengths != NULL) { free(diagonal_lengths); diagonal_lengths = NULL; }
  if (dihedral_angles != NULL) { free(dihedral_angles); dihedral_angles = NULL; }

  printf("done\n");

  /*************************************/

  exit(0);
}
