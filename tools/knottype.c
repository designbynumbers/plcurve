/*******************************************************************

  knottype.c : Uses plc_homfly code to identify knots.

	   ************************************************************/

#include"plCurve.h"
#include"plcTopology.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/*#include<../utilib/mangle.h>
  %#include<../utilib/ordie.h>*/
#include<argtable2.h>

#include<polynomials.h>

/* Global variables live here. */

struct arg_lit *print_homfly; // output homfly polynomial
struct arg_int *timeout; // timeout for homfly computation
struct arg_lit *latex;  // output in LaTeX format
struct arg_int *seed;         // specify a seed for the random number generator (for debugging)
struct arg_lit *mathematica; // output in Mathematica format
struct arg_lit *print_pdcode;  // display pdcode
struct arg_lit *classify; // attempt to determine knot type
struct arg_lit *verify; // check whether knot type matches filename
struct arg_lit *composite; // treat multiple input files as summands in composite knot

struct arg_lit *verbose;
struct arg_file *corefile;
struct arg_lit  *help;
struct arg_lit  *quiet;

struct arg_end *end;
struct arg_end *helpend;

plCurve *core;
FILE    *infile_fptr,*outfile_fptr;

int    QUIET=0;
int    tmout = 30;

gsl_rng *rng; /* The global random number generator */

/****************************** Main procedure ********************************/

int main(int argc,char *argv[])
{
  int         infilenum,nerrors;

  void *argtable[] = 
    {
      print_pdcode = arg_lit0("p","print-pdcode","print pd-code for vect file"),
      timeout = arg_int0("t","timeout","<int>","maximum number of seconds to spend computing HOMFLY (default 30)"),
      latex = arg_lit0("l","latex","print output in LaTeX format"), 
      mathematica = arg_lit0("m","mathematica","print output in Mathematica format"),
      print_homfly = arg_lit0("h","homfly","print HOMFLYPT polynomial"),
      classify = arg_lit0("k","classify","classify the knot type"),
      composite = arg_lit0("#","composite","treat multiple inputs as summands in composite knot"),
      verify = arg_lit0("v","verify","verify knot type matches filename"),
      seed = arg_int0(NULL, "seed", "<int>", "seed for random number generator (optional)"),
      verbose = arg_lit0(NULL,"verbose","print debugging information"),
      corefile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
      quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"), 
      help = arg_lit0(NULL,"help","display help message"),
      end = arg_end(20)
    };
  
  void *helptable[] = {help,helpend = arg_end(20)};

  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("knottype: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
			   first set of errors was probably more
			   helpful, so we display it. */

      fprintf(stderr,"knottype (" PACKAGE_STRING ") compiled " __DATE__ " " __TIME__ "\n");
      arg_print_errors(stdout,end,"knottype");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      fprintf(stderr,"knottype (" PACKAGE_STRING ") compiled " __DATE__ " " __TIME__ "\n");
      printf("knottype classifies VECT file, mostly using HOMFLYPT polynomial.\n"
	     "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  const gsl_rng_type *rng_T;
  
  gsl_rng_env_setup();
  rng_T = gsl_rng_default;
  rng = gsl_rng_alloc(rng_T);

  int seedi;

  if (seed->count > 0) { seedi = seed->ival[0]; }
  else { seedi = (int) (time(0)); }
  gsl_rng_set(rng, seedi);

  QUIET = quiet->count > 0; /* Register if we're in batch mode */
 
  if (!QUIET) {
    printf("knottype (" PACKAGE_STRING ") compiled " __DATE__ " " __TIME__ "\n");
    printf("with %s random number gen, seeded with %d.\n", gsl_rng_name(rng), seedi);
  }

  fflush(stdout);
 
  /* Now we have parsed the arguments and are ready to work. */

  monomial *composite_poly = NULL;
  int ncomposite_monos = 0;

  if (timeout->count > 0) {

    tmout = timeout->ival[0];

  }
  
  for(infilenum = 0;infilenum < corefile->count;infilenum++) {
    
    /* We begin by loading the current core link. */
    
    infile_fptr = fopen(corefile->filename[infilenum],"r");
    
    if (infile_fptr == NULL) {
      
      fprintf(stderr,"knottype: Couldn't open file %s.\n",
	      corefile->filename[infilenum]);
      continue;  /* Try the next file */
      
    }
    
    int plr_error_num;
    char plr_error_str[1024];
    
    core = plc_read(infile_fptr,
		    &plr_error_num,plr_error_str,sizeof(plr_error_str));
    
    /* We now demonstrate the plCurve library's error handling protocol: */
    
    if (plr_error_num > 0) {   /* This is the signal for an error. */
      
      fprintf(stderr,"knottype: link reading error\n%s\n",plr_error_str);
      continue;  /* Try the next file */
      
    }
    
    fclose(infile_fptr);
     
    /* Report to the user */
    
    if (!QUIET) {
      
      printf("%s\n",corefile->basename[infilenum]);
      fflush(stdout);
    }
    
    /* Now we actually do the work. */
    
    if (print_pdcode->count > 0) {

      pd_code_t *core_pd = pd_code_from_plCurve(rng, core);
      pd_code_t *simplified_pd = pd_simplify(core_pd);
      pd_write_KnotTheory(stdout,simplified_pd);
      pd_code_free(&core_pd);
      pd_code_free(&simplified_pd);
    }
    
    if (print_homfly->count > 0) {
      
      char *chomfly;

      pd_code_t *core_pd = pd_code_from_plCurve(rng,core);
      chomfly = pd_homfly_timeout(core_pd,tmout);
      pd_code_free(&core_pd);
      
      printf("Homfly polynomial:(%s)\n",chomfly);
      
      if (verbose->count > 0) {
	
	char *compare;
	compare = lmpoly_check(chomfly);
	printf("Homfly check     :(%s)\n",compare);
	
	if (strcmp(chomfly,compare)) {
	  
	  fprintf(stderr,"Homfly and check homfly DIFFERENT. Aborting run.\n");
	  exit(1);
	  
	}

	free(compare);

      }

      free(chomfly);

    }

    if (latex->count > 0) {

      char *chomfly, *converted;

      pd_code_t *core_pd = pd_code_from_plCurve(rng,core);
      chomfly = pd_homfly_timeout(core_pd,tmout);
      pd_code_free(&core_pd);
      
      converted = lmpoly_to_latex(chomfly);
      printf("Homfly polynomial (latex): %s\n",converted);
      free(converted);
      free(chomfly);

    }

    if (mathematica->count > 0) {

      char *chomfly, *converted;

      pd_code_t *core_pd = pd_code_from_plCurve(rng,core);
      chomfly = pd_homfly_timeout(core_pd,tmout);
      pd_code_free(&core_pd);
      
      converted = lmpoly_to_mathematica(chomfly);
      printf("Homfly polynomial (mathematica): %s\n",converted);
      free(converted);
      free(chomfly);

    }     
     
    if (composite->count > 0) {

      if (infilenum == 0) {  /* At first, we just load the polynomial into composite_poly. */

	char *chomfly;

	pd_code_t *core_pd = pd_code_from_plCurve(rng,core);
	chomfly = pd_homfly_timeout(core_pd,tmout);
	pd_code_free(&core_pd);
	
	composite_poly = lmpoly_to_polynomial(chomfly,&ncomposite_monos);
	free(chomfly);

      } else { /* Afterwards, we take the product of the current homfly with the current composite_poly. */

	monomial *product,*poly;
	char *chomfly;
	int poly_monos,product_monos;

	pd_code_t *core_pd = pd_code_from_plCurve(rng,core);
	chomfly = pd_homfly_timeout(core_pd,tmout);
	pd_code_free(&core_pd);chomfly = pd_homfly_timeout(core_pd,tmout);

	poly = lmpoly_to_polynomial(chomfly,&poly_monos);
	product = product_polynomial(poly,poly_monos,composite_poly,ncomposite_monos,&product_monos);
	free(poly);
	free(composite_poly);

	composite_poly = product;
	ncomposite_monos = product_monos;

      }

    }
	 
    if (classify->count > 0) {

      plc_knottype *kt;
      int nposs;

      kt = plc_classify(rng,core,&nposs);
      
      if (kt == NULL) {
	
	 printf("Unable to classify\n");
	 
      } else {
	
	int i,j;
	 
	printf("Possible knot type(s):\n");
	 
	for(i=0;i<nposs;i++) {
	   
	  printf("|");
	   
	  for(j=0;j<kt[i].nf;j++) {
	     
	    printf("<%d_%d %s>",kt[i].cr[j],kt[i].ind[j],kt[i].sym[j]);
	     
	  }
	   
	  printf("|\n");
	   
	}

	free(kt);

      }

    }

    if (verify->count > 0) {

      plc_knottype *kt;
      int nposs;

      kt = plc_classify(rng,core,&nposs);
       
      if (kt == NULL) {
	
	printf("Unable to compute knot type from file.\n");
	printf("RESULT: Maybe ok\n");
	
      } else {
	 
	int i,j;
	
	printf("Possible knot type(s):\n");
	
	for(i=0;i<nposs;i++) {
	  
	  printf("|");
	  
	  for(j=0;j<kt[i].nf;j++) {
	    
	    printf("<%d_%d %s>",kt[i].cr[j],kt[i].ind[j],kt[i].sym[j]);
	    
	  }
	  
	  printf("|\n");
	  
	}
	
	printf("Knot type from name:\n");
	int cr,ind,comps;
	
	if (sscanf(corefile->basename[infilenum],"kl_%d_%d_%d",&cr,&comps,&ind) == 3) {

	  printf("%d^%d_%d\n",cr,comps,ind);
	  printf("Warning: Can't verify knot type of links in this version.\n");
	  
	  printf("RESULT: Maybe ok\n");
	  
	} else if (sscanf(corefile->basename[infilenum],"kl_%d_%d_",&cr,&ind) == 2) {
	  
	  int i;
	  bool ok = true;
	  bool maybe = false;
	  
	  printf("%d_%d\n",cr,ind);
	  
	  for(i=0;i<nposs;i++) {
	    
	    if (kt[i].cr[0] == cr && kt[i].ind[0] == ind) { maybe = true; }
	    if (kt[i].cr[0] != cr || kt[i].ind[0] != ind) { ok = false; }
	    
	  }
	  
	  if (ok) {
	    
	    printf("RESULT: Ok\n");
	    
	  } else if (maybe) {
	    
	    printf("RESULT: Maybe ok\n");
	    
	  } else {
	    
	    printf("RESULT: Not ok\n");
	    
	  }
	  
	}
	
      }
      
      free(kt);
      
    }
    
    /* Now we cleanup memory. */

    plc_free(core);
  
  }
  
  /* We now output the composite polynomial, if it exists. */
  
  if (print_homfly->count > 0 && composite->count > 0) {
    
    char *chomfly;
    chomfly = polynomial_to_lmpoly(composite_poly,ncomposite_monos);
    printf("\nComposite Homfly polynomial:(%s)\n",chomfly);
    
    if (verbose->count > 0) {
      
      char *compare;
      compare = lmpoly_check(chomfly);
      printf("Composite Homfly check     :(%s)\n",compare);
      
      if (strcmp(chomfly,compare)) {
	
	fprintf(stderr,"Homfly and check homfly DIFFERENT. Aborting run.\n");
	exit(1);
	
      }
      
      free(compare);
      
    }
    
    free(chomfly);
    
  }
  
  if (latex->count > 0) {
    
    char *converted;
    converted = poly_to_latex(composite_poly,&ncomposite_monos);
    printf("Composite Homfly polynomial (latex): %s\n",converted);
    free(converted);
    
  }
  
  if (mathematica->count > 0) {
    
    char *converted;
    converted = poly_to_mathematica(composite_poly,&ncomposite_monos);
    printf("Composite Homfly polynomial (mathematica): %s\n",converted);
    free(converted);
    
  }
  
  if (!QUIET) { printf("clearing memory..."); fflush(stdout); }

  gsl_rng_free(rng);
  arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

  if (!QUIET) { printf("done\n"); }
  
  return 0;

}




