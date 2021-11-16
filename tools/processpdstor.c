/*

  This program spits out the pdstor files in Ken Millett's ccode format
  in order to allow checking against Eric R's knot identification code.

*/

#include <plCurve.c>
#include <plcTopology.h>
#include <pd_multidx.h>
#include <pd_cyclic.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include <argtable2.h>
#include <assert.h>
#include <string.h>
#include <sys/stat.h>

//#include <Python.h>
//#include <pdcode/diagram_api.h>

// Turn asserts ON.
#define DEBUG 1

/* Global variables live here. */

struct arg_file *infile;  //
struct arg_file *outfile; // optional outfile override
struct arg_int  *numcomps;

struct arg_lit  *allcrossings;
struct arg_lit  *allorientations;
struct arg_lit  *orbitreps;
struct arg_lit  *knotsonly;
struct arg_lit  *verbose;
struct arg_lit  *help;
struct arg_lit  *quiet;
struct arg_lit  *KnotTheory;
struct arg_lit  *simplify;
struct arg_lit  *cantsimplify;
struct arg_lit  *knottype;
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
      allcrossings = arg_lit0("e","generate-all-crossing-signs","store a version of the pdcode with all crossing sign choices"),
      allorientations = arg_lit0(NULL,"generate-all-orientations","store a version of the pdcode with all orientation choices"),
      orbitreps = arg_lit0(NULL,"orbit-reps","when generating versions with all crossing sign/orientation choices, only store 1 representative for each orbit of the symmetry group"), 
      knotsonly = arg_lit0("k","knots-only","only process pdcodes for knots (1 component)\n"),
      numcomps = arg_int0("c","number-of-components","<n>","only process pdcodes for links with <n> components\n"), 
      verbose = arg_lit0(NULL,"verbose","print debugging information"),
      quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"),
      KnotTheory = arg_lit0("K","KnotTheory","print pd codes in the style of knottheory (WARNING: WILL NOT HANDLE SPLIT LINKS)"),
      simplify = arg_lit0("s","simplified-diagrams","output simplified diagrams (unless they have zero crossings)"),
      cantsimplify = arg_lit0(NULL,"no-simplifications-possible","don't output diagrams which can be simplified"),
      knottype = arg_lit0("","knottypes","output HOMFLY/knottype instead of ccode"),
      help = arg_lit0(NULL,"help","display help message"),
      end = arg_end(20)
    };

  void *helptable[] = {help,helpend = arg_end(20)};
  void *helpendtable[] = {helpend};

  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("processpdstor: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */

    nerrors = arg_parse(argc,argv,helptable);

    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      fprintf(stderr,"processpdstor compiled " __DATE__ " " __TIME__ "\n");
      arg_print_errors(stdout,end,"processpdstor");

      printf("usage\n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(1);

    } else {  /* The help table matched, which means we
		 asked for help or gave nothing */

      fprintf(stderr,"processpdstor compiled " __DATE__ " " __TIME__ "\n");
      printf("processpdstor converts a file of pd_codes in the pdstor text file format\n"
	     "to one of three outputs:\n"
	     "\n"
	     "  1) the Millett/Ewing ccode format (default) \n"
	     "  2) the KnotTheory PD code format (if the -K flag is set) \n"
	     "  3) the HOMFLY polynomial or knot type (if the --knottypes flag is set) \n"
	     "\n"
	     "If\n"
	     "\n"
	     "--generate-all-crossing-signs is set, generate results with all assignments of crossing signs\n"
	     "\n"
	     "--generate-all-orientations is set, generate results with all component orientations\n"
	     "\n"
	     "  (Both options may be set at the same time.)\n"
	     "\n"
	     "If --orbit-reps is set, then all 2^n assignments of crossings to an n-crossing diagram are\n"
	     "generated, but if the diagram has a symmetry group, only one representative of each symmetry\n"
	     "orbit is stored. For the 'standard trefoil diagram'\n"
	     "\n"
	     "       ---------                \n"
	     "       |       |                \n"
	     "   ----+-------+----            \n"
	     "   |   |       |   |            \n"
	     "   ----+--------   |            \n"
	     "       |           |            \n"
	     "       -------------            \n"
	     "\n"
	     "there would be 8 assignments of +/- to these three crossings, but \n"
	     "there is a 6-fold symmetry group which permutes the crossings \n"
	     "(rotate by 120 or flip) so there are only 4 orbit representatives: \n"
	     "+++, ++-, +--, and ---.\n"
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

  /* If the simplify flag is set we need to intialize
     the python libraries*/
  //  if(simplify->count != 0) {
  //  printf("Initializing Python libraries.\n");
  //  Py_Initialize();
  //  import_libpl__pdcode__diagram();
  //}

  printf("processpdstor %s crossing code generator for pdstor files\n",PACKAGE_VERSION);

  char svntag[1024];
  sprintf(svntag,"%s",SVNVERSION);
  if (!strstr("exported",svntag)) {  /* We were built from svn */
      printf("svn version %s\n",SVNVERSION);
  }
  
  printf("Built %s, %s.\n", __DATE__ , __TIME__ );

  if (knottype->count > 0) {
    printf("generating knot types for");
  } else {
    printf("generating crossing codes for");
  }

  printf(" %d pdstor files\n",infile->count);

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

      if (KnotTheory->count == 0 && knottype->count == 0) { 

	ofname = mangle(infile->basename[i],".pdstor",".ccode");

      } else if (KnotTheory->count > 0) {

	ofname = mangle(infile->basename[i],".pdstor",".KnotTheory");

      } else { /* We must be generating knot types */

	ofname = mangle(infile->basename[i],".pdstor",".knottypes");

      }
      
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

    if (knottype->count == 0) {
      
      printf("writing crossing codes...\n");

    } else {

      printf("writing knot types...\n");

    }
    
    int j;
    int codes_written = 0;
    for(j=0;!feof(in);j++) {

      if (verbose->count > 0) {

	PD_VERBOSE = 50;

      } else if (quiet->count > 0) {

	PD_VERBOSE = 0;

      } else {

	PD_VERBOSE = 10;

      }

      pd_code_t *inpd;
      inpd = pd_read(in);
      pd_uid_t inuid = inpd->uid;
      pd_idx_t incross = inpd->ncross;

      if (simplify->count > 0 || cantsimplify->count > 0) {

	pd_code_t *workpd = pd_simplify(inpd);
	assert(pd_ok(workpd));
	pd_code_free(&inpd);
	inpd=workpd;

      }

      if (
	  (cantsimplify->count==0 && inpd->ncross > 0) ||
	  (cantsimplify->count > 0 && inpd->ncross == incross) ||
	  (knottype->count > 0) ) {
      
	assert(inpd != NULL);

	if (!pd_ok(inpd)) {

	  pd_printf("pdcode read from file does not pass pd_ok\n",inpd);
	  exit(1);

	}

	if ((knotsonly->count == 1 && inpd->ncomps == 1) ||
	    (numcomps->count == 1 && numcomps->ival[0] == inpd->ncomps) ||
	    (numcomps->count == 0 && knotsonly->count == 0)) {

	  pd_or_t *component_orientations;
	  component_orientations = calloc(inpd->ncomps,sizeof(pd_or_t));
	  int i;
	  for(i=0;i<inpd->ncomps;i++) {
	    component_orientations[i] = PD_POS_ORIENTATION;
	  }
	  if (allcrossings->count > 0) {
	    for(i=0;i<inpd->ncross;i++) {
	      inpd->cross[i].sign = PD_NEG_ORIENTATION;
	    }
	  } // Otherwise, we'll want to keep the original crossing signs.

	  pd_stor_t *expansions = NULL;

	  if (orbitreps->count > 0) {
	    
	    expansions = pd_new_pdstor();

	  }

	  pd_multidx_t *orientation_idx;
	  pd_idx_t*    *comp_orientations;
	  pd_idx_t      one = 1, two = 2;
	  comp_orientations = calloc(inpd->ncomps,sizeof(pd_idx_t *));

	  if (allorientations->count > 0) {

	    int i;
	    for(i=0;i<inpd->ncomps;i++) { comp_orientations[i] = &two; }

	  } else {

	    int i;
	    for(i=0;i<inpd->ncomps;i++) { comp_orientations[i] = &one; }

	  }

	  orientation_idx = pd_new_multidx(inpd->ncomps,
					   (void **)(comp_orientations),
					   cyclic_ops);

	  unsigned int norientations = pd_multidx_nvals(orientation_idx);
	  unsigned int orcount;

	  for(orcount=0;
	      orcount < norientations;
	      orcount++,pd_increment_multidx(orientation_idx)) {

	    pd_multidx_t *crossing_idx;
	    pd_idx_t*    *comp_crossings;

	    comp_crossings = calloc(inpd->ncross,sizeof(pd_idx_t *));

	    if (allcrossings->count > 0) {

	      int i;
	      for(i=0;i<inpd->ncross;i++) { comp_crossings[i] = &two; }

	    } else {

	      int i;
	      for(i=0;i<inpd->ncross;i++) { comp_crossings[i] = &one; }

	    }

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
		  component_orientations[i] = PD_NEG_ORIENTATION;
		} else {
		  component_orientations[i] = PD_POS_ORIENTATION;
		}
	      }

	      if (allcrossings->count > 0) {

		for(i=0;i<working_pd->ncross;i++) {
		  if (crossing_idx->i[i] == 1) {
		    working_pd->cross[i].sign = PD_NEG_ORIENTATION;
		  } else {
		    working_pd->cross[i].sign = PD_POS_ORIENTATION;
		  }
		}

	      } else { /* We're only doing one sign assignment for the crossings; if we've
			  GOT crossing signs, keep them; otherwise, we may as well set them 
		          all to "positive" as we're outputting diagrams (always). */

		for(i=0;i<working_pd->ncross;i++) {
		  if (working_pd->cross[i].sign == PD_UNSET_ORIENTATION) {
		    working_pd->cross[i].sign = PD_POS_ORIENTATION;
		  }
		}

	      }

	      if (orbitreps->count > 0) { /* We are not outputting the results yet */

		pd_addto_pdstor(expansions,working_pd,DIAGRAM_ISOTOPY); /* Store, deleting duplicates. */

	      } else {

		/* We want to write out the results while we still have the component orientations 
		   and crossing signs recorded. */

		codes_written++;
	      
		fprintf(out,"# %d crossing pd from original UID %lu, orientations ",working_pd->ncross,inuid);
		for(i=0;i<working_pd->ncomps;i++) {
		  fprintf(out,"%c",pd_print_or(component_orientations[i]));
		}

		fprintf(out, " crossings ");
		for(i=0;i<working_pd->ncross;i++) {
		  fprintf(out,"%c",pd_print_or(working_pd->cross[i].sign));
		}
	      
		fprintf(out,"\n"); 
	      
		/*If the KnotTheory flag is set we use
		  pd_write_KnotTheory to write codes*/
		if(KnotTheory->count != 0) {
		  
		  pd_write_KnotTheory(out,working_pd);

		} else if (knottype->count != 0) {

		  if (working_pd->ncomps == 1) {

		    plc_knottype *knottype;
		    int nposs,i;
		    
		    knottype = pd_classify(working_pd,&nposs);
		    for(i=0;i<nposs;i++) {
		      plc_write_knottype(out,knottype[i]);
		    }
		    free(knottype);

		  } else {

		    char *homfly;
		    homfly = pd_homfly_timeout(working_pd,120);
		    fprintf(out,"%s",homfly);
		    free(homfly);

		  }
		    
		} else {
	      
		/*If KnotTheory and knottype not flagged not set then
		  we write codes as ccodes*/
	      
		  char *ccode = pdcode_to_ccode(working_pd);
		  fprintf(out,"%s\n",ccode);
		  free(ccode);
		}

	      }
	    
	      pd_code_free(&working_pd);

	    }

	    pd_free_multidx(&crossing_idx);
	    free(comp_crossings);

	  }

	  pd_free_multidx(&orientation_idx);
	  free(comp_orientations);
	  free(component_orientations);

	  /* Now, if we stored orbit-reps, we're going to have to write out the elements in the pdstor
	     (since we haven't written anything yet). */
	
	  if (orbitreps->count > 0) { 
	
	    pd_code_t *thispd;
	    
	    for (thispd = pd_stor_firstelt(expansions);thispd != NULL;thispd = pd_stor_nextelt(expansions)) {
	      
	      codes_written++;
	      
	      fprintf(out,"# %d crossing pd from original UID %lu, crossings ",
		      thispd->ncross,inuid);
	      for(i=0;i<thispd->ncross;i++) {
		fprintf(out,"%c",pd_print_or(thispd->cross[i].sign));
	      }
	      
	      /* fprintf(out," orientation "); */
	      /* for(i=0;i<thispd->ncomps;i++) { */
	      /*   fprintf(out,"%c",pd_print_or(component_orientations[i])); */
	      /* } */
	      fprintf(out,"\n"); 
	      
	      /*If the KnotTheory flag is set we use
		pd_write_KnotTheory to write codes*/
	      if(KnotTheory->count != 0) {

		pd_write_KnotTheory(out,thispd);

	      } else if (knottype->count > 0) {
		
	        if (thispd->ncomps == 1) {

		    plc_knottype *knottype;
		    int nposs,i;
		    
		    knottype = pd_classify(thispd,&nposs);
		    for(i=0;i<nposs;i++) {
		      plc_write_knottype(out,knottype[i]);
		    }
		    free(knottype);

		  } else {

		    char *homfly;
		    homfly = pd_homfly_timeout(thispd,120);
		    fprintf(out,"%s",homfly);
		    free(homfly);

		  }
		
	      } else {

		/*If KnotTheory flag and knottype flags are both not set then
		  we write codes as ccodes*/
	     
		char *ccode = pdcode_to_ccode(thispd);
		fprintf(out,"%s\n",ccode);
		free(ccode);
		
	      }
	      
	      pd_code_free(&thispd);
	      
	    }
	    
	    pd_free_pdstor(&expansions);

	  }
	  
	  pd_code_free(&inpd);
	
	}

      }
      
    }
    
    if (knottype->count > 0) {
      
      printf("done (wrote %d knot types or homfly polynomials).\n",codes_written);

    } if (KnotTheory->count > 0) {

      printf("done (wrote %d KnotTheory pd codes).\n",codes_written);

    } else {

      printf("done (wrote %d crossing codes).\n",codes_written);

    }

    if (outfile->count == 0) {
      printf("closing output file %s...",ofname);
      fclose(out);
    }

  }

  if(simplify->count != 0) {
      //Py_Finalize();
  }

  printf("done\n");

  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
  arg_freetable(helpendtable,sizeof(helpendtable)/sizeof(helpendtable[0]));

  printf("done\n");

  exit(0);
}
