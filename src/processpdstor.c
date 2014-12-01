/*

  This program spits out the pdstor files in Ken Millett's ccode format
  in order to allow checking against Eric R's knot identification code.

*/

#include<plCurve.c>
#include<plcTopology.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <string.h>

#include <argtable2.h>
#include <assert.h>
#include <string.h>
#include <sys/stat.h>

// Turn asserts ON.
#define DEBUG 1 

/* Global variables live here. */

struct arg_file *infile;  // 
struct arg_file *outfile; // optional outfile override

struct arg_lit  *allcrossings;
struct arg_lit  *allorientations;
struct arg_lit  *knotsonly;
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
      allcrossings = arg_lit0("e","generate-all-crossing-signs","store a version of the pdcode with all crossing sign choices"),
      allorientations = arg_lit0(NULL,"generate-all-orientations","store a version of the pdcode will all orientation choices"),
      knotsonly = arg_lit0("k","knots-only","only process pdcodes for knots (1 component)\n"),
      verbose = arg_lit0(NULL,"verbose","print debugging information"),
      quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"), 
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
	     "to the Millett/Ewing ccode format. For each pdcode, we either store the \n"
	     "ccode format -OR IF-\n"
	     "\n"
	     "--generate-all-crossing-signs is set, generate ccodes with all assignments of crossing signs\n"
	     "\n"
	     "--generate-all-orientations is set, generate ccodes with all component orientations\n"
	     "\n"
	     "Both options may be set at the same time.\n"
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
     
  printf("processpdstor crossing code generator for pdstor files\n");   
  printf("generating crossing codes for %d pdstor files\n",infile->count);

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
      ofname = mangle(infile->basename[i],".pdstor",".ccode");
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

    printf("writing crossing codes...\n");
    int j;
    int codes_written = 0;
    for(j=0;!feof(in);j++) { 
      
      PD_VERBOSE = 50;

      pd_code_t *inpd;
      inpd = pd_read(in);

      assert(inpd != NULL);

      if (!pd_ok(inpd)) {

	  pd_printf("pdcode read from file does not pass pd_ok\n",inpd);
	  exit(1);

      }

      if (knotsonly->count == 0 || inpd->ncomps == 1) { 

	pd_or_t *component_orientations;
	component_orientations = calloc(inpd->ncomps,sizeof(pd_or_t));
	int i;
	for(i=0;i<inpd->ncomps;i++) { 
	  component_orientations[i] = PD_POS_ORIENTATION;
	}
	for(i=0;i<inpd->ncross;i++) { 
	  inpd->cross[i].sign = PD_NEG_ORIENTATION;
	}
	
	/** INSERT COMPONENT SWITCHING STUFF HERE **/
	/** INSERT CROSSING SWITCHING STUFF HERE **/
	
	fprintf(out,"# %d crossing pd with UID %lu, crossings ",inpd->ncross,inpd->uid);
	for(i=0;i<inpd->ncross;i++) { 
	  fprintf(out,"%c",pd_print_or(inpd->cross[i].sign));
	}
	
	fprintf(out," orientation ");
	for(i=0;i<inpd->ncomps;i++) { 
	  fprintf(out,"%c",pd_print_or(component_orientations[i]));
	}
	fprintf(out,"\n");
	
	char *ccode = pdcode_to_ccode(inpd);
	fprintf(out,"%s\n",ccode);
      
	free(ccode);
	free(component_orientations);

	codes_written++;

      }

      pd_code_free(&inpd);
	
    }
    printf("done (wrote %d crossing codes).\n",codes_written);

    if (outfile->count == 0) { 
      printf("closing output file %s...",ofname);
      fclose(out);
    }

  }
    
  printf("done\n");

  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
  arg_freetable(helpendtable,sizeof(helpendtable)/sizeof(helpendtable[0]));

  printf("done\n");

  exit(0);
}
