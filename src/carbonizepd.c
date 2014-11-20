/*

  This program just acts as an interface to convert pdcodes from files to compilable
  c code. The purpose is mainly for testing and to allow programs to compile in a 
  version of the crossing database in order to simplify distributions. 

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

struct arg_lit  *verbose;
struct arg_lit  *weakchecking;
struct arg_lit  *help;
struct arg_lit  *quiet;
struct arg_end  *end;
struct arg_end  *helpend;

pd_code_t *pd_read_weak(FILE *infile)
{
  return NULL; /* This is a stub */
}

void pd_write_c_weak(FILE *outfile, pd_code_t *pd, char *pdname)
/* Writes a c procedure which recreates the pd code pd.
   The procedure will be called pd_create_(pdname). Unlike
   the library version, this custom version DOES NOT check
   pd_ok at the end of the run. */
{
  fprintf(outfile,"pd_code_t *pd_create_%s() { \n\n",pdname);
  fprintf(outfile,
	  "/* This procedure is machine generated by pd_write_c */\n"
	  "/* and probably shouldn't be hand-edited. */\n\n");

  fprintf(outfile,
	  "pd_code_t *pd;\n"
	  "pd = pd_code_new(%d);\n"
	  "assert(pd != NULL);\n",pd->MAXVERTS);

  fprintf(outfile,
	  "pd->ncross = %d;\n"
	  "pd->nedges = %d;\n"
	  "pd->ncomps = %d;\n"
	  "pd->nfaces = %d;\n"
	  "sprintf(pd->hash,\"%%s\",\"%s\");\n",
	  pd->ncross,pd->nedges,pd->ncomps,pd->nfaces,pd->hash);

  pd_idx_t i,j;

  fprintf(outfile,"\n/* Crossing data. */\n\n");

  /* Now rebuild the crossing buffer. */

  for(i=0;i<pd->ncross;i++) {

    for(j=0;j<4;j++) {

      fprintf(outfile,"pd->cross[%d].edge[%d] = %d;\n",i,j,pd->cross[i].edge[j]);

    }

    fprintf(outfile,"pd->cross[%d].sign = %d;\n\n",i,pd->cross[i].sign);

  }

  /* The edge buffer... */

  fprintf(outfile,"\n/* Edge data */\n\n");

  for(i=0;i<pd->nedges;i++) {

    fprintf(outfile,
	    "pd->edge[%d].head = %d;\n"
	    "pd->edge[%d].headpos = %d;\n"
	    "pd->edge[%d].tail = %d;\n"
	    "pd->edge[%d].tailpos = %d;\n\n",
	    i, pd->edge[i].head, i, pd->edge[i].headpos,
	    i, pd->edge[i].tail, i, pd->edge[i].tailpos);

  }

  /* Component data */

  fprintf(outfile,"\n/* Component Data */\n\n");

  for(i=0;i<pd->ncomps;i++) {

    fprintf(outfile,
	    "pd->comp[%d].nedges = %d;\n"
	    "pd->comp[%d].tag = '%c';\n\n",
	    i,pd->comp[i].nedges,i,pd->comp[i].tag);

    fprintf(outfile,
	    "pd->comp[%d].edge = calloc(pd->comp[%d].nedges,sizeof(pd_idx_t));\n"
	    "assert(pd->comp[%d].edge != NULL);\n\n",
	    i,i,i);

    for(j=0;j<pd->comp[i].nedges;j++) {

      fprintf(outfile,
	      "pd->comp[%d].edge[%d] = %d;\n",i,j,pd->comp[i].edge[j]);

    }

    fprintf(outfile,"\n");

  }

  /* Face data. */

  fprintf(outfile,"\n/* Face data */\n\n");

  for(i=0;i<pd->nfaces;i++) {

    fprintf(outfile,
	    "pd->face[%d].nedges = %d;\n"
	    "pd->face[%d].edge = calloc(pd->face[%d].nedges,sizeof(pd_idx_t));\n"
	    "pd->face[%d].or = calloc(pd->face[%d].nedges,sizeof(pd_or_t));\n"
	    "assert(pd->face[%d].edge != NULL);\n"
	    "assert(pd->face[%d].or != NULL);\n\n",
	    i,pd->face[i].nedges,i,i,i,i,i,i);

    for(j=0;j<pd->face[i].nedges;j++) {

      fprintf(outfile,
	      "pd->face[%d].edge[%d] = %d;\n"
	      "pd->face[%d].or[%d] = %d;\n\n",
	      i,j,pd->face[i].edge[j],
	      i,j,pd->face[i].or[j]);

    }

  }

  fprintf(outfile,
	  "\n/* End of data. */\n\n"
	  "assert(pd_ok(pd));\n"
	  "return pd;\n\n"
	  "}\n\n");

}

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


int main(int argc,char *argv[]) {

  int nerrors;

  void *argtable[] = 
    {
      
      infile = arg_filen(NULL,NULL,"<filename>",1,1000,"text format pd_code file"),
      outfile = arg_filen("o","outfile","<filename>",0,1,"filename for (single) output file"),
      weakchecking = arg_lit0(NULL,"weak-checking","do not test for pd_ok (don't do this unless you understand why your pd code is not passing pd_ok)"),
      verbose = arg_lit0(NULL,"verbose","print debugging information"),
      quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"), 
      help = arg_lit0(NULL,"help","display help message"),
      end = arg_end(20)
    };
  
  void *helptable[] = {help,helpend = arg_end(20)};
  void *helpendtable[] = {helpend};

  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("carbonizepd: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      fprintf(stderr,"carbonizepd compiled " __DATE__ " " __TIME__ "\n");
      arg_print_errors(stdout,end,"carbonizepd");

      printf("usage\n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      fprintf(stderr,"carbonizepd compiled " __DATE__ " " __TIME__ "\n");
      printf("carbonizepd converts pd codes in the pd_code text file format\n"
	     "to compileable C code which regenerates the same pd_code in memory\n"
	     "The purpose of this is to automatically generate pd codes which can\n"
	     "be baked into executables for testing and convenience purposes.\n"
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
     
  printf("carbonizepd C code generator for pdstor files\n");   
  printf("generating C code for %d pdstor files\n",infile->count);

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
      ofname = mangle(infile->basename[i],".pdstor",".c");
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

    printf("writing pd_create_X functions...\n");
    int j;
    for(j=0;!feof(in);j++) { 
      
      PD_VERBOSE = 50;

      pd_code_t *inpd;

      if (weakchecking->count > 0) { 

	inpd = pd_read_weak(in);

      } else {

        inpd = pd_read(in);

      }

      assert(inpd != NULL);

      if (weakchecking->count == 0) { 

	if (!pd_ok(inpd)) {

	  pd_printf("pdcode read from file does not pass pd_ok\n",inpd);
	  exit(1);

	}

      }

      char name[4096];
      char *temp_name = mangle(infile->basename[i],".pdstor","");
      sprintf(name,"%s_%d",temp_name,j);
      free(temp_name);

      if (weakchecking->count > 0) { 
	
	pd_write_c_weak(out,inpd,name);

      } else {

	pd_write_c(out,inpd,name);
	
      }

      pd_code_free(&inpd);

      printf("\t wrote pd_create_%s\n",name);

    }
    printf("done (wrote %d functions).\n",j);

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
