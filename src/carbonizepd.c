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

pd_code_t *pd_read_err_weak(FILE *infile, int *err)

/* Reads an (ASCII) pd code written by pd_write. Return a pointer if we succeed, NULL if we fail. */

{
  /* pd   <hash> <uid> */ /* Remember that the hash is base64 encoded using a-zA-Z;: character set */

  unsigned long int input_temp,input_temp2,input_temp3,input_temp4;
  char hash[PD_HASHSIZE];
  pd_uid_t uid;

  /* The first thing we do is spin through whitespace until we encounter a "p" */
  int read_char;
  for(read_char = fgetc(infile);isspace(read_char);read_char = fgetc(infile));
  ungetc(read_char,infile);

  /* We now read the rest of the line, to see if it contains a hash
     and uid. We're going to bet on hashes (and hence, lines) no
     longer than 4096 characters here. */

  char pd_line[4096];

  if (fgets(pd_line,4096,infile) == NULL) {

    pd_error(SRCLOC,"infile is already at EOF-- can't pd_read from it\n",NULL);
    return NULL;

  }

  /* Now we try the alternatives. If we're going to read a hash, we have
     to do a little dodge to construct the pattern for sscanf, because we
     want the sscanf pattern to reflect PD_HASHSIZE (which might change) */

  char hash_template[1024];
  sprintf(hash_template," pd %%%d[a-zA-Z0-9;:] %%lu ",PD_HASHSIZE);

  int sscanf_result = sscanf(pd_line,hash_template,hash,&input_temp);

  if (sscanf_result != 2) {

    /* Check for an incorrect hash. */

    char bogus_hash[4096];
    if (sscanf(pd_line," pd %4096s %lu ",bogus_hash,&input_temp) == 2) {

      pd_error(SRCLOC,
	       "first line of pdcode file\n"
	       "%s"
	       "is in format\n"
	       "pd <hash> <uid>\n"
	       "but hash == %s, which violates a hash rule\n"
	       "the hash may be at most %d characters long\n"
	       "and should only contain characters [a-zA-Z0-9;:]\n"
	       "\n"
	       "If you don't know the hash, it should be omitted\n"
	       "and this line should read\n"
	       "\n"
	       "pd\n",NULL,pd_line,bogus_hash,PD_HASHSIZE);
      return NULL;

    } else if (sscanf(pd_line," pd %4096s %lu ",bogus_hash,&input_temp) == 1) {

      pd_error(SRCLOC,
	       "first line of pdcode file appears to be in format\n"
	       "pd <hash> <uid>\n"
	       "and includes hash %s, but does not include uid (an unsigned integer)\n",NULL,bogus_hash);
      return NULL;

    }

    /* Well, it looks like we didn't even TRY to provide a hash and uid. */
    /* In this case, we should match "(whitespace)pd(whitespace)" */

    char pd_string[32];

    if (sscanf(pd_line," %32s ",pd_string) != 1) {

      pd_error(SRCLOC,
	       "first line of pdcode file is neither in format\n"
	       "pd <hash> <uid>\n"
	       "nor\n"
	       "pd\n",NULL);
      return NULL;

    }

    if (strcmp(pd_string,"pd") != 0) {

      pd_error(SRCLOC,
	       "first line of pdcode file contains\n"
	       "%s\n"
	       "which is neither in format\n"
	       "pd <hash> <uid>\n"
	       "nor\n"
	       "pd\n",NULL,pd_string);
      return NULL;

    }

    /* We DO actually match this. This means that we should set the
       hash and the uid ourselves. */

    uid = PD_UNSET_UID;
    sprintf(hash,"unset");

  } else { /* We matched the pd <hash> <uid> format */

    uid = (pd_uid_t)(input_temp);

  }

  /* nv   <nverts> */

  int cross,edge,comp,pos,face;
  if (fscanf(infile," nv %lu ",&input_temp) != 1) {

    pd_error(SRCLOC,
	     "second line of pdcode file should be in format\n"
	     "nv <nverts>\n"
	     "where <nverts> == number of crossings in pd_code\n",NULL);

    return NULL;

  }

  /* We now know the number of crossings to allocate space for,
     and we can call pd_new. We make sure that we have space for
     at least two more crossings than we see now (to prevent zero
     crossing diagrams from crashing everything). */

  if ((pd_idx_t)(input_temp) == 0) {

    /* If we're a zero-crossing unknot, there is one piece of 
       information that can actually matter, which is the tag 
       of the (single) component. We're going to discard the 
       rest of the file, but first we'll run through and extract
       the "tag" (if present). */

    char rolling_buffer[4];
    pd_tag_t tag = 'A';
    
    rolling_buffer[0] = (char)(getc(infile));
    rolling_buffer[1] = (char)(getc(infile));
    rolling_buffer[2] = (char)(getc(infile));
    rolling_buffer[3] = 0;

    for(;!feof(infile);) {

      if (!strcmp(rolling_buffer,"tag")) { 

	fscanf(infile," %c ",&tag); 

      }

      rolling_buffer[0] = rolling_buffer[1]; rolling_buffer[1] = rolling_buffer[2];
      rolling_buffer[2] = (char)(getc(infile));

    }

    pd_code_t *outcode;
    outcode = pd_build_unknot(0);
    outcode->comp[0].tag = tag;
      
    return outcode;

  }

  pd_code_t *pd = pd_code_new((pd_idx_t)(input_temp));

  /* We now have to copy the hash into the new pd code. */

  strncpy(pd->hash,hash,PD_HASHSIZE);
  pd->uid = uid;
  pd->ncross = (pd_idx_t)(input_temp);

  if (pd->ncross > pd->MAXVERTS) {

    fprintf(stderr,
	    "%s (%d): Reading pd code which appears to be valid but has %d crossings. "
	    "         We allocated space for pd->MAXVERTS = %d.\n",__FILE__,__LINE__,
	    pd->ncross,pd->MAXVERTS);
    exit(1);

  }

  /* <nv lines of crossing information in the format edge edge edge edge> */
  /* There's an optional crossing sign at the end of the line (+ or -) */

  for(cross=0;cross<pd->ncross;cross++) {

    for(pos=0;pos<4;pos++) {

      if(fscanf(infile," %lu ",&input_temp) != 1) {

	pd_error(SRCLOC,"error on crossing %d (of %d), in pd (so far) %PD",pd,cross,pd->ncross);
	pd_code_free(&pd);
	return NULL;

      }

      pd->cross[cross].edge[pos] = (pd_idx_t)(input_temp);

    }

    int peek;
    peek = fgetc(infile);
    if (peek == '+') {

      pd->cross[cross].sign = PD_POS_ORIENTATION;

    } else if (peek == '-') {

      pd->cross[cross].sign = PD_NEG_ORIENTATION;

    } else {

      ungetc(peek,infile);
      pd->cross[cross].sign = PD_UNSET_ORIENTATION;

    }

  }

  /* ne   <nedges> */

  if (fscanf(infile," ne %lu ",&input_temp) != 1) {

    pd_error(SRCLOC,
	     "after crossing lines, should have\n"
	     "ne <nedges>\n"
	     "but this file does not\n",pd);
    pd_code_free(&pd);
    return NULL;

  }
  pd->nedges = (pd_idx_t)(input_temp);

  if (pd->nedges > pd->MAXEDGES) {

    pd_error(SRCLOC,
	    "%s (%d): Reading pd code which appears to be valid but has %d edges. "
	    "         We allocated space for pd->MAXEDGES = %d.\n",
	    pd,pd->nedges,pd->MAXEDGES);
    pd_code_free(&pd);
    return NULL;

  }

  /* <ne lines of crossing information in the format tail, tailpos -> head, headpos> */

  for(edge=0;edge<pd->nedges;edge++) {

    if(fscanf(infile," %lu,%lu -> %lu,%lu ",
	      &input_temp,&input_temp2,&input_temp3,&input_temp4) != 4) {

      pd_error(SRCLOC,
	       "edge %d is not in the format\n"
	       "<crossing>,<pos> -> <crossing>,<pos>\n"
	       "where <crossing> and <pos> are positive integers\n",pd,edge);
      pd_code_free(&pd); return NULL;

    }

    pd->edge[edge].tail = (pd_idx_t)(input_temp);
    pd->edge[edge].tailpos = (pd_pos_t)(input_temp2);

    pd->edge[edge].head = (pd_idx_t)(input_temp3);
    pd->edge[edge].headpos = (pd_pos_t)(input_temp4);

  }

  /* nc   <ncomps> */

  if (fscanf(infile,"nc %lu ",&input_temp) != 1) {

    pd_error(SRCLOC,
	     "after edge lines, pdcode file should have a line\n"
	     "nc <ncomps>\n"
	     "where <ncomps> is a positive integer\n"
	     "but this one doesn't",pd);
    pd_code_free(&pd); return NULL;

  }
  pd->ncomps = (pd_idx_t)(input_temp);

  if (pd->ncomps > pd->MAXCOMPONENTS) {

    pd_error(SRCLOC,
	    "Reading pd code which appears to be valid but has %d components. "
	    "We allocated space for pd->MAXCOMPONENTS = %d in %PD\n",
	    pd,pd->ncomps,pd->MAXCOMPONENTS);
    pd_code_free(&pd);
    return NULL;

  }

  /* <nc lines, each containing nedges : +/- edge +/- edge .... +/- edge> */

  pd_tag_t next_tag = 'A';

  for(comp=0;comp<pd->ncomps;comp++) {

    if (fscanf(infile," %lu : ",&input_temp) != 1) {

      pd_error(SRCLOC,
	       "component %d should start with\n"
	       "<nedges> : \n"
	       "where <nedges> is a positive integer\n"
	       "but this file doesn't\n",pd,comp);
      pd_code_free(&pd);
      return false;

    }

    pd->comp[comp].nedges = (pd_idx_t)(input_temp);

    /* We haven't allocated space in the component. */

    pd->comp[comp].edge = calloc(pd->comp[comp].nedges,sizeof(pd_idx_t));
    assert(pd->comp[comp].edge != NULL);

    if (pd->comp[comp].nedges > pd->MAXEDGES) {

      fprintf(stderr,
	    "%s (%d): Reading component which appears to be valid but has %d edges. "
	    "         We expect only a maximum of pd->MAXEDGES = %d.\n",
	      __FILE__,__LINE__,pd->nedges,pd->MAXEDGES);
      exit(1);

    }

    for(edge=0;edge<pd->comp[comp].nedges;edge++) {

      if(fscanf(infile," %lu ",&input_temp) != 1) {

	pd_error(SRCLOC,"edge %d of component %d should be a positive integer\n"
		 "but isn't in this file\n",pd,edge,comp);
	pd_code_free(&pd); return NULL;

      }

      pd->comp[comp].edge[edge] = (pd_idx_t)(input_temp);


    }

    /* Now we're at the end, and the next character is either "tag X" or a newline. */
    /* Skip whitespace until we find a character... */

    char testchar;
    fscanf(infile," %c",&testchar);

    if (testchar == 't') { /* for tag */

      ungetc((int)(testchar),infile);
      if (fscanf(infile," tag %c ",&(pd->comp[comp].tag)) != 1) {

	pd_error(SRCLOC,
		 "component %d of pd_code has extra characters after %d edges\n"
		 "which start with t (so we're assuming a tag) but they don't match\n"
		 "tag <tag>\n"
		 "where <tag> is a single ASCII character\n",
		 pd,comp,edge);
	pd_code_free(&pd);
	return NULL;
      }

    } else {

      ungetc((int)(testchar),infile);
      pd->comp[comp].tag = next_tag++;

    }

  }

  /*nf   <nfaces> */

  if (fscanf(infile," nf %lu ",&input_temp) != 1) {

    pd_error(SRCLOC,
	     "after component data, pdfile should contain\n"
	     "nf <nfaces>\n"
	     "where <nfaces> is a positive integer\n"
	     "but this one doesn't\n",pd);
    pd_code_free(&pd); return NULL;

  }
  pd->nfaces = (pd_idx_t)(input_temp);

  if (pd->nfaces > pd->MAXFACES) {

    pd_error(SRCLOC,
	    "Reading pd code which appears to be valid but has %d faces. "
	     "We expected a maximum of pd->MAXFACES = %d in %PD\n",pd,
	     pd->nfaces,pd->MAXFACES);
    return NULL;

  }

  /* <nf lines, each in the format nedges edge edge
     ... edge giving face info counterclockwise> */

  for(face=0;face<pd->nfaces;face++) {

    if (fscanf(infile," %lu : ",&input_temp) != 1) {

      pd_error(SRCLOC,
	       "face record %d (of %d) is expected to begin\n"
	       "<nedges> : \n"
	       "where <nedges> is a positive integer\n"
	       "but this one doesn't in %PD\n",pd,
	       face,pd->nfaces);

      pd_code_free(&pd);
      return NULL;

    }

    pd->face[face].nedges = (pd_idx_t)(input_temp);

    if (pd->face[face].nedges > pd->MAXEDGES) {

      pd_error(SRCLOC,
	       "Reading face which appears to be valid but has %d edges. "
	       "We expected a maximum of pd->MAXEDGES = %d.\n",
	       pd,pd->face[face].nedges,pd->MAXEDGES);
      return NULL;

    }

    pd->face[face].edge = calloc(pd->face[face].nedges,sizeof(pd_idx_t));
    pd->face[face].or = calloc(pd->face[face].nedges,sizeof(pd_or_t));
    assert(pd->face[face].edge != NULL && pd->face[face].or != NULL);

    for(edge=0;edge<pd->face[face].nedges;edge++) {

      char orientation[2];
      if(fscanf(infile," %1[+-]s ",orientation) != 1) {

	pd_error(SRCLOC,
		 "edge %d of face %d in pdcode file\n"
		 "doesn't have a sign which matches [+-]\n",
		 pd,edge,face);
	pd_code_free(&pd);
	return NULL;

      }

      if (fscanf(infile," %lu ",&input_temp) != 1) {

	pd_error(SRCLOC,
		 "edge %d of face %d in pdcode file\n"
		 "has sign (+/-) but doesn't match\n"
		 "<+-> <edge>\n"
		 "where <edge> is a positive integer\n",
		 pd,edge,face);

	pd_code_free(&pd); return NULL;
      }

      pd->face[face].edge[edge] = (pd_idx_t)(input_temp);
      pd->face[face].or[edge] = (orientation[0] == '+') ? PD_POS_ORIENTATION : PD_NEG_ORIENTATION;

    }

  }

  return pd;

}

pd_code_t *pd_read_weak(FILE *infile)
{
  return pd_read_err_weak(infile,NULL);
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
