/*
 * struts.c : This test program reads link information from a VECT file and
 *            outputs the strut information in various (hopefully useful)
 *            formats. 
 */

/* Copyright 2004 The University of Georgia. */

/* This file is part of liboctrope.
   
liboctrope is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

liboctrope is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with liboctrope; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#ifdef HAVE_MATH_H
#include <math.h>
#endif
#ifdef HAVE_STRING_H
#include <string.h>
#endif

#include "octrope.h"
#include "octrope_internal.h"
#include"argtable2.h"

#ifndef min 
  #define min(A,B) ((A < B) ? A : B)
#endif

int compare_mrlocs(const void *a, const void *b)
{
  octrope_mrloc *A,*B;

  A = (octrope_mrloc *)(a);
  B = (octrope_mrloc *)(b);

  if (A->mr - B->mr > 0) {

    return 1;

  } else {

    return -1;

  }

}


void analyze_and_print_mrlocs(plCurve *L,int n_mr,octrope_mrloc *mrlist)

{

  printf("mrlocs: Sorting %d mrlocs by minrad...",n_mr);
  qsort(mrlist,n_mr,sizeof(octrope_mrloc),compare_mrlocs);
  printf("done.\n\n");

  int i;

  printf("Cmp Vert Minrad \n"
	 "-------------------------\n");

  for(i=0;i<n_mr;i++) {

    printf("%3d %4d %.16g\n",mrlist[i].component,mrlist[i].vert,mrlist[i].mr);

  }

  printf("\n\n");

}
    

int main(int argc,char *argv[]) 

{

  FILE *infile_fptr,*outfile_fptr;
  plCurve *L;

  int    mr_size,n_mr;
  octrope_mrloc *mrlist;
  double minrad;
  double eps = 0;
  double cutoff = 0.5;
  
  int nerrors;
  char outfile_name[1000];
  
  struct arg_file *infile  = arg_file1(NULL,NULL,"<file>", "input file");
  struct arg_int  *levels  = arg_int0("l", "levels","<n>", "number of octree levels");
  struct arg_int  *debuglevel = arg_int0("v", "verbosity","<n>", "level of debugging information to print");
  struct arg_lit  *stats   = arg_lit0("s","stats","display minrad loc set statistics");
  struct arg_dbl  *epsilon = arg_dbl0("e","epsilon","<x>","find mrlocs within this tolerance of minimum length");
  struct arg_dbl  *lambda = arg_dbl0("l","lambda","<x>","find mrlocs with radius of curvature < x");
  struct arg_lit  *show    = arg_lit0("S","show","Show mrlocs");
  struct arg_int  *maxmr   = arg_int0("M","maxmr","<n>","maximum number of mrlocs to find");
  
  struct arg_lit  *help    = arg_lit0("h","help","display help message");
  struct arg_lit  *no_file = arg_lit0("n","nofile","Don't write out the mrloc file");
  struct arg_end  *end     = arg_end(20);

  void *argtable[] = {help,epsilon,lambda,levels,stats,show,maxmr,
                      no_file,debuglevel,infile,end};
  
  struct arg_end  *helpend = arg_end(20);

  void *helptable[] = {help, helpend};
  char revision[20] = "$Revision: 1.18 $";
  char *dollar;

  int i;
    
  dollar = strchr(&revision[1],'$');
  dollar[0] = '\0';
  printf("Mrloc v%s, liboctrope v%s\n",&revision[11],PACKAGE_VERSION);
  /* The PACKAGE_VERSION preprocessor symbol is defined in "config.h" by autoconf. */

  printf("  Find the minrad locations for a polygonal knot.\n");
  
  /* We start by parsing the command-line arguments with argtable. */
  
  if (arg_nullcheck(argtable) != 0 || arg_nullcheck(helptable) != 0)
    printf("error: insufficient memory\n");

  nerrors = arg_parse(argc,argv,argtable);
  
  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */

    nerrors = arg_parse(argc,argv,helptable);

    if (nerrors > 0) {  /* The help table didn't match either-- the first set of */
                        /* errors was probably more helpful, so we display it. */

      arg_print_errors(stdout,helpend,"helptable");
      arg_print_errors(stdout,end,"mrlocs");
      exit(1);

    } else {  /* The help table matched, which means we asked for help or gave nothing */
  
      printf("mrlocs computes the set of minimum radius corners of a Geomview VECT file.\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }
  
  /* We have now parsed the command-line args, and we slot them into the proper globals. */
  
  if (levels->count > 0) {
    
    octrope_set_levels(*(levels->ival));
    
  }
  
  if (debuglevel->count > 0) {
    
    octrope_set_debug(*(debuglevel->ival));
    printf("Debug level set to %d.\n",*(debuglevel->ival));
    
  }

  if (epsilon->count > 0) {

    eps = *(epsilon->dval);

  }

  if (lambda->count > 0) {

    cutoff = *(lambda->dval)/2.0;

  }
  
  /* Next, we mess around with opening files and such. */
  /* Notice that "argtable" will have already objected if no filename was passed */
  
  infile_fptr = fopen(*(infile->filename),"r");
  
  if (infile_fptr == NULL) {
    
    fprintf(stderr,"mrlocs: Couldn't open file %s.\n",*(infile->filename));
    exit(0);

  }
  
  octrope_error_num = 0;
  L = plc_read(infile_fptr,&octrope_error_num,octrope_error_str,80);
  
  /* We now demonstrate the octrope library's error handling protocol: */
  
  if (octrope_error_num > 0) {   /* This is the signal for an error. */
    
    fprintf(stderr,"mrlocs: link reading error\n%s\n",octrope_error_str);
    exit(octrope_error_num);
    
  }
  
  fclose(infile_fptr);
  
  /* We now compute the mr set. */
  
  if (maxmr->count > 0) {
    mr_size = *(maxmr->ival);
  } else {
    /* We allocate an extra-large buffer. */
    mr_size = 2*plc_num_edges(L)+1;      
  }
  mrlist = (octrope_mrloc *)(calloc(mr_size,sizeof(octrope_mrloc)));

  if (epsilon->count > 0) {

    minrad = octrope_minrad(L,0,eps,mrlist,mr_size,&n_mr);

  } else {

     minrad = octrope_minrad(L,cutoff,0,mrlist,mr_size,&n_mr);

  }
  
  if (octrope_error_num > 0) {
    
    fprintf(stderr,"%s\n",octrope_error_str);
    exit(1);
    
  }
  
  /* We now go ahead and check whether we should output human-readable data. */
  
  if (stats->count > 0) {
    analyze_and_print_mrlocs(L,n_mr,mrlist);
  }   
  
  /* We now prepare to output the data to a file. */   
  
  if (no_file->count == 0) {
  
    if (strlen(*(infile->basename)) > sizeof(outfile_name)-20) {
      
      fprintf(stderr,"mrlocs: Ridiculously long input filename can't be parsed.\n");
      exit(1);
      
    }
    
    /* Create filename for output. */
    
    sprintf(outfile_name,"%s",*(infile->basename));
    
    if (strstr(outfile_name,".vect") != NULL) {
  
      sprintf(strstr(outfile_name,".vect"),".mrlocs");
  
    } else {
  
      strcat(outfile_name,".mrlocs");
  
    }
    
    outfile_fptr = fopen(outfile_name,"w");
    
    if (outfile_fptr == NULL) {
      
      fprintf(stderr,"mrlocs: Couldn't open %s for writing.\n",outfile_name);
      exit(1);
      
    }
    
    /* Now we write the mrloc set to a file. */

    fprintf(outfile_fptr,"MRLOCS\n");
    
    for(i=0;i<n_mr;i++) {

      fprintf(outfile_fptr,"%d %d %.16g\n",
	      mrlist[i].component,mrlist[i].vert,mrlist[i].mr);

    }
    
    /* Now we close files and free memory. */
    
    printf("mrlocs: %d mrlocs written to %s.\n",n_mr,outfile_name);
    fclose(outfile_fptr);   
      
  } else {
    printf("mrlocs: %d mrlocs found.\n",n_mr);
    printf("mrlocs: minrad %g\n",minrad);
  }
  free(mrlist);
  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

  exit(0);
}
