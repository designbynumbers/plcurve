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

struct edge_record {
  int comp;
  int lead_vert;
  int n_struts;
};

typedef struct augmented_strut_type {
    const octrope_strut *strut;
    double compression;
    double s,t;
    double end;
} augmented_strut;

/* typedef struct plot_box_type {
    int pos[2];
    int nstruts;    
    } plot_box; */

#define pline_edges(P) (((P).open) ? (P).nv-1 : (P).nv)

int PD_VERBOSE = 0;

int compare_augmented_struts(const void *a, const void *b)
{
  augmented_strut *A,*B;

  A = (augmented_strut *)(a);
  B = (augmented_strut *)(b);

  if (A->s < B->s) {
    return -1; 
  } else {
    return 1;
  }
}
        
int compare_edgerecs(const void *a, const void *b)
{
  struct edge_record *A,*B;

  A = (struct edge_record *)(a);
  B = (struct edge_record *)(b);

  if (A->n_struts - B->n_struts != 0) {

    return B->n_struts - A->n_struts;

  } 

  if (A->comp - B->comp != 0) {

    return A->comp - B->comp;

  }

  return A->lead_vert - B->lead_vert;

}

int compare_struts_by_length(const void *a,const void *b)
{
  octrope_strut *A,*B;
  
  A = (octrope_strut *)(a);
  B = (octrope_strut *)(b);

  if (A->length > B->length) return 1;
  return -1;
}

void analyze_and_print_struts(plCurve *L, int n_struts, 
                              octrope_strut *strutlist)

{
  struct edge_record *edgedata;
  int    *offsets;

  int    i,j;

  fprintf(stderr,"\nstruts: Building edge database.\n");

  /* We first build a table of edge records. */

  offsets    =  (int *)(calloc(L->nc,sizeof(int)));
  offsets[0] = 0;
  
  for (i=1;i<L->nc;i++) {
    offsets[i] = pline_edges(L->cp[i-1]) + offsets[i-1];
  }
    
  edgedata = (struct edge_record *)(calloc(plc_num_edges(L),
					   sizeof(struct edge_record)));

  for(i=0; i < L->nc; i++) {
    for(j=0; j < pline_edges(L->cp[i]); j++) {
      edgedata[offsets[i] + j].comp = i;
      edgedata[offsets[i] + j].lead_vert = j;
      edgedata[offsets[i] + j].n_struts = 0;
    }
  }

  /* Now we go through and count the number of struts in which each edge is *
   * involved.                                                              */

  for(i=0;i<n_struts;i++) {
    edgedata[offsets[strutlist[i].component[0]] + strutlist[i].lead_vert[0]].n_struts++;
    edgedata[offsets[strutlist[i].component[1]] + strutlist[i].lead_vert[1]].n_struts++;
  }

  /* Sort by that count ... */
  fprintf(stderr,"struts: Sorting %d edges by number of incident struts (%d struts total)...\n\n",
          plc_num_edges(L),n_struts);

  qsort(edgedata,plc_num_edges(L),sizeof(struct edge_record),compare_edgerecs);

  /* And print out the results. */

  for(i=0; i< min(plc_num_edges(L),10);i++) {
    printf("Comp: %d  Lead Vert: %d  # of struts: %d.\n",
           edgedata[i].comp,edgedata[i].lead_vert,edgedata[i].n_struts);
  }

  /* Now we free some of our memory. */

  free(offsets);
  free(edgedata);

  /* Now we display some information about the shortest struts. */

  fprintf(stderr,"\n\nstruts: Sorting %d struts by strut length...\n\n",n_struts);

  qsort(strutlist,n_struts,sizeof(octrope_strut),compare_struts_by_length);

  for(i=0; i < min(n_struts,10);i++) {
    printf("%10g %3d:%3d  (%5g) - %3d:%3d (%5g)\n",
           strutlist[i].length,
	   strutlist[i].component[0], strutlist[i].lead_vert[0], strutlist[i].position[0],
	   strutlist[i].component[1], strutlist[i].lead_vert[1], strutlist[i].position[1]);
  }

}

int main(int argc,char *argv[]) {

  FILE *infile_fptr,*outfile_fptr;
  plCurve *L;

  int    sl_size,n_struts,s_cnt;
  octrope_strut *strutlist;
  double shortest;
  double eps = 0.01;
  double cutoff = 1.0;

  int nerrors;
  char outfile_name[1000];
  
  struct arg_file *infile  = arg_file1(NULL,NULL,"<file>", "input file");
  struct arg_int  *levels  = arg_int0("l", "levels","<n>", "number of octree levels");
  struct arg_int  *debuglevel = arg_int0("v", "verbosity","<n>", "level of debugging information to print");
  struct arg_lit  *stats   = arg_lit0("s","stats","display strut set statistics");
  struct arg_dbl  *epsilon = arg_dbl0("e","epsilon","<x>","find struts within this tolerance of minimum length");
  struct arg_dbl  *tube_radius = arg_dbl0("r","radius","<x>","find struts shorter than 2*x");
  struct arg_lit  *show    = arg_lit0("S","show","Show struts");
  struct arg_int  *maxst   = arg_int0("M","maxstruts","<n>","maximum number of struts to find");
  struct arg_lit  *strut_v = arg_lit0(NULL,"sv","Show struts by their endpoints.  Implies -n.");
  struct arg_lit  *help    = arg_lit0("h","help","display help message");
  struct arg_lit  *no_file = arg_lit0("n","nofile","Don't write out the strut file");
  struct arg_end  *end     = arg_end(20);

  void *argtable[] = {help,epsilon,tube_radius,levels,stats,show,maxst,
                      no_file,strut_v,debuglevel,infile,end};
  
  struct arg_end  *helpend = arg_end(20);

  void *helptable[] = {help, helpend};
  char revision[20] = "$Revision: 1.18 $";
  char *dollar;

  plc_vector se[2];
  int i;
    
  dollar = strchr(&revision[1],'$');
  dollar[0] = '\0';
  printf("Struts v%s, liboctrope v%s\n",&revision[11],PACKAGE_VERSION);
  /* The PACKAGE_VERSION preprocessor symbol is defined in "config.h" by autoconf. */

  printf("  Find the struts for a polygonal knot.\n");
  
  /* We start by parsing the command-line arguments with argtable. */
  
  if (arg_nullcheck(argtable) != 0 || arg_nullcheck(helptable) != 0)
    printf("error: insufficient memory\n");

  nerrors = arg_parse(argc,argv,argtable);
  
  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */

    nerrors = arg_parse(argc,argv,helptable);

    if (nerrors > 0) {  /* The help table didn't match either-- the first set of */
                        /* errors was probably more helpful, so we display it. */

      arg_print_errors(stdout,helpend,"helptable");
      arg_print_errors(stdout,end,"struts");
      exit(1);

    } else {  /* The help table matched, which means we asked for help or gave nothing */
  
      printf("struts computes the set of self-contacts of a Geomview VECT file.\n"
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

  if (tube_radius->count > 0) {

    cutoff = 2*tube_radius->dval[0];

  }

  /* Next, we mess around with opening files and such. */
  /* Notice that "argtable" will have already objected if no filename was passed */
  
  infile_fptr = fopen(*(infile->filename),"r");
  
  if (infile_fptr == NULL) {
    
    fprintf(stderr,"struts: Couldn't open file %s.\n",*(infile->filename));
    exit(0);

  }
  
  octrope_error_num = 0;
  L = plc_read(infile_fptr,&octrope_error_num,octrope_error_str,80);
  
  /* We now demonstrate the octrope library's error handling protocol: */
  
  if (octrope_error_num > 0) {   /* This is the signal for an error. */
    
    fprintf(stderr,"struts: link reading error\n%s\n",octrope_error_str);
    exit(octrope_error_num);
    
  }
  
  fclose(infile_fptr);
  
  /* We now compute the strut set. */
  
  if (maxst->count > 0) {
    sl_size = *(maxst->ival);
  } else {
    /* We allocate an extra-large buffer. */
    sl_size = 20*plc_num_edges(L);      
  }
  strutlist = (octrope_strut *)(calloc(sl_size,sizeof(octrope_strut)));

  if (epsilon->count > 0) {

    n_struts = octrope_struts(L,0,eps,strutlist,sl_size,&shortest,NULL,0);

  } else if (tube_radius->count > 0) {

    n_struts = octrope_struts(L,cutoff,0,strutlist,sl_size,&shortest,NULL,0);

  } else {

    n_struts = octrope_struts(L,1.0,0,strutlist,sl_size,&shortest,NULL,0);

  } 
  
  if (octrope_error_num > 0) {
    
    fprintf(stderr,"%s\n",octrope_error_str);
    exit(1);
    
  }
  
  /* We now go ahead and check whether we should output human-readable data. */
  
  if (stats->count > 0) {
    analyze_and_print_struts(L,n_struts,strutlist);
  }   
  
  if (show->count > 0) {
    for (s_cnt = 0; s_cnt < n_struts; s_cnt++) {
      printf("%d:%d - %d:%d %15.15f %15.15f\n",
        strutlist[s_cnt].component[0],strutlist[s_cnt].lead_vert[0],
        strutlist[s_cnt].component[1],strutlist[s_cnt].lead_vert[1],
        strutlist[s_cnt].position[0],strutlist[s_cnt].position[1]
      );
    }
  }

  if (strut_v->count > 0) {
    for(i=0;i<n_struts;i++) {
      octrope_strut_ends(L,&strutlist[i],se);
      printf("%15.15g %15.15g %15.15g to %15.15g %15.15g %15.15g\n",
        se[0].c[0],se[0].c[1],se[0].c[2],
        se[1].c[0],se[1].c[1],se[1].c[2]);
    }
  }
  
  /* We now prepare to output the data to a SKEL file. */   
  
  if (no_file->count == 0 && strut_v->count == 0) {
  
    if (strlen(*(infile->basename)) > sizeof(outfile_name)-20) {
      
      fprintf(stderr,"struts: Ridiculously long input filename can't be parsed.\n");
      exit(1);
      
    }
    
    /* Create filename for SKEL output. */
    
    sprintf(outfile_name,"%s",*(infile->basename));
    
    if (strstr(outfile_name,".vect") != NULL) {
  
      sprintf(strstr(outfile_name,".vect"),".struts.skel");
  
    } else {
  
      strcat(outfile_name,".struts.skel");
  
    }
    
    outfile_fptr = fopen(outfile_name,"w");
    
    if (outfile_fptr == NULL) {
      
      fprintf(stderr,"struts: Couldn't open %s for writing.\n",outfile_name);
      exit(1);
      
    }
    
    /* Now we write the strut set to a Geomview "SKEL" data structure. */
    
    fprintf(outfile_fptr,"SKEL\n");
    fprintf(outfile_fptr,"%d %d \n",2*n_struts,n_struts);
    
    for(i=0;i<n_struts;i++) {
      
      octrope_strut_ends(L,&strutlist[i],se);
      fprintf(outfile_fptr,"%g %g %g \n",se[0].c[0],se[0].c[1],se[0].c[2]);
      fprintf(outfile_fptr,"%g %g %g \n",se[1].c[0],se[1].c[1],se[1].c[2]);
      
    }
    
    for(i=0;i<n_struts;i++) {
      
      fprintf(outfile_fptr,"2 %d %d \n",2*i,2*i+1);
      
    }
    
    /* Now we close files and free memory. */
    
    printf("struts: %d struts written to %s.\n",n_struts,outfile_name);
    fclose(outfile_fptr);   
      
  } else {
    printf("struts: %d struts found.\n",n_struts);
  }
  free(strutlist);
  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

  exit(0);
}
