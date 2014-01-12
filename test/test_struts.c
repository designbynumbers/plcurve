/*
 * test_struts.c : This test program checks the "octrope_struts" call.
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


struct edge_record {
  int comp;
  int lead_vert;
  int n_struts;
};

int TESTS_PASSED = 0;

#ifndef min 
  #define min(A,B) ((A < B) ? A : B)
#endif

/* Forward function declarations. */

void analyze_and_print_struts(plCurve *L, int n_struts, 
                              octrope_strut *strutlist);

/**********************  Program Testing Code *****************************/

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


void test_case(plCurve *L,octrope_strut expected_strut)

{
  octrope_strut strutlist[20];
  int n_struts;
  double shortest;

  n_struts = octrope_struts(L,0,100,strutlist,sizeof(strutlist),&shortest,
                            NULL,0);

  if (octrope_error_num > 0) {  /* There was an error in the octrope_struts call. */

    fprintf(stderr,"test_struts: octrope_struts reports %s.\n",octrope_error_str);
    exit(1);

  }

  if (n_struts == 1) {

    if ((strutlist[0].component[0] == expected_strut.component[0] &&
        strutlist[0].lead_vert[0] == expected_strut.lead_vert[0] &&
         fabs(strutlist[0].position[0] - expected_strut.position[0]) < 1e-10) 
        &&
        (strutlist[0].component[1] == expected_strut.component[1] &&
         strutlist[0].lead_vert[1] == expected_strut.lead_vert[1] &&
         fabs(strutlist[0].position[1] - expected_strut.position[1]) < 1e-10)) 
      
      {

        TESTS_PASSED++;
        return;
        
      }

    if ((strutlist[0].component[1] == expected_strut.component[0] &&
        strutlist[0].lead_vert[1] == expected_strut.lead_vert[0] &&
        fabs(strutlist[0].position[1] - expected_strut.position[0]) < 1e-10) 
        &&
        (strutlist[0].component[0] == expected_strut.component[1] &&
        strutlist[0].lead_vert[0] == expected_strut.lead_vert[1] &&
         fabs(strutlist[0].position[0] - expected_strut.position[1]) < 1e-10)) 

      {

        TESTS_PASSED++;
        return;
        
      }


  }
    
  /* We failed the "struts" test. Go ahead and print the results. */

  fprintf(stderr,"Strut test %d FAILED.\n\n",TESTS_PASSED);

  fprintf(stderr,"Expected strut from (%d, %d, %g) to (%d, %d, %g) "
	  "(comp,vert,pos).\n",
          expected_strut.component[0],
	  expected_strut.lead_vert[0],expected_strut.position[0],
          expected_strut.component[1],expected_strut.lead_vert[1],
	  expected_strut.position[1]);

  fprintf(stderr,"Now printing list of found struts...\n");
  analyze_and_print_struts(L,n_struts,strutlist);

  exit(1);

}

#define pline_edges(P) (((P).open) ? (P).nv-1 : (P).nv)

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
  fprintf(stderr,"struts: Sorting %d edges by number of incident "
      "struts (%d struts total)...\n\n", plc_num_edges(L),n_struts);

  qsort(edgedata,plc_num_edges(L),sizeof(struct edge_record),compare_edgerecs);

  /* And print out the results. */

  for(i=0; i< min(plc_num_edges(L),10);i++) {
    printf("Comp: %d  Lead Vert: %d  # of struts: %d.\n",
           edgedata[i].comp,edgedata[i].lead_vert,edgedata[i].n_struts);
  }

  /* Now we free some of our memory. */

  free(offsets);
  free(edgedata);
}


int main() 

     /* Performs a set of basic confidence checks on the strut finder. */

{
  plCurve  *L;
  octrope_strut exp_strut;
  int nv[2] = {3,3}, ccarray[2] = {0,0}; 
  bool open[2] = {true,true};

  L = plc_new(2,nv,open,ccarray);

  /* Now we go ahead and make the tests, moving the verts of each component to produce struts. */

  /* This should be an edge-edge strut. */

  L->cp[0].vt[0].c[0] = -2;    L->cp[0].vt[1].c[0] = -1;    L->cp[0].vt[2].c[0] =  1;
  L->cp[0].vt[0].c[1] = -2;    L->cp[0].vt[1].c[1] = -1;    L->cp[0].vt[2].c[1] =  1;
  L->cp[0].vt[0].c[2] =  1;    L->cp[0].vt[1].c[2] =  1;    L->cp[0].vt[2].c[2] =  1;
 
  L->cp[1].vt[0].c[0] = -2;    L->cp[1].vt[1].c[0] = -1;    L->cp[1].vt[2].c[0] =  1;
  L->cp[1].vt[0].c[1] =  2;    L->cp[1].vt[1].c[1] =  1;    L->cp[1].vt[2].c[1] = -1;
  L->cp[1].vt[0].c[2] =  0;    L->cp[1].vt[1].c[2] =  0;    L->cp[1].vt[2].c[2] =  0;

  exp_strut.component[0] = 0; exp_strut.lead_vert[0] = 1; exp_strut.position[0] = 0.5;
  exp_strut.component[1] = 1; exp_strut.lead_vert[1] = 1; exp_strut.position[1] = 0.5;

  test_case(L,exp_strut);

  /* This should be an vertex-vertex strut. */

  L->cp[0].vt[0].c[0] = -2;    L->cp[0].vt[1].c[0] = 0;    L->cp[0].vt[2].c[0] =  2;
  L->cp[0].vt[0].c[1] = -2;    L->cp[0].vt[1].c[1] = 0;    L->cp[0].vt[2].c[1] =  2;
  L->cp[0].vt[0].c[2] =  2;    L->cp[0].vt[1].c[2] = 1;    L->cp[0].vt[2].c[2] =  2;
 
  L->cp[1].vt[0].c[0] = -2;    L->cp[1].vt[1].c[0] = 0;    L->cp[1].vt[2].c[0] =  2;
  L->cp[1].vt[0].c[1] =  2;    L->cp[1].vt[1].c[1] = 0;    L->cp[1].vt[2].c[1] = -2;
  L->cp[1].vt[0].c[2] =  -1;    L->cp[1].vt[1].c[2] = 0;    L->cp[1].vt[2].c[2] =  -1;

  exp_strut.component[0] = 0; exp_strut.lead_vert[0] = 1; exp_strut.position[0] = 0;
  exp_strut.component[1] = 1; exp_strut.lead_vert[1] = 1; exp_strut.position[1] = 0;

  test_case(L,exp_strut);

  /* This should be an vertex-edge strut. */

  L->cp[0].vt[0].c[0] = -2;    L->cp[0].vt[1].c[0] = 0;    L->cp[0].vt[2].c[0] =  2;
  L->cp[0].vt[0].c[1] = -2;    L->cp[0].vt[1].c[1] = 0;    L->cp[0].vt[2].c[1] =  2;
  L->cp[0].vt[0].c[2] =  2;    L->cp[0].vt[1].c[2] = 1;    L->cp[0].vt[2].c[2] =  2;
 
  L->cp[1].vt[0].c[0] = -2;    L->cp[1].vt[1].c[0] =-1;    L->cp[1].vt[2].c[0] =  1;
  L->cp[1].vt[0].c[1] =  0;    L->cp[1].vt[1].c[1] = 0;    L->cp[1].vt[2].c[1] =  0;
  L->cp[1].vt[0].c[2] =  0;    L->cp[1].vt[1].c[2] = 0;    L->cp[1].vt[2].c[2] =  0;

  exp_strut.component[0] = 0; exp_strut.lead_vert[0] = 1; exp_strut.position[0] = 0;
  exp_strut.component[1] = 1; exp_strut.lead_vert[1] = 1; exp_strut.position[1] = 0.5;

  test_case(L,exp_strut);

  /* This should be an edge-vertex strut. */

  L->cp[0].vt[0].c[0] = -2;    L->cp[0].vt[1].c[0] =-1;    L->cp[0].vt[2].c[0] =  1;
  L->cp[0].vt[0].c[1] =  0;    L->cp[0].vt[1].c[1] = 0;    L->cp[0].vt[2].c[1] =  0;
  L->cp[0].vt[0].c[2] =  0;    L->cp[0].vt[1].c[2] = 0;    L->cp[0].vt[2].c[2] =  0;

  L->cp[1].vt[0].c[0] = -2;    L->cp[1].vt[1].c[0] = 0;    L->cp[1].vt[2].c[0] =  2;
  L->cp[1].vt[0].c[1] = -2;    L->cp[1].vt[1].c[1] = 0;    L->cp[1].vt[2].c[1] =  2;
  L->cp[1].vt[0].c[2] =  2;    L->cp[1].vt[1].c[2] = 1;    L->cp[1].vt[2].c[2] =  2;

  exp_strut.component[0] = 0; exp_strut.lead_vert[0] = 1; exp_strut.position[0] = 0.5;
  exp_strut.component[1] = 1; exp_strut.lead_vert[1] = 1; exp_strut.position[1] = 0;

  test_case(L,exp_strut);

  /* We now start to play games with the closing edges. */

  L->cp[0].open = false; L->cp[1].open = false;

  /* This should be a closing edge-closing edge strut. */
  
  L->cp[0].vt[0].c[0] = -2;    L->cp[0].vt[1].c[0] = 0;     L->cp[0].vt[2].c[0] =  2;
  L->cp[0].vt[0].c[1] = -2;    L->cp[0].vt[1].c[1] = 0;     L->cp[0].vt[2].c[1] =  2;
  L->cp[0].vt[0].c[2] =  1;    L->cp[0].vt[1].c[2] = 21;    L->cp[0].vt[2].c[2] = 1;

  L->cp[1].vt[0].c[0] = -2;    L->cp[1].vt[1].c[0] =  0;    L->cp[1].vt[2].c[0] =  2;
  L->cp[1].vt[0].c[1] =  2;    L->cp[1].vt[1].c[1] =  0;    L->cp[1].vt[2].c[1] =  -2;
  L->cp[1].vt[0].c[2] =  0;    L->cp[1].vt[1].c[2] =  -20;  L->cp[1].vt[2].c[2] =  0;

  exp_strut.component[0] = 0; exp_strut.lead_vert[0] = 2; exp_strut.position[0] = 0.5;
  exp_strut.component[1] = 1; exp_strut.lead_vert[1] = 2; exp_strut.position[1] = 0.5;

  test_case(L,exp_strut);

  /* We now introduce a closed curve and look for self-struts. */

  plc_free(L);

  nv[0] = 4; open[0] = 0;
  L = plc_new(1,nv,open,ccarray);

  /* A closed triangle with 3 vertices. */

  L->cp[0].vt[0].c[0] =  0;    L->cp[0].vt[1].c[0] = 10;   L->cp[0].vt[2].c[0] = 2;
  L->cp[0].vt[0].c[1] =  0;    L->cp[0].vt[1].c[1] = 0;    L->cp[0].vt[2].c[1] = 2;
  L->cp[0].vt[0].c[2] =  0;    L->cp[0].vt[1].c[2] = 0;    L->cp[0].vt[2].c[2] = 0;

  octrope_strut strutlist[20];
  int n_struts;
  double shortest;

  L->cp[0].nv = 3;

  n_struts = octrope_struts(L,0,100,strutlist,sizeof(strutlist),&shortest,
                            NULL,0);

  if (n_struts != 0) {

    fprintf(stderr,"test_struts: Expected no struts on closed triangle. Got %d.\n",
	    n_struts);
    exit(1);

  }

  /* An open triangle with 4 vertices and a 0-length strut. */

  L->cp[0].vt[3].c[0] =  0;    
  L->cp[0].vt[3].c[1] =  0;    
  L->cp[0].vt[3].c[2] =  0;    

  L->cp[0].open = true;
  L->cp[0].nv = 4;

  exp_strut.component[0] = 0; exp_strut.lead_vert[0] = 0; exp_strut.position[0] = 0;
  exp_strut.component[1] = 0; exp_strut.lead_vert[1] = 2; exp_strut.position[1] = 1.0;

  test_case(L,exp_strut);

  /* Now we exit the checker. */

  exit(0);

} 
