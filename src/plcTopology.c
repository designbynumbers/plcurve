/*
 *  Routines to work with plCurves as topological objects
 *
 */

/* Copyright 2009 The University of Georgia. */

/* This file is part of plCurve.

plCurve is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

plCurve is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with plCurve; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#include <config.h>
#include <plCurve.h>
#include <homfly.h>

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#ifdef HAVE_MATH_H
  #include <math.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_STDARG_H
#include <stdarg.h>
#endif
#ifdef HAVE_CTYPE_H
#include <ctype.h>
#endif
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#ifdef HAVE_ASSERT_H
#include <assert.h>
#endif
#ifdef HAVE_FLOAT_H
#include <float.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

char *plc_lmpoly(char *ccode,int timeout); // This is in pllmpoly02.c.

/* Convert plCurve to Millett/Ewing crossing code (documented in plCurve.h) */

/* Plan: 

   Find crossings of each component with itself, then crossings 
   of pairs of components, add to a container of crossings. */

typedef struct crossing_struct { 

  // The crossing struct encodes parameter values along each 
  // component, as well as component numbers.

  int lower_cmp;
  int upper_cmp;
  double lower_s;
  double upper_s;
  
} crossing;

typedef struct crossing_container_struct { 

  // A simple array which resizes. 

  int size;
  int used;
  crossing *buf;

} crossing_container;

int crossing_compare(const void *A, const void *B) { 

  crossing *a = (crossing *)(A);
  crossing *b = (crossing *)(B);

  int acmp,bcmp; 
  double as,bs;

  if (a->lower_cmp < a->upper_cmp) { 

    acmp = a->lower_cmp; as = a->lower_s;
    
  } else {

    acmp = a->upper_cmp; as = a->upper_s;

  }

  if (b->lower_cmp < b->upper_cmp) { 

    bcmp = b->lower_cmp; as = a->lower_s;

  } else {

    bcmp = b->upper_cmp; bs = b->upper_s;

  }

  if (a->cmp != b->cmp) { 

    return acmp - bcmp;

  } else { 

    if (as < bs) { return -1; } 
    else { return +1; }

  }

}

crossing_container *crossing_container_new(int size)
{
  
  crossing_container *ret;

  ret = (crossing_container *)(calloc(1,sizeof(crossing_container)));
  if (ret == NULL) { exit(1); }
  
  ret->size = size;
  ret->used = 0;
  ret->buf = calloc(size,sizeof(crossing));
  
  if (ret->bug == NULL) { exit(1); }

  return ret;
}


void crossing_container_free(crossing_container **c) 
{
  if (c != NULL) { 

    if (*c != NULL) { 

      if ((*c)->buf != NULL) { free((*c)->buf); (*c)->buf = NULL; }
      free(*c);

    }

    *c = NULL;

  }

}

void crossing_container_add(crossing_container *c,crossing cross) { 

  if (c->size == c->used) { 

    c->size = 2*c->size;
    c->buf = realloc(c->buf,c->size*sizeof(crossing));
    if (c->buf == NULL) { exit(1); }

  }

  c->buf[c->used++] = cross;

}

/* We now have the main crossing primitive. The basic error handling
   strategy here is this. First, we'll choose a random rotation of 
   the curve. The projection direction will always be the z-axis.

   Then, we'll generate a "flat" version of the curve, and run through 
   the planar crossing algorithm to figure out where the crossings are. 
   Each crossing will either be able to be classified, or the edges
   involved will be "tagged-as-trouble". 

   If there is a crossing, and we can figure out where, we'll return true
   and fill in the fields in the crossing record c. 

   If there is a special case of some kind, we'll just tag the offending
   vertices and continue. We'll then perturb just those vertices in space
   and re-run the process, hoping to converge. 
*/

int  segment_side(plc_vector A[2],plc_vector B) { /* Classify which side of the segment B is on */

  plc_vector a0a1,a0b,cross;

  a0a1 = plc_vect_diff(A[1],A[0]);
  a0b = plc_vect_diff(b,A[0]);
  cross = plc_cross_prod(a0a1,a0b);

  if (cross.c[2] > 1e-10) { return +1; } 
  else if (cross.c[2] < -1e-10) { return -1; }
  return 0;

}

void locate_crossing(plc_vector A[2],plc_vector B[2],double *sA,double *sB,bool *ok)
{

  /* Locating the crossing point. 

       The plan is to consider the quadrilateral A[0] B[0] A[1] B[1].
       We know that there's a crossing inside this thing:


 B[1]  *---------*  A[1]
       |\       /|
       | \     / |
       |1 \   /  |
       |   \ /   |
       |  2 * X  |
       |   / \   |
       |  /   \  |
       |3/     \ |
       |/       \|
 A[0]  *---------*  B[0]


       The idea is to compute angles 1 and 3 using plc_angle and then 
       compute angle 2. Once we have all the angles of the triangle,
       we can determine the distances from A[0] and B[1] to X using 
       the law of sines. */

    double angle1, angle2, angle3;
    double pi = 3.141592653589793;
    bool   ok1,ok3;

    angle1 = plc_angle(plc_vect_diff(A[0],B[1]),plc_vect_diff(B[0],B[1]),&ok1);
    angle3 = plc_angle(plc_vect_diff(B[1],A[0]),plc_vect_diff(A[1],A[0]),&ok3);

    if (!ok1 || !ok3) { // We can't classify the crossing

      *sA = 0; *sB = 0; *ok = false;
      return;

    } 

    if (angle1 < 1e-8 || angle3 < 1e-8) { *sA = 0; *sB = 0; *ok = false; return; }

    angle2 = pi - angle1 - angle3;

    double s1,s2,s3;
    s2 = plc_distance(A[0],B[1]);
 
    if (s2 < 1e-8 || angle2 < 1e-8) { *sA = 0; *sB = 0; *ok = false; return; }

    /* Now the law of sines says that s2/sin(angle2) = s1/sin(angle1) = s3/sin(angle3),
       so we can solve for the unknown s1 and s3 reasonably easily. */
    
    s1 = sin(angle1) * s2/sin(angle2); 
    s3 = sin(angle3) * s2/sin(angle2);

    /* We now return the fractional distance along the edge. */

    *sA = s1/plc_distance(A[0],A[1]);
    *sB = 1 - (s3/plc_distance(B[0],B[1]));
    *ok = true;
    
}
   
bool tag_as_trouble(plc_vector A[2], plc_vector B[2]) {

  if (plc_distance(A[0],A[1]) < 1e-8 || plc_distance(B[0],B[1]) < 1e-8) { 

    return true;

  } 

  int ssAB0,ssAB1,ssBA0,ssBA1; // segment side (of segment) A (of point) B[0], and so forth.

  ssAB0 = segment_side(A[0],B[0]); ssAB1 = segment_side(A[0],B[1]); 
  ssBA0 = segment_side(B[0],A[0]); ssBA1 = segment_side(B[0],A[1]);

  if (ssAB0*ssAB1 == -1 && ssBA0*ssBA1 == -1) { /* Edges definitely DO cross */

    double sA, sB;
    bool ok;

    locate_crossing(A,B,&sA,&sB,&ok);

    if (ok) { return false; } 
    else {return true; }

  } else if ((ssAB0*ssAB1 == 0 && ssBA0*ssBA1 != +1) || /* An endpoint of B is colinear with A and B is not separated from A */
	     (ssAB0*ssAB1 != +1 && ssBA0*ssBA1 == 0)) { / * An endpoint of A is colinear with B and A is not separated from B */

      return true;

  } 

  return false;  /* No trouble here; we can either separate the segments or resolve the crossing. */
							   
}

/* We now make various passes through the curve in an attempt to make modifications
   which don't change topology. We'll fold in the octrope library here, because that's
   the best way to know what radius modifications are acceptable without changing knot type. */



/* Now we need to isolate the crossings and populate the crossing container. */

crossing_container *findcrossings(plCurve *L) {

  crossing_container *cc = crossing_container_new(1024);
  int cmp1,cmp2,vt1,vt2;

  /* Make sure all components are closed */

  for(cmp1=0;cmp1<L->nc;cmp1++) { 
    
    if(L->cp[cmp1].open == true) { 

      fprintf(stderr,"plcurve: Can't compute crossing code for open components.\n");
      exit(1);

    }
    
  }

  /* Main loop */

  for(cmp1=0;cmp1<L->nc;cmp1++) { 

    for(cmp2=cmp1;cmp2<L->nc;cmp2++) { 

      for(vt1=0;vt1<L->nc;vt1++) { 

	for(vt2=vt2+2;vt2<L->nc;vt2++) { 

	  /* Check for a crossing between 

	     edge L->cp[cmp1].vt[vt1] -- L->cp[cmp1].vt[vt1+1]
	     edge L->cp[cmp2].vt[vt2] -- L->cp[cmp2].vt[vt2+1]

	     If the crossing occurs at the _far_ end of either
	     edge, lookahead to try to resolve the crossing.
	     
	     If the crossing occurs at the _near_ end of either
	     edge, ignore it. */

}
   



char *plc_ccode( plCurve *L )
{
  char *code;

  if (L->nc != 1) {

    printf("plc_ccode: L has %d components, but this version only works for 1 comp.\n",
	   L->nc);
    return NULL;

  }

  /* Step 0. Change to a new temporary working directory to avoid stomping
     any user files */

  char *template,*newdir;

  template = calloc(32,sizeof(char));
  sprintf(template,"codedirXXXXXX");

  newdir = mkdtemp(template);
  assert(newdir != NULL);

  int result;
  result = chdir(newdir);
  assert(result != -1);

  /* Step 1. Write a res.pts file in compEric (21.17f) floating point format. */

  int vt;
  FILE *resfile;

  resfile = fopen("res.pts","w");
  assert(resfile != NULL);


  for(vt=0;vt<L->cp[0].nv;vt++) {

    fprintf(resfile,"  %21.17f  %21.17f  %21.17f\n",
	    plc_M_clist(L->cp[0].vt[vt]));

  }

  fclose(resfile);

  /* Step 2. Open a pipe to plcompEric and pass in answers
     to the various questions to make it go. */

  int edges;
  FILE *CEpipe;

  edges = plc_num_edges(L);

  CEpipe = popen("plcompEric > /dev/null","w");
  assert(CEpipe != NULL);

  /* plcompEric expects

 	# of edges
 	# of reps  (needs to be at least 1)
 	use restart? 1 yes, 0 no (we use 1 because we want to use res.pts)
 	save points in zpoints? 1/0 (no, so 0)
 	random seed, integer (useless, we set to 12345)
 	ropelength (useless here, we set to 30, but don't know why)
 	radius of perturbation (set to 0 to use this config)

	Note: comp516 makes the zmatrix file (also zkdata and zpoints, don't need) */

  fprintf(CEpipe,"%d\n1\n1\n1\n12345\n30\n0.0\n",edges);
  pclose(CEpipe);

  /* Step 4. Allocate space for the crossing code. */

  int csize = 256*32;
  code = calloc( csize, sizeof(char));

  /* Step 3. Copy the zmatrix file generated by compEric. */

  FILE *zmatrix;
  int   spos = 0;

  zmatrix = fopen("zmatrix","r");
  assert(zmatrix != NULL);

  int zmchar;
  while((zmchar = fgetc(zmatrix)) != EOF) {

    code[spos++] = (char)(zmchar);

    if (spos > csize-10) {

      csize *= 2;
      code = realloc(code,csize*sizeof(char));
      assert(code != NULL);

    }

  }

  fclose(zmatrix);

  code[spos] = 0; /* Add the terminating 0 manually */

  /* Step 4. Check and return. */

  //printf("Processed %d edge knot successfully. Crossing code is \n %s\n",
  //	 edges,code);


  /* Step 5. Cleanup the working directory. */

  remove("res.pts");
  remove("zkdata");
  remove("zmatrix");
  remove("zpoints");

  result = chdir("..");
  assert(result != -1);

  remove(newdir);
  free(template);

  return code;

}

char *plc_homfly( plCurve *L )
/* Compute homfly polynomial by calling the hidden plc_lmpoly function */
/* By default, this version times out after 60 seconds. */

{
  char *ccode;
  char *homfly;

  ccode = plc_ccode(L);
  homfly = plc_lmpoly(ccode,60);
  free(ccode);

  if (homfly == NULL) {
    homfly = calloc(128,sizeof(char));
    sprintf(homfly,"ccode unable to create crossing code for L");
  }

  return homfly;
}

int homcmp(const void *A, const void *B)
{
  plc_knottype *a,*b;

  a = (plc_knottype *)(A);
  b = (plc_knottype *)(B);

  return strcmp(a->homfly,b->homfly);
}

/* Find the knot type of a plCurve */
/* Sets nposs to the number of possible knottypes found for the curve. If we cannot
   classify the knot, return 0 for nposs and NULL for the buffer of knot types. */
plc_knottype *plc_classify( plCurve *L, int *nposs)

{
  plc_knottype kt = { 1, { 0 }, { 1 }, { "(1,1,e)" }, { "[[1]]N " } },*ktmatch,*ret;
  char *homfly;

  plc_knottype unknot = {1, { 0 }, { 1 }, { "(1,1,e)" }, { "[[1]]N " } };

  /* We have to work around a bug in lmpoly here, where it horks on one or two
     crossing diagrams. */

  char *ccode,*cptr;
  ccode = plc_ccode(L);  /* The number of crossings is 2 + the number of \n's in ccode. */

  int Ncount = 0;
  for(cptr = strchr(ccode,'\n');cptr != NULL;cptr = strchr(cptr+1,'\n'),Ncount++);

  if (Ncount == 3 || Ncount == 4) { /* 1 or 2 crossings, we know this must be the unknot */

    ret = calloc(1,sizeof(plc_knottype));
    assert(ret != NULL);
    ret[0] = unknot;
    *nposs = 1;

    free(ccode);
    return ret;

  }


  /* Now we know that the ccode has at least 3 crossings, compute the homfly (lmpoly should work) */

  homfly = plc_lmpoly(ccode,60);
  free(ccode);

  if (homfly == NULL) {
      homfly = calloc(128,sizeof(char));
      sprintf(homfly,"ccode unable to create crossing code for L");
  }

  if (homfly == NULL) { *nposs = 0; return NULL; }
  strncpy(kt.homfly,homfly,MAXHOMFLY);
  free(homfly);

  /* Now search for matching knottypes in ktdb */

  ktmatch = bsearch(&kt,ktdb,KTDBSIZE,sizeof(plc_knottype),homcmp);

  if (ktmatch == NULL) { /* No matching homfly was found */

    *nposs = 0;
    return NULL;

  }

  /* If there is more than one match, we might have landed anywhere in the collection of matching types */

  plc_knottype *mlo,*mhi;

  for(mlo = ktmatch;!strcmp(kt.homfly,mlo->homfly) && mlo > ktdb;mlo--);  // search down until we don't match anymore or start of buffer
  if (strcmp(kt.homfly,mlo->homfly)) { mlo++; }                           // if we don't match (we might if we went to the start of buffer)
  assert(!strcmp(kt.homfly,mlo->homfly));

  for(mhi = ktmatch;!strcmp(kt.homfly,mhi->homfly) && mhi < &(ktdb[KTDBSIZE]);mhi++);
  if (strcmp(kt.homfly,mhi->homfly)) { mhi--; }
  assert(!strcmp(kt.homfly,mhi->homfly));

  /* Now we know the interval between mlo and mhi (inclusive) is the set of knot types which match this homfly */

  *nposs = (mhi - mlo) + 1;
  ret = calloc(*nposs,sizeof(plc_knottype));
  assert(ret != NULL);

  int i;
  for(ktmatch=mlo,i=0;ktmatch <= mhi;i++,ktmatch++) {ret[i] = *ktmatch;} // Copy the matches to the output buffer
  return ret;

}
