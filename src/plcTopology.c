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
#include"plCurve.h"
#include"plcTopology.h"
#include"homfly.h"
#include"octrope.h"


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

#include"polynomials.h"

/*
   Convert to a Millett/Ewing representation and then compute HOMFLY
   using the venerable lmpoly code.

   The Millett/Ewing representation of a knot diagram numbers the
   crossings from 1 to ncrossings and then stores for each crossing
   the crossing connected to each arc coming from the crossing in the
   order

       a
       ^
       |
       |
   d<--|-->b (meaning that this can go either way, depending on the orientation of cr)
       |
       |
       c

   Note that the crossings are ordered CLOCKWISE, and that the +/- convention is
   the same as our usual one (see the Millett/Ewing paper, goddamn it!).

       d                    a
       ^                    ^
       |                    |
  c--------->a        d-----|----->b
       |                    |
       |                    |
       b                    c
   + crossing            - crossing

   So a crossing code representation of a plCurve is a char buffer
   containing lines of the form

   17+2b10c11c31a

   meaning that crossing 17 is a positive crossing

   connected in the a position to the b position of crossing 2,
   connected in the b position to the c position of crossing 10,
   connected in the c position to the c position of crossing 11 and
   connected in the d position to the a position of crossing 31.

   In order to simplify communication with the lmpoly code of Ewing
   and Millett, we store the crossing code as a standard (0
   terminated) string, including newlines. We will read from
   that string using a replacement version of the "read" primitive.

   We first need a primitive which tells us which of the four positions
   a given edge comes into or out of a crossing at: */

void pd_millet_ewing_positions(pd_code_t *pd, pd_idx_t edge, char *headpos, char *tailpos)

{
  pd_check_edge(SRCLOC,pd,edge);
  pd_check_notnull(SRCLOC,"headpos",headpos);
  pd_check_notnull(SRCLOC,"tailpos",tailpos);

  /* First, we get the head position */

  pd_idx_t oip,oop; /* overstrand-in-pos, overstrand-out-pos */
  pd_idx_t uip,uop; /* understrand-in-pos, understrand-out-pos */

  pd_overstrand_pos(pd,pd->edge[edge].head,&oip,&oop);
  pd_understrand_pos(pd,pd->edge[edge].head,&uip,&uop);

  if (pd->cross[pd->edge[edge].head].edge[oip] == edge) {

    /* We found it coming IN as an overstrand -- this is position c */
    *headpos = 'c';

  } else if (pd->cross[pd->edge[edge].head].edge[uip] == edge) {

    /* We found it coming IN as an understrand. This is either b or d,
       depending on the sign of the crossing. */

    if (pd->cross[pd->edge[edge].head].sign == PD_POS_ORIENTATION) {

      *headpos = 'b';

    } else {

      *headpos = 'd';

    }

  } else {

    pd_error(SRCLOC,
	     "in pd \n"
	     "%PD \n"
	     "%EDGE has head %d (%CROSS), \n"
	     "but doesn't appear in the overstrand in position (%d) \n"
	     "or understrand in position (%d).",
	     pd, edge, pd->edge[edge].head, pd->edge[edge].head,
	     oip,uip);

  }

  /* Now we deal with the tail position */

  pd_overstrand_pos(pd,pd->edge[edge].tail,&oip,&oop);
  pd_understrand_pos(pd,pd->edge[edge].tail,&uip,&uop);

  if (pd->cross[pd->edge[edge].tail].edge[oop] == edge) {

    /* We found it going OUT as an overstrand -- this is position a */
    *tailpos = 'a';

  } else if (pd->cross[pd->edge[edge].tail].edge[uop] == edge) {

    /* We found it going OUT as an understrand. This is either b or d,
       depending on the sign of the crossing. */

    if (pd->cross[pd->edge[edge].tail].sign == PD_POS_ORIENTATION) {

      *tailpos = 'd';

    } else {

      *tailpos = 'b';

    }

  } else {

    pd_error(SRCLOC,
	     "in pd \n"
	     "%PD \n"
	     "%EDGE has tail %d (%CROSS), \n"
	     "but doesn't appear in the overstrand out position (%d) \n"
	     "or understrand out position (%d).",
	     pd, edge, pd->edge[edge].tail, pd->edge[edge].tail,
	     oop,uop);

  }

}

/* We can now write the translator: */

char *pdcode_to_ccode(pd_code_t *pd)

{
  char *ccode;
  int   ccode_size = 256*pd->ncross;

  ccode = calloc(ccode_size,sizeof(char)); /* Get room for a big string */

  pd_idx_t cr;
  char    *ccode_buf = ccode;
  int      thiscr_used,total_used=0;

  for(cr=0;cr<pd->ncross;cr++) {

    char or = pd->cross[cr].sign == PD_POS_ORIENTATION ? '+':'-';

    thiscr_used = snprintf(ccode_buf,ccode_size-total_used,"%d%c",cr+1,or);
    total_used += thiscr_used;
    ccode_buf += thiscr_used;

    pd_idx_t a_edge, b_edge, c_edge, d_edge;
    bool b_is_in; /* Keep track of whether the b position is coming in or going out */
    pd_idx_t a_to, b_to, c_to, d_to;
    char a_to_pos, b_to_pos, c_to_pos, d_to_pos;
    char scratch;

    pd_overstrand(pd,cr,&c_edge,&a_edge);

    a_to = pd->edge[a_edge].head + 1;  /* The +1 is because M/E codes are 1-based */
    pd_millet_ewing_positions(pd,a_edge,&a_to_pos,&scratch);

    c_to = pd->edge[c_edge].tail + 1;  /* The +1 is because M/E codes are 1-based */
    pd_millet_ewing_positions(pd,c_edge,&scratch,&c_to_pos);

    b_is_in = (pd->cross[cr].sign == PD_POS_ORIENTATION);

    if (b_is_in) {

      pd_understrand(pd,cr,&b_edge,&d_edge);

      b_to = pd->edge[b_edge].tail + 1;
      pd_millet_ewing_positions(pd,b_edge,&scratch,&b_to_pos);

      d_to = pd->edge[d_edge].head + 1;
      pd_millet_ewing_positions(pd,d_edge,&d_to_pos,&scratch);

    } else { /* b is going OUT */

      pd_understrand(pd,cr,&d_edge,&b_edge);

      b_to = pd->edge[b_edge].head + 1;
      pd_millet_ewing_positions(pd,b_edge,&b_to_pos,&scratch);

      d_to = pd->edge[d_edge].tail + 1;
      pd_millet_ewing_positions(pd,d_edge,&scratch,&d_to_pos);

    }

    /* Having gathered all the information, we can write it to
       the buffer now. */

    thiscr_used = snprintf(ccode_buf,ccode_size-total_used,"%d%c%d%c%d%c%d%c\n",
			   a_to,a_to_pos,
			   b_to,b_to_pos,
			   c_to,c_to_pos,
			   d_to,d_to_pos
			   );

    total_used += thiscr_used;
    ccode_buf += thiscr_used;

    if (total_used >= ccode_size) { /* This shouldn't result in a write to unsafe memory because
				       we used snprintf, but is a sign that something's gone
				       horribly wrong. */

      pd_error(SRCLOC,"conversion to millet-ewing code failed for %PD at Millet-Ewing code\n %s",pd,ccode_buf);

    }

  }

  thiscr_used = snprintf(ccode_buf,ccode_size-total_used,"\n\n");

  total_used += thiscr_used;
  ccode_buf += thiscr_used;

  return ccode;

}


char *plc_lmpoly(char *ccode,int timeout); // This is in pllmpoly02.c.
char *pd_homfly_timeout(pd_code_t *pd, int timeout)
{
  char *ccode;
  char *homfly_lmpoly;

  ccode = pdcode_to_ccode(pd);

  if (PD_VERBOSE > 50) {

    pd_printf("pd_homfly: Converted pd code\n %PD \n to crossing code \n %s \n",
	     pd,ccode);

  }

  homfly_lmpoly = plc_lmpoly(ccode,timeout);

  if (PD_VERBOSE > 50) {

    printf("got (raw) HOMFLY string %s\n",homfly_lmpoly);

  }

  if (homfly_lmpoly == NULL) {
      return NULL; // NULL signifies timeout or other error in plc_lmpoly
  }

  /* Now transform to (more readable) latex form */

  char *homfly = lmpoly_to_latex(homfly_lmpoly);
  free(homfly_lmpoly);

  if (PD_VERBOSE > 50) {

    printf("got finished HOMFLY string %s\n",homfly);

  }

  free(ccode);
  ccode = NULL;

  return homfly;

}
char *pd_homfly(pd_code_t *pd) {
  
  pd_code_t *pds = pd_simplify(pd);
  char *pd_homfly_result = pd_homfly_timeout(pd, 300); /* Maximum of 5 minutes */
  pd_code_free(&pds);
  
  return pd_homfly_result;
}

char *plc_homfly( gsl_rng *rng, plCurve *L)

{
  pd_code_t *pd = pd_code_from_plCurve(rng,L);
  if (pd == NULL) { return NULL; } // Prevent segfault
  char *homfly = pd_homfly(pd);
  pd_code_free(&pd);
  return homfly;
}

/* Convert plCurve to pd_code_t, projecting along a direction and
   resolving geometric degeneracies as needed. */

/* Plan:

   Find crossings of each component with itself, then crossings
   of pairs of components, add to a container of crossings. */

typedef struct crossing_struct {

  // The crossing struct encodes parameter values along each
  // component, as well as component numbers. This is the object
  // stored in the master list of crossings.

  int lower_cmp;
  int upper_cmp;
  double lower_s;
  double upper_s;
  int sign;

} crossing;

typedef struct crossing_reference {

  // This crossing struct stores an index into the
  // the main array of crossings. There are actually
  // two crossing_reference entries for each crossing,
  // as they are stored at the upper and lower position
  // and component.

  int crossnum;
  double s;

} crossing_reference;

typedef struct crossing_container_struct {

  // A simple array which resizes.

  int size;
  int used;
  crossing *buf;

} crossing_container;

typedef struct crossing_reference_container_struct {

  int size;
  int used;
  crossing_reference *buf;

} crossing_reference_container;

int crossing_compare(const void *A, const void *B) {

  crossing *a = (crossing *)(A);
  crossing *b = (crossing *)(B);

  int acmp,bcmp;
  double as = -1.0,bs = -1.0;

  if (a->lower_cmp < a->upper_cmp) {

    acmp = a->lower_cmp; as = a->lower_s;

  } else {

    acmp = a->upper_cmp; as = a->upper_s;

  }

  if (b->lower_cmp < b->upper_cmp) {

    bcmp = b->lower_cmp; bs = b->lower_s;

  } else {

    bcmp = b->upper_cmp; bs = b->upper_s;

  }

  if (acmp != bcmp) {

    return acmp - bcmp;

  } else {

    assert(as > 0 && bs > 0);
    if (as < bs) { return -1; }
    else { return +1; }

  }

}


int crossing_reference_compare(const void *A,const void *B) {

  crossing_reference *a = (crossing_reference *)(A);
  crossing_reference *b = (crossing_reference *)(B);

  if (a->s < b->s) { return -1; }
  else { return 1; }

}

crossing_container *crossing_container_new(int size)
{

  crossing_container *ret;

  ret = (crossing_container *)(calloc(1,sizeof(crossing_container)));
  if (ret == NULL) { exit(1); }

  ret->size = size;
  ret->used = 0;
  ret->buf = calloc(size,sizeof(crossing));

  if (ret->buf == NULL) { exit(1); }

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
  a0b = plc_vect_diff(B,A[0]);
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

  int i,j;

  for(i=0;i<2;i++) {
    for(j=0;j<2;j++) {
      if (plc_distance(A[i],B[j]) < 1e-8) { return true; }
    }
  }

  int ssAB0,ssAB1,ssBA0,ssBA1; // segment side (of segment) A (of point) B[0], and so forth.

  ssAB0 = segment_side(A,B[0]); ssAB1 = segment_side(A,B[1]);
  ssBA0 = segment_side(B,A[0]); ssBA1 = segment_side(B,A[1]);

  if (ssAB0*ssAB1 == -1 && ssBA0*ssBA1 == -1) { /* Edges definitely DO cross */

    double sA, sB;
    bool ok;

    locate_crossing(A,B,&sA,&sB,&ok);

    if (ok) { return false; }
    else {return true; }

  } else if ((ssAB0*ssAB1 == 0 && ssBA0*ssBA1 != +1) ||
    /* An endpoint of B is colinear with A and B is not separated from A */
	     (ssAB0*ssAB1 != +1 && ssBA0*ssBA1 == 0)) {
    /* An endpoint of A is colinear with B and A is not separated from B */

      return true;

  }

  return false;  /* No trouble here; we can either separate the segments or resolve the crossing. */

}

/* We now make various passes through the curve in an attempt to make
   modifications which don't change topology. This returns a
   z-axis-clean version of plCurve with the same topology, or NULL if
   it can't find one (probably because the curve has a
   self-intersection in space). */

plCurve *make_zprojection_generic(gsl_rng *rng, plCurve *L) {

  plCurve *Lcopy,*Lprojection;
  Lcopy = plc_copy(L);

  int  cmpA,cmpB,vtA,vtB;
  bool projection_clean = false;
  int  attempt_counter;

  for(attempt_counter = 0; attempt_counter < 5 && !projection_clean; attempt_counter++) {

    /* We are going to run through the edge pairs in L, attempting to
       perturb our way out of any trouble that we run into. To do so,
       we'll need to know how far we dare perturb any vertex.

       This data is provided by octrope. We cap the perturbation at 1e-2,
       because a very large perturbation could run us into trouble
       in some unexpected way. */

    double safe_radius = (0.3)*octrope_thickness(Lcopy,NULL,0,0);
    safe_radius = (safe_radius > 1e-2) ? 1e-2 : safe_radius;

    if (safe_radius < 1e-16) {

      fprintf(stderr,"plc_ccode: plCurve is numerically singular. Can't compute projection.\n");
      return NULL;

    }

    Lprojection = plc_copy(Lcopy);
    plc_project(Lprojection,plc_build_vect(0,0,1));

    projection_clean = true;

    for(cmpA=0;cmpA < Lcopy->nc;cmpA++) {

      for(cmpB=cmpA;cmpB < Lcopy->nc;cmpB++) {

	for(vtA=0;vtA < Lcopy->cp[cmpA].nv;vtA++) {

	  int firstBvert, lastBvert;

	  /* If we're on the same component, then we start the B count _after_ A
	     and end before the last vertex to avoid comparing the same pair of
	     edges twice. We also have to be careful never to check adjacent edges,
	     which will always have a degenerate crossing where they meet. Remember
	     that there's a last edge which goes from nv-1 to nv (which is really
	     the same as vertex 0 because of wraparound).

	     If we're on different components, we need to run the A and B counts
	     across the entire component. */

	  if (cmpA == cmpB) {
	    firstBvert = vtA+2;
	    if (vtA==0) {
	      lastBvert = L->cp[cmpB].nv-1; // Do second-to-last, but not bridge edge.
	    } else {
	      lastBvert = L->cp[cmpB].nv; // Compare with bridge edge.
	    }
	  }
	  else { firstBvert = 0; lastBvert = L->cp[cmpB].nv; }

	  for(vtB=firstBvert;vtB < lastBvert;vtB++) {

	    if (tag_as_trouble(&(Lprojection->cp[cmpA].vt[vtA]),
			       &(Lprojection->cp[cmpB].vt[vtB]))) {

	      projection_clean = false;

	      plc_perturb_vect(rng,&Lcopy->cp[cmpA].vt[vtA],safe_radius);
	      plc_perturb_vect(rng,&Lcopy->cp[cmpA].vt[vtA+1],safe_radius);

	      plc_perturb_vect(rng,&Lcopy->cp[cmpB].vt[vtB],safe_radius);
	      plc_perturb_vect(rng,&Lcopy->cp[cmpB].vt[vtB+1],safe_radius);

	      plc_fix_wrap(Lcopy);

	    }

	  }

	}

      }

    }

    /* Now that we've compared all pairs of vertices (agonizingly
       slowly, but what other option do we really have?) we should
       have checked for and eliminated all the nongeneric crossings in
       this projection. We sweep again to make sure. */

    plc_free(Lprojection);

  }

  if (projection_clean == false) {

    fprintf(stderr,"plc_ccode: Couldn't resolve singularities in %d passes. plCurve is probably singular.\n",attempt_counter);
    return NULL;

  }

  return Lcopy;

}

/* Now we need to isolate the crossings and populate the crossing container. */

crossing_container *findcrossings(plCurve *L) {

  crossing_container *cc = crossing_container_new(1024);
  int cmpA,cmpB,vtA,vtB;

  /* Make sure all components are closed */

  for(cmpA=0;cmpA<L->nc;cmpA++) {

    if(L->cp[cmpA].open == true) {

      fprintf(stderr,"plc_ccode Can't compute crossing code for open components.\n");
      exit(1);

    }

  }

  plCurve *Lproj;
  Lproj = plc_copy(L);
  plc_project(Lproj,plc_build_vect(0,0,1));

  /* Main loop */

  for(cmpA=0;cmpA < L->nc;cmpA++) {

    for(cmpB=cmpA;cmpB < L->nc;cmpB++) {

      for(vtA=0;vtA < L->cp[cmpA].nv;vtA++) {

	int firstBvert, lastBvert;

	/* If we're on the same component, then we start the B count _after_ A
	   and end before the last vertex to avoid comparing the same pair of
	   edges twice. We also have to be careful never to check adjacent edges,
	   which will always have a degenerate crossing where they meet. Remember
	   that there's a last edge which goes from nv-1 to nv (which is really
	   the same as vertex 0 because of wraparound).

	   If we're on different components, we need to run the A and B counts
	   across the entire component. */

	if (cmpA == cmpB) {
	  firstBvert = vtA+2;
	  if (vtA==0) {
	    lastBvert = L->cp[cmpB].nv-1; // Do second-to-last, but not bridge edge.
	  } else {
	    lastBvert = L->cp[cmpB].nv; // Compare with bridge edge.
	  }
	}
	else { firstBvert = 0; lastBvert = L->cp[cmpB].nv; }

	for(vtB=firstBvert;vtB < lastBvert;vtB++) { // We have to avoid comparing the last edge with the first.

	  /* Check for a crossing between

	     edge L->cp[cmpA].vt[vtA] -- L->cp[cmpA].vt[vtA+1]
	     edge L->cp[cmpB].vt[vtB] -- L->cp[cmpB].vt[vtB+1]

	     Since we've already made the projection generic,
	     this shouldn't occur at either end.
	  */

	  int ssAB0,ssAB1,ssBA0,ssBA1; // segment side (of segment) A (of point) B[0], and so forth.

	  plc_vector *A,*B;
	  A = &(Lproj->cp[cmpA].vt[vtA]); B = &(Lproj->cp[cmpB].vt[vtB]);

	  ssAB0 = segment_side(A,B[0]); ssAB1 = segment_side(A,B[1]);
	  ssBA0 = segment_side(B,A[0]); ssBA1 = segment_side(B,A[1]);

	  if (ssAB0*ssAB1 == -1 && ssBA0*ssBA1 == -1) { /* Edges definitely DO cross */

	    double sA, sB;
	    bool ok;

	    locate_crossing(A,B,&sA,&sB,&ok);

	    if (ok) {

	      double htA, htB;
	      htA = (1-sA)*L->cp[cmpA].vt[vtA].c[2] + (sA)*L->cp[cmpA].vt[vtA+1].c[2];
	      htB = (1-sB)*L->cp[cmpB].vt[vtB].c[2] + (sB)*L->cp[cmpB].vt[vtB+1].c[2];
	      crossing cross;

	      plc_vector Atan, Btan;
	      Atan = plc_vect_diff(Lproj->cp[cmpA].vt[vtA+1],Lproj->cp[cmpA].vt[vtA]);
	      Btan = plc_vect_diff(Lproj->cp[cmpB].vt[vtB+1],Lproj->cp[cmpB].vt[vtB]);
	      plc_vector vcross;

	      if (htA > htB) {

		cross.upper_cmp = cmpA; cross.upper_s = vtA + sA;
		cross.lower_cmp = cmpB; cross.lower_s = vtB + sB;

		vcross = plc_cross_prod(Atan,Btan);
		cross.sign = vcross.c[2] > 0 ? +1 : -1;

	      } else {

		cross.lower_cmp = cmpA; cross.lower_s = vtA + sA;
		cross.upper_cmp = cmpB; cross.upper_s = vtB + sB;

		vcross = plc_cross_prod(Btan,Atan);
		cross.sign = vcross.c[2] > 0 ? +1 : -1;

	      }

	      crossing_container_add(cc,cross);

	    } else {

	      fprintf(stderr,"plc_ccode: Unexpected error computing "
		      "crossings for (supposedly) generic L.\n");
	      crossing_container_free(&cc);
	      return NULL;

	    }

	  }

	}

      }

    }

  }

  /* Now we've filled the crossing container. We can free the projection curve
     and contine to assembling the pd_code from this data. */

  plc_free(Lproj);

  qsort(cc->buf,cc->used,sizeof(crossing),crossing_compare);
  return cc;

}

crossing_reference_container *divide_crossings_by_component(crossing_container *cc,plCurve *L) {

  crossing_reference_container *crc = calloc(L->nc,sizeof(crossing_reference_container));
  int cr,cmp;

  for(cr=0;cr<cc->used;cr++) {

    crc[cc->buf[cr].lower_cmp].size++; crc[cc->buf[cr].upper_cmp].size++;

  }

  for(cmp=0;cmp<L->nc;cmp++) {

    if (crc[cmp].size == 0) {

      if (L->nc == 1) {

	// We have a single unknotted component with no crossings. Add a "virtual" self-crossing
	// in order to generate a valid pd_code. Remember that this crossing will be
	// added twice to the component record, so we set size to 2, not 1.

	crc[cmp].size = 2;

	crossing cross;
	cross.lower_cmp = cmp;
	cross.upper_cmp = cmp;
	cross.lower_s = 0.33;
	cross.upper_s = 0.66;
	cross.sign = +1;

	crossing_container_add(cc,cross);

      } else {

	printf("pd_code_from_plCurve: This build does not handle split components of links.\n");
	exit(1);

      }

    }

    crc[cmp].buf = calloc(crc[cmp].size,sizeof(crossing_reference));
    crc[cmp].used = 0;

  }

  for(cr=0;cr<cc->used;cr++) {

    cmp = cc->buf[cr].lower_cmp;
    crc[cmp].buf[crc[cmp].used].crossnum = cr;
    crc[cmp].buf[crc[cmp].used].s = cc->buf[cr].lower_s;
    crc[cmp].used++;

    cmp = cc->buf[cr].upper_cmp;
    crc[cmp].buf[crc[cmp].used].crossnum = cr;
    crc[cmp].buf[crc[cmp].used].s = cc->buf[cr].upper_s;
    crc[cmp].used++;

  }

  for(cmp=0;cmp<L->nc;cmp++) {

    if (crc[cmp].used != crc[cmp].size) { /* There's a problem! */

      fprintf(stderr,"pd_code_from_plCurve: Expected component %d to have %d crossings, but actually had %d.\n",
	      cmp,crc[cmp].size,crc[cmp].used);
      exit(1);

    }

    qsort(crc[cmp].buf,crc[cmp].size,sizeof(crossing_reference),crossing_reference_compare);

  }

  return crc;

}

enum height_t { lower, upper };

typedef struct crossing_strand_struct {

  pd_idx_t crossingnumber;
  pd_idx_t cmp;
  double   s;
  enum height_t height;

} crossing_strand_t;

int crossing_strand_cmp(const void *A,const void *B)
{
  crossing_strand_t *a = (crossing_strand_t *)(A);
  crossing_strand_t *b = (crossing_strand_t *)(B);

  if (a->cmp != b->cmp) {

    return b->cmp-a->cmp;

  }

  if (fabs(a->s - b->s) < 1e-8) { return 0; } /* These are the same guy! */

  return (a->s < b->s) ? -1 : 1;

}

pd_code_t *assemble_pdcode(plCurve *L,crossing_reference_container crc[],
			   crossing_container *cc) {

  /* We now use the component-wise crossing references and the
     crossing container's crossing numbering scheme to assemble a
     pd_code_t by walking around every component filling in the
     details.

     Let's review what we have:

     crc[i]->buf identifies all the crossings encountered along
     component i, in arclength order cc->buf identifies all the
     crossings by lower and upper component number and arclength and
     gives a sign for each crossing.

     We want to condense from this a list of oriented edges and
     crossings. We can make some elementary observations:

     1) Each entry in crc[i]->buf defines an edge uniquely. We can say
        that the head of the edge is at the arclength position given
        in crc[i]->buf[j].s, and at the crossing given by
        crc[i]->buf[s].crossnum.

     2) We can use these s values to identify the corresponding edges
        by their entries in the the cc->buf.

     So the first part of our process is going to be to replace the
     crc and cc by edge records records, telling us where the head and
     tail of each edge are, and whether they are upper or lower at
     that point.

     This is going to require some searching and sorting. The first
     thing we can do is to break up the crossing records in the
     crossing container into easier-to-sort-and-search pieces.

  */

  crossing_strand_t *cs = calloc(2*cc->used,sizeof(crossing_strand_t));
  pd_idx_t i;

  for(i=0;i<cc->used;i++) {

    cs[2*i].crossingnumber = i;
    cs[2*i].cmp = cc->buf[i].lower_cmp;
    cs[2*i].s = cc->buf[i].lower_s;
    cs[2*i].height = lower;

    cs[2*i+1].crossingnumber = i;
    cs[2*i+1].cmp = cc->buf[i].upper_cmp;
    cs[2*i+1].s = cc->buf[i].upper_s;
    cs[2*i+1].height = upper;

  }

  /* void qsort (void *array, size_t count, size_t size,
                 comparison_fn_t compare) */
  qsort(cs,2*cc->used,sizeof(crossing_strand_t),crossing_strand_cmp);

  /* Now we're going to look for evidence of a double crossing... */

  for(i=1;i<2*cc->used;i++) {

    if (crossing_strand_cmp(&(cs[i]),&(cs[i-1])) == 0) { /* Same component and s! */

      /* We've found a double crossing. Quit gracefully-- we'll
	 need another projection to make this work. */

      free(cs);
      return NULL;

    }

  }

  /* Now that we've made the list of crossings searchable, we're going
     to assemble edge data. Again, we'll first put our buffer of edges
     into an "expanded" form, then translate that into the information
     we need. */

  typedef struct pd_temp_edge_struct {

    pd_idx_t head;
    pd_idx_t tail;

    enum height_t head_height;
    enum height_t tail_height;

    pd_idx_t cmp;

  } pd_temp_edge_t;

  pd_idx_t total_edges = 0;
  pd_idx_t cmp,edge;

  for(cmp=0;cmp<L->nc;cmp++) {

    total_edges += crc[cmp].size;

  }

  pd_temp_edge_t *tempedges = calloc(total_edges,sizeof(pd_temp_edge_t));

  for(edge=0,cmp=0;cmp<L->nc;cmp++) {

    for(i=0;i<crc[cmp].size;i++,edge++) {

      /* We're now entering a record for the edge with HEAD at this
	 crossing. We know the component number and s for this edge,
	 but not the head crossing or over/under data.  However, we
	 can get this by matching the data we have to the cs array,
	 using

	 void * bsearch (const void *key, const void *array,
	 size_t count, size_t size, comparison_fn_t compare)

      */

      assert(edge < total_edges); /* Make sure we don't overrun the buffer */

      crossing_strand_t key;
      key.cmp = cmp;
      crossing_strand_t *csrec;

      /* First, we search for the head data (which is easier, because
	 this is the head) */

      key.s = crc[cmp].buf[i].s;
      csrec = bsearch(&key,cs,2*cc->used,
		      sizeof(crossing_strand_t),crossing_strand_cmp);
      assert(csrec != NULL);
      assert(csrec - cs < 2*cc->used);

      tempedges[edge].head = csrec->crossingnumber;
      tempedges[edge].head_height = csrec->height;

      /* Next, we search for the tail data (this is harder because we
	 need to wrap around if i == 0). Note: The tail of the edge is
	 part of the same component, so we don't need to reset
	 key.cmp. */

      key.s =
	(i == 0) ? crc[cmp].buf[crc[cmp].size-1].s : crc[cmp].buf[i-1].s;
      csrec =
	bsearch(&key,cs,2*cc->used,sizeof(crossing_strand_t),
		crossing_strand_cmp);
      assert(csrec != NULL);
      assert(csrec - cs < 2*cc->used);

      tempedges[edge].tail = csrec->crossingnumber;
      tempedges[edge].tail_height = csrec->height;

      /* Now we fill in the component: */

      tempedges[edge].cmp = cmp;

    }

  }

  free(cs); cs = NULL;

  /* We now need to assign positions for the edges at each crossing,
     and come up with crossing records which are compatible with the
     given headpos and tailpos positions and with the signs of the
     crossings.  We're going to use something like the KnotTheory
     conventions to do this in a consistent way: note that we can
     determine position uniquely from the combination of crossing sign,
     strand height (lower/upper) and head/tail.

            3                 3
            ^		      ^
            |		      |
      4 ----------> 2   4<-----------2
            |		      |
            |                 |
            1		      1

         positive         negative
         crossing         crossing

  */

  pd_code_t *pd = pd_code_new(cc->used > 1 ? cc->used : 2);
  /* Clear as much space as we need */
  assert(pd != NULL);

  pd->ncross = cc->used;
  pd->nedges = 2*pd->ncross;
  assert(pd->nedges == total_edges);
  assert(pd->nedges <= pd->MAXEDGES);
  assert(pd->ncross <= pd->MAXVERTS);

  /* We start by setting the crossing signs */

  for(i=0;i<pd->ncross;i++) {

    pd->cross[i].sign
      = (cc->buf[i].sign > 0) ? PD_POS_ORIENTATION : PD_NEG_ORIENTATION;

  }

  /* Now we can go ahead and assign positions. */

  for(i=0;i<pd->nedges;i++) {

    pd->edge[i].head = tempedges[i].head;
    pd->edge[i].tail = tempedges[i].tail;

    pd_idx_t head = pd->edge[i].head;
    pd_idx_t tail = pd->edge[i].tail;

    /* Recall our rules:

            2                 2
            ^		      ^
            |		      |
      3 ----------> 1  3 <----------- 1
            |		      |
            |                 |
            0		      0

         positive         negative
         crossing         crossing
    */

    if (tempedges[i].head_height == lower) {

      pd->edge[i].headpos = 0;
      pd->cross[head].edge[0] = i;

    } else { /* The head_height is upper */

      if (pd->cross[head].sign == PD_POS_ORIENTATION) {

	pd->edge[i].headpos = 3;
	pd->cross[head].edge[3] = i;

      } else {

	pd->edge[i].headpos = 1;
	pd->cross[head].edge[1] = i;

      }

    }

     /* Recall our rules:

            2                 2
            ^		      ^
            |		      |
      3 ----------> 1  3 <----------- 1
            |		      |
            |                 |
            0		      0

         positive         negative
         crossing         crossing
    */

    if (tempedges[i].tail_height == lower) {

      pd->edge[i].tailpos = 2;
      pd->cross[tail].edge[2] = i;

    } else { /* The tail_height is upper */

      if (pd->cross[tail].sign == PD_POS_ORIENTATION) {

	pd->edge[i].tailpos = 1;
	pd->cross[tail].edge[1] = i;

      } else {

	pd->edge[i].tailpos = 3;
	pd->cross[tail].edge[3] = i;

      }

    }

  }

  /* At this point, we've shed the tempedges, too. */
  free(tempedges);
  tempedges = NULL;

  pd_regenerate(pd);

  /* Now we know how many components we should have found from the
     plCurve data. */

  assert(pd->ncomps == L->nc);

  /* This is the most important test: */

  if (!pd_ok(pd)) {

      fprintf(stderr,"pd_code_from_plCurve: Algorithm produced "
	      "the invalid pd_code\n");
      pd_write(stderr,pd);
      fprintf(stderr,"\n"
	      "from plCurve input. This is a bug in the library,"
	      "or a truly singular example polygon.\n");
      exit(1);

  }

  return pd;

}

/* We can now finally assemble everything. We use a source of randomness to defeat nongeneric configurations. */

/* There are two "secret" functions used for debugging */

bool pd_code_from_plCurve_verbose = false;
bool pd_code_from_plCurve_debug = false;

void set_pd_code_from_plCurve_verbose(bool val) {

  pd_code_from_plCurve_verbose = val;

}

void set_pd_code_from_plCurve_debug(bool val) {

  pd_code_from_plCurve_debug = val;

}

pd_code_t *pd_code_from_plCurve(gsl_rng *rng, plCurve *L) {

  int        attempt_number;
  bool       pd_created = false;
  plCurve    *workingL = plc_copy(L);
  pd_code_t  *pd = NULL;
  plc_vector new_axis = {{0,0,1}};

  for(attempt_number=0;attempt_number < 9 && !pd_created;attempt_number++) {

    /* Step 1. Rotate the whole thing randomly. */

    if (!pd_code_from_plCurve_debug) { new_axis = plc_random_vect(); }

    plc_random_rotate(workingL,new_axis);

    if (pd_code_from_plCurve_verbose) {

      printf("pd_code_from_plCurve\n"
	     "------------------------------------------------\n"
	     "Starting run on %d vertex, %d component plCurve.\n",plc_num_verts(L),L->nc);
      printf("Rotating so that (%g,%g,%g) becomes the z-axis...done\n",plc_M_clist(new_axis));
      printf("Perturbing to make crossing-clean...");
    }

    /* Step 2. Perturb to make it crossing-clean. */

    plCurve *genericL = make_zprojection_generic(rng,L);

    if (pd_code_from_plCurve_verbose) {

      if (genericL != NULL) {

	printf("success.\n");

      } else {

	printf("failure on attempt %d.\n",attempt_number);

      }

    }

    if (genericL == NULL) { /* We didn't succeed. Randomly rotate, try again. */

      pd_created = false;
      continue;

    }

    /* Step 3. Populate crossing container. */

    if (pd_code_from_plCurve_verbose) {

      printf("Searching for crossings in generic version of plCurve...");

    }

    crossing_container *cc = findcrossings(genericL);

    if (pd_code_from_plCurve_verbose) {

      if (cc != NULL) {

	printf("success "
	       "(%d crossings)\n",cc->used);

      } else {

	printf("failure on attempt %d.\n",attempt_number);

      }

    }

    if (cc == NULL) {

      pd_created = false;
      plc_free(genericL);
      continue;

    }

    /* Step 4. Sort the crossings by component. */

    if (pd_code_from_plCurve_verbose) {

      printf("Splitting crossings up by component...");

    }

    crossing_reference_container *crc = divide_crossings_by_component(cc,genericL);

    if (pd_code_from_plCurve_verbose) {

      printf("done.\n"
	     "\t%d crossings are divided among %d components as:\n",cc->used,L->nc);

      int i,total=0;
      for(i=0;i<L->nc;i++) {

	printf("\t%d crossings on component %d (%d vertices)\n",crc[i].used,i,L->cp[i].nv);
	total += crc[i].used;

      }

      printf("total number of crossings on components %d ",total);

      if (total != 2*cc->used) {

	printf("!= expected number %d = 2 * %d crossings total\n",
	       2*cc->used,cc->used);
	exit(1);

      } else {

	printf("== expected number %d = 2 * %d crossings total\n",
	       2*cc->used,cc->used);

      }

    }

    /* Step 5. Assemble and test the pd_code */

    if (pd_code_from_plCurve_verbose) {

      printf("assembling pd code from crossing and crossing reference data...");

    }

    pd = assemble_pdcode(L,crc,cc);

    if (pd_code_from_plCurve_verbose) {

      if (pd != NULL) {

	printf("success\n");
	printf("generated %d crossing pd code from crossing container of %d crossings\n"
	       "---------------------------------------------------------------------\n\n",
	       pd->ncross,cc->used);
      } else {

	printf("failure on attempt %d.\n",attempt_number);

      }

    }

    if (pd == NULL) {

      pd_created = false;

    } else {

      pd_created = true;

    }

    /* Step 6. Housekeeping! */

    int cmp;
    for(cmp=0;cmp<L->nc;cmp++) {

      if (crc[cmp].buf != NULL) { free(crc[cmp].buf); crc[cmp].buf = NULL; }

    }
    free(crc);
    crossing_container_free(&cc);
    plc_free(genericL);

  }

  plc_free(workingL);

  if (pd_created) {

    return pd;

  } else {

    fprintf(stderr,"pd_code_from_plCurve: Couldn't create pd_code from %d vertex, %d component plCurve after %d attempts.\n",
	    plc_num_verts(L),L->nc,attempt_number);
    return NULL;

  }

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
plc_knottype *plc_classify(gsl_rng *rng, plCurve *L, int *nposs)

{
  plc_knottype kt = { 1, { 0 }, { 1 }, { "(1,1,e)" }, { "[[1]]N " } },*ktmatch,*ret;
  char *homfly;

  plc_knottype unknot = {1, { 0 }, { 1 }, { "(1,1,e)" }, { "[[1]]N " } };

  /* We have to work around a bug in lmpoly here, where it horks on one or two
     crossing diagrams. */

  char *ccode,*cptr;
  pd_code_t *pdC;

  pdC = pd_code_from_plCurve(rng,L);
  if (pdC == NULL) { *nposs = 0; return NULL; } /* For a singular polygon, don't classify. */

  ccode = pdcode_to_ccode(pdC);  /* The number of crossings is 2 + the number of \n's in ccode. */
  free(pdC); /* We're not going to use this again, may as well free it now */

  if (ccode == NULL) { *nposs = 0; return NULL; } // Otherwise you WILL segfault

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

plc_knottype *pd_classify(pd_code_t *pdC, int *nposs)

{
  plc_knottype kt = { 1, { 0 }, { 1 }, { "(1,1,e)" }, { "[[1]]N " } },*ktmatch,*ret;
  char *homfly;

  plc_knottype unknot = {1, { 0 }, { 1 }, { "(1,1,e)" }, { "[[1]]N " } };

  /* We have to work around a bug in lmpoly here, where it horks on one or two
     crossing diagrams. */

  char *ccode,*cptr;

  if (pdC == NULL) { *nposs = 0; return NULL; } /* For a singular polygon, don't classify. */

  if (pdC->ncross == 0 || pdC->ncross == 1 || (pdC->ncomps == 1 && pdC->ncross == 2)) {

    ret = calloc(1,sizeof(plc_knottype));
    assert(ret != NULL);
    ret[0] = unknot;
    *nposs = 1;
    return ret;

  }

  ccode = pdcode_to_ccode(pdC);  /* The number of crossings is 2 + the number of \n's in ccode. */
  free(pdC); /* We're not going to use this again, may as well free it now */

  if (ccode == NULL) { *nposs = 0; return NULL; } // Otherwise you WILL segfault

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
void plc_write_knottype(FILE *out, plc_knottype kt)
/* Prints the knot type (or types) in a neatly formatted human-readable version. */
{
  int i;
  fprintf(out,"%d_%d",kt.cr[0],kt.ind[0]);
  for(i=1;i<kt.nf;i++) {
    fprintf(out,"#%d_%d",kt.cr[i],kt.ind[i]);
  }
  fprintf(out," (%s) \n",kt.homfly);
}


/*** Over and under information **/

/* All of this will depend on an internal function. This will be called more than once in typical use,
   but it's so fast that the time penalty is negligible in practice and more than offset by the coding
   clarity of having the under and over information returned separately! */

 /* The convention used to determine sign is this:

            ^
            |
       ----------->
            |
            |

      positive crossing
      (upper tangent vector) x (lower tangent vector) points OUT of screen.

            ^
            |
       -----|----->
            |
            |

      negative crossing
      (upper tangent vector) x (lower tangent vector) points INTO screen.

 */

void pd_overunder_internal(pd_code_t *pd, pd_idx_t cr,
			   pd_idx_t *over_in_pos, pd_idx_t *over_out_pos,
			   pd_idx_t *under_in_pos, pd_idx_t *under_out_pos)

{
    pd_check_cr(SRCLOC,pd,cr);
    pd_check_notnull(SRCLOC,"over_in_pos",over_in_pos); pd_check_notnull(SRCLOC,"over_out_pos",over_out_pos);
    pd_check_notnull(SRCLOC,"under_in_pos",under_in_pos); pd_check_notnull(SRCLOC,"under_out_pos",under_out_pos);

    if (pd->edge[pd->cross[cr].edge[0]].head == cr &&
        pd->edge[pd->cross[cr].edge[0]].headpos == 0) { /* The edge at position 0 is incoming */

        if (pd->edge[pd->cross[cr].edge[1]].head == cr &&
            pd->edge[pd->cross[cr].edge[1]].headpos == 1) { /* The edge at position 1 is incoming */

            if (pd->cross[cr].sign == PD_POS_ORIENTATION) {

                // We're in the case:
                //
                //        3
                //        ^
                //        |
                //  0 -----------> 2
                //        |
                //        |
                //        1
                //
                //

                *over_in_pos = 0; *over_out_pos = 2;
                *under_in_pos = 1; *under_out_pos = 3;

            } else if (pd->cross[cr].sign == PD_NEG_ORIENTATION) { /* The sign of the crossing is negative */

                // We're in the case:
                //
                //        3
                //        ^
                //        |
                //  0 ----|------> 2
                //        |
                //        |
                //        1
                //
                //

                *over_in_pos = 1; *over_out_pos = 3;
                *under_in_pos = 0; *under_out_pos = 2;

            } else if (pd->cross[cr].sign == PD_UNSET_ORIENTATION) { /* The sign of the crossing is unset */

                // We're in the case:
                //
                //       3
                //       ^
                //       |
                // 0 ----+------> 2
                //       |
                //       |
                //       1
                //
                // For an unset crossing... the over strand obeys "positive" sign rules
                //

                *over_in_pos = 0; *over_out_pos = 2;
                *under_in_pos = 1; *under_out_pos = 3;

            }

        } else { /* The edge at position 1 is outgoing */

            if (pd->cross[cr].sign == PD_POS_ORIENTATION) {

                // We're in the case:
                //
                //        3
                //        |
                //        |
                //  0 ----|------> 2
                //        |
                //        |
                //        v
                //        1
                //
                //

                *over_in_pos = 3; *over_out_pos = 1;
                *under_in_pos = 0; *under_out_pos = 2;

            } else if (pd->cross[cr].sign == PD_NEG_ORIENTATION) { /* The sign of the crossing is negative */
                // We're in the case:
                //
                //       3
                //       |
                // 0 ---------> 2
                //       |
                //       v
                //       1

                *over_in_pos = 0; *over_out_pos = 2;
                *under_in_pos = 3; *under_out_pos = 1;

            } else if (pd->cross[cr].sign == PD_UNSET_ORIENTATION) {
                // We're in the case:
                //
                //      3
                //      |
                // 0 ---+---> 2
                //      |
                //      v
                //      1

                *over_in_pos = 3; *over_out_pos = 1;
                *under_in_pos = 0; *under_out_pos = 2;

            }

        }

    } else { /* the edge at position 0 is outgoing! */

        if (pd->edge[pd->cross[cr].edge[1]].head == cr &&
            pd->edge[pd->cross[cr].edge[1]].headpos == 1) { /* The edge at position 1 is incoming */

            if (pd->cross[cr].sign == PD_POS_ORIENTATION) {
                // We're in the case:
                //
                //       3
                //       ^
                //       |
                // 0 <---|---- 2
                //       |
                //       1

                *over_in_pos = 1; *over_out_pos = 3;
                *under_in_pos = 2; *under_out_pos = 0;

            } else if (pd->cross[cr].sign == PD_NEG_ORIENTATION) { /* The sign of the crossing is negative */
                // We're in the case:
                //
                //        3
                //        ^
                //        |
                //  0 <------- 2
                //        |
                //        1

                *over_in_pos = 2; *over_out_pos = 0;
                *under_in_pos = 1; *under_out_pos = 3;

            } else if (pd->cross[cr].sign == PD_UNSET_ORIENTATION) {
                // We're in the case:
                //
                //       3
                //       ^
                //       |
                // 0 <---+---- 2
                //       |
                //       1

                *over_in_pos = 1; *over_out_pos = 3;
                *under_in_pos = 2; *under_out_pos = 0;

            }

        } else { /* The edge at position 1 is outgoing */

            if (pd->cross[cr].sign == PD_POS_ORIENTATION) {
                // We're in the case:
                //
                //       3
                //       |
                // 0 <-------- 2
                //       |
                //       v
                //       1

                *over_in_pos = 2; *over_out_pos = 0;
                *under_in_pos = 3; *under_out_pos = 1;

            } else if (pd->cross[cr].sign == PD_NEG_ORIENTATION) { /* The sign of the crossing is negative */
                // We're in the case:
                //
                //       3
                //       |
                // 0 <---|--- 2
                //       |
                //       v
                //       1

                *over_in_pos = 3; *over_out_pos = 1;
                *under_in_pos = 2; *under_out_pos = 0;

            } else if (pd->cross[cr].sign == PD_UNSET_ORIENTATION) {
                // We're in the case:
                //
                //       3
                //       |
                // 0 <---+---- 2
                //       |
                //       v
                //       1

                *over_in_pos = 2; *over_out_pos = 0;
                *under_in_pos = 3; *under_out_pos = 1;

            }

        }

    }

}

void pd_overstrand_pos(pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgepos, pd_idx_t *outgoing_edgepos)

/* Returns the position in crossing cr of pd (that is, a number in 0..3) of the incoming and outgoing
   edges of the strand going over at this crossing, using the crossing sign data
   and edge orientations in order to compute the answer. */
{
  pd_idx_t oip,oop;
  pd_idx_t uip,uop;

  pd_overunder_internal(pd,cr,&oip,&oop,&uip,&uop);
  *incoming_edgepos = oip;
  *outgoing_edgepos = oop;
}

void pd_understrand_pos(pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgepos, pd_idx_t *outgoing_edgepos)

/* Returns the position in crossing cr of pd (that is, a number in 0..3) of the incoming and outgoing
   edges of the strand going under at this crossing, using the crossing sign data
   and edge orientations in order to compute the answer. */

{
  pd_idx_t oip,oop;
  pd_idx_t uip,uop;

  pd_overunder_internal(pd,cr,&oip,&oop,&uip,&uop);
  *incoming_edgepos = uip;
  *outgoing_edgepos = uop;
}

void pd_overstrand(pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgenum, pd_idx_t *outgoing_edgenum)

/*  Returns the edge number  (that is, a number in 0..pd->nedges)
    of the incoming and outgoing edges of the strand going OVER at crossing cr of pd,
    using the sign of the crossing to determine. */

{
  pd_idx_t oip, oop;

  pd_overstrand_pos(pd,cr,&oip,&oop);
  *incoming_edgenum = pd->cross[cr].edge[oip];
  *outgoing_edgenum = pd->cross[cr].edge[oop];
}

void pd_understrand(pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgenum, pd_idx_t *outgoing_edgenum)

/*  Returns the edge number  (that is, a number in 0..pd->nedges)
    of the incoming and outgoing edges of the strand going UNDER at crossing cr of pd,
    using the sign of the crossing to determine. */
{
  pd_idx_t uip, uop;

  pd_understrand_pos(pd,cr,&uip,&uop);
  *incoming_edgenum = pd->cross[cr].edge[uip];
  *outgoing_edgenum = pd->cross[cr].edge[uop];
}


int pd_linking_number(pd_code_t *pd,pd_idx_t c1,pd_idx_t c2)

/* Computes the linking number of two components of the pdcode. */
/* Requires that crossing signs be set; otherwise, fails out. */

{
  pd_check_cmp(SRCLOC,pd,c1);
  pd_check_cmp(SRCLOC,pd,c2);

  /* First, we check that all the crossing orientations are set. */

  pd_idx_t i;

  for(i=0;i<pd->ncross;i++) {

    assert(pd->cross[i].sign != PD_UNSET_ORIENTATION);

  }

  /* We'll loop over the shorter component. */

  pd_component_t *shortCmp;
  pd_idx_t longCmpNum;

  if (pd->comp[c1].nedges < pd->comp[c2].nedges) {

    shortCmp = &(pd->comp[c1]);
    longCmpNum = c2;

  } else {

    shortCmp = &(pd->comp[c2]);
    longCmpNum = c1;

  }

  /* Now we loop over shortCmp, identifying the other
     component at each crossing and seeing if it's
     longCmpNum. */

  int lk = 0;

  pd_idx_t cmpEdge,otherEdgeNum;
  pd_edge_t *e;
  pd_crossing_t *cr;

  for(cmpEdge=0;cmpEdge<shortCmp->nedges;cmpEdge++) {

    e = &(pd->edge[shortCmp->edge[cmpEdge]]);
    cr = &(pd->cross[e->head]);
    otherEdgeNum = cr->edge[(e->headpos+1)%4]; /* Get the intersecting edge */

    pd_idx_t otherComp,otherPos;
    pd_component_and_pos(pd,otherEdgeNum,&otherComp,&otherPos);

    if (otherComp == longCmpNum) { /* This crossing counts! */

      if(cr->sign == PD_UNSET_ORIENTATION) {

	lk++;

      } else if (cr->sign == PD_POS_ORIENTATION) {

	lk++;

      } else {

	lk--;

      }

    }

  }

  lk /= 2;
  return lk;

}

unsigned int pd_unsigned_linking_number(pd_code_t *pd,pd_idx_t c1,pd_idx_t c2)

/* Computes the linking number of two components of the pdcode. */
/* Requires that crossing signs be set; otherwise, fails out. */

{
  pd_check_cmp(SRCLOC,pd,c1);
  pd_check_cmp(SRCLOC,pd,c2);

  /* We'll loop over the shorter component. */

  pd_component_t *shortCmp;
  pd_idx_t longCmpNum;

  if (pd->comp[c1].nedges < pd->comp[c2].nedges) {

    shortCmp = &(pd->comp[c1]);
    longCmpNum = c2;

  } else {

    shortCmp = &(pd->comp[c2]);
    longCmpNum = c1;

  }

  /* Now we loop over shortCmp, identifying the other
     component at each crossing and seeing if it's
     longCmpNum. */

  int lk = 0;

  pd_idx_t cmpEdge,otherEdgeNum;
  pd_edge_t *e;
  pd_crossing_t *cr;

  for(cmpEdge=0;cmpEdge<shortCmp->nedges;cmpEdge++) {

    e = &(pd->edge[shortCmp->edge[cmpEdge]]);
    cr = &(pd->cross[e->head]);
    otherEdgeNum = cr->edge[(e->headpos+1)%4]; /* Get the intersecting edge */

    pd_idx_t otherComp,otherPos;
    pd_component_and_pos(pd,otherEdgeNum,&otherComp,&otherPos);

    if (otherComp == longCmpNum) { /* This crossing counts! */
	lk++;
    }

  }

  lk /= 2;
  return lk;

}
