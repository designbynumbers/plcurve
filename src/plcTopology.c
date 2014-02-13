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
#include <plcTopology.h>
#include <homfly.h>
#include <octrope.h>


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

  int ssAB0,ssAB1,ssBA0,ssBA1; // segment side (of segment) A (of point) B[0], and so forth.

  ssAB0 = segment_side(A,B[0]); ssAB1 = segment_side(A,B[1]); 
  ssBA0 = segment_side(B,A[0]); ssBA1 = segment_side(B,A[1]);

  if (ssAB0*ssAB1 == -1 && ssBA0*ssBA1 == -1) { /* Edges definitely DO cross */

    double sA, sB;
    bool ok;

    locate_crossing(A,B,&sA,&sB,&ok);

    if (ok) { return false; } 
    else {return true; }

  } else if ((ssAB0*ssAB1 == 0 && ssBA0*ssBA1 != +1) || /* An endpoint of B is colinear with A and B is not separated from A */
	     (ssAB0*ssAB1 != +1 && ssBA0*ssBA1 == 0)) { /* An endpoint of A is colinear with B and A is not separated from B */

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

       This data is provided by octrope. We cap the perturbation at 1e-4,
       because a very large perturbation could run us into trouble
       in some unexpected way. */

    double safe_radius = (0.3)*octrope_thickness(Lcopy,NULL,0,0); 
    safe_radius = (safe_radius > 1e-4) ? 1e-4 : safe_radius;
    
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
	    
	    if (tag_as_trouble(&(Lprojection->cp[cmpA].vt[vtA]),&(Lprojection->cp[cmpB].vt[vtB]))) {

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

	      fprintf(stderr,"plc_ccode: Unexpected error computing crossings for (supposedly) generic L.\n");
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

      fprintf(stderr,"pdcode_from_plCurve: Expected component %d to have %d crossings, but actually had %d.\n",
	      cmp,crc[cmp].size,crc[cmp].used);
      exit(1);

    } 

    qsort(crc[cmp].buf,crc[cmp].size,sizeof(crossing_reference),crossing_reference_compare);

  }

  return crc;

}

pd_code_t *assemble_pdcode(plCurve *L,crossing_reference_container crc[], crossing_container *cc) {

  // We now use the component-wise crossing references and the crossing container's crossing
  // numbering scheme to assemble a pd_code_t by walking around every component filling in 
  // the details. 

  pd_code_t *pd = pd_code_new(cc->used); /* Clear as much space as we need */
  assert(pd != NULL);

  pd->ncross = cc->used;
  pd->nedges = 2*pd->ncross;
  pd->ncomps = L->nc;

  int cmp,cr,edge;

  edge = 0;

  for(cmp=0;cmp<L->nc;cmp++) { 

    for(cr=0;cr<crc[cmp].used;cr++,edge++) { 

      // The edge data is always going to be ordered starting at upper in. 
      // This isn't required the by the pd_code_t spec, but it's going to 
      // help keep us organized.

      int pd_crossing = crc[cmp].buf[cr].crossnum; // This is the index into cc->buf (the master crossing list)
      
      pd->cross[pd_crossing].sign = cc->buf[pd_crossing].sign > 0 ? PD_POS_ORIENTATION : PD_NEG_ORIENTATION;

      if (cmp == cc->buf[pd_crossing].upper_cmp && fabs(crc[cmp].buf[cr].s - cc->buf[pd_crossing].upper_s) < 1e-13) { // This is the upper arc

	pd->cross[pd_crossing].edge[0] = (cr == 0) ? (edge+crc[cmp].used-1) : edge-1;  // Incoming edge is the one BEFORE this one, must wrap if we just started cmp
	pd->cross[pd_crossing].edge[2] = edge; // This (outgoing) edge is always the current edge number.

      } else { // This is the lower arc, which means that we have to pay attention to sign.

	if (pd->cross[pd_crossing].sign == PD_POS_ORIENTATION) { // Upper out to lower out is ccw, meaning that upper in to lower in is ccw
  
	  pd->cross[pd_crossing].edge[1] = (cr == 0) ? (edge+crc[cmp].used-1) : edge-1;  // Incoming edge is the one BEFORE this one, must wrap if we just started cmp
	  pd->cross[pd_crossing].edge[3] = edge; // This (outgoing) edge is always the current edge number.
	       
	} else { // Upper out to lower out is cw, meaning that upper in to lower out is ccw

	  pd->cross[pd_crossing].edge[3] = (cr == 0) ? (edge+crc[cmp].used-1) : edge-1;  // Incoming edge is the one BEFORE this one, must wrap if we just started cmp
	  pd->cross[pd_crossing].edge[1] = edge; // This (outgoing) edge is always the current edge number.
	 
	}

      }
      
    }

  }

  /* We should have assembled a (skeletal) pd_code now, containing only crossing information. */

  pd_regenerate(pd);

  /* This is the most important test: */

  if (!pd_ok(pd)) { 

    fprintf(stderr,"pd_code_from_plCurve: Algorithm produced the invalid pd_code\n");
    pd_write(stderr,pd);
    fprintf(stderr,"\n"
	    "from plCurve input. This is a bug in the library, or a truly singular example polygon.\n");
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
  pd_code_t  *pd;
  plc_vector new_axis = {{0,0,1}};
  
  for(attempt_number=0;attempt_number < 3 && !pd_created;attempt_number++) { 
  
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
    
    
      
      
    
    
  
  
  
  
      

/* We now have a convenience function to convert the resulting pdcode into 
   a Millett-Ewing crossing code. */

/* The Millett/Ewing representation of a knot diagram numbers the
   crossings from 1 to ncrossings and then stores for each crossing
   the crossing connected to each arc coming from the crossing in the
   order
   
         a
         |
         |
     b---|-->d
         |
         V
         c

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

  */

char *ccode_from_pd_code( pd_code_t *pdC) {

  /* This is just a stub for now. */

  return NULL;

}   

char *old_plc_ccode( plCurve *L )
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




char *plc_homfly(gsl_rng *rng, plCurve *L )
/* Compute homfly polynomial by calling the hidden plc_lmpoly function */
/* By default, this version times out after 60 seconds. */

{
  pd_code_t *pdC;
  char *ccode;
  char *homfly;

  pdC = pd_code_from_plCurve(rng,L);
  ccode = ccode_from_pd_code(pdC);

  homfly = plc_lmpoly(ccode,60);
  free(ccode);
  free(pdC);

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
  ccode = ccode_from_pd_code(pdC);  /* The number of crossings is 2 + the number of \n's in ccode. */
  free(pdC); /* We're not going to use this again, may as well free it now */

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
