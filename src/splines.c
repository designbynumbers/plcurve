/*
 * Routines to create, destroy, and convert spline equivalents of plCurves
 *
 * $Id: splines.c,v 1.31 2007-09-21 21:06:50 cantarel Exp $
 *
 * This code generates refinements of plCurves, component by component, using
 * the Numerical Recipes spline code for interpolation.
 *
 */

/* Copyright 2004 The University of Georgia. */

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

#include <math.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>

static inline void spline_strand_new(/*@out@*/ plc_spline_strand *Pl,
                                     int  ns, bool  open, int  cc) {

  assert(ns >= 1);

  Pl->open = open;
  Pl->ns = ns;

  Pl->svals = (double *)malloc((ns+2)*sizeof(double));
  assert(Pl->svals != NULL);
  Pl->svals++; /* so that Pl->svals[-1] is valid. */

  Pl->vt = (plc_vector *)malloc((ns+2)*sizeof(plc_vector));
  assert(Pl->vt != NULL);
  Pl->vt++; /* so that Pl->vt[-1] is a valid space */

  Pl->vt2 = (plc_vector *)malloc((ns+2)*sizeof(plc_vector));
  assert(Pl->vt2 != NULL);
  Pl->vt2++; /* so that Pl->vt2[-1] is a valid space */

  Pl->cc = cc;
  Pl->clr = (plc_color *)malloc(cc*sizeof(plc_color));
  assert(Pl->clr != NULL);
}

/*
 * Procedure allocates memory for a new spline. The number of
 * components is given by components. The number of data samples in each
 * component shows up in the buffer pointed to by ns.  The closed/open
 * nature of each strand is given in the array pointed to by open.
 *
 */
plc_spline *plc_spline_new(const int          components,
                                   const int  * const ns,
                                   const bool * const open,
                                   const int  * const cc)
{
  plc_spline *L;
  int i;

  /* First, we check to see that the input values are reasonable. */

  assert(components >= 1);
  assert(ns != NULL);
  assert(open != NULL);
  assert(cc != NULL);

  /* Now we attempt to allocate space for these components. */

  L = malloc(sizeof(plc_spline));
  assert(L != NULL);
  L->nc = components;
  L->cp = malloc(L->nc*sizeof(plc_spline_strand));
  assert(L->cp != NULL);

  for (i = 0; i < L->nc; i++) {
    spline_strand_new(&L->cp[i],ns[i],open[i],cc[i]);
  }

  L->cst = NULL;

  /*@-compdef@*/ /* Splint isn't sure that everything gets defined */
  return L;
  /*@=compdef@*/
} /* plc_spline_new */

/*
 * Free the memory used to hold vertices in a given strand (not the memory of
 * the strand itself).  We then set all the values in the strand data structure
 * to reflect the fact that the memory has been freed.  We can call
 * spline_strand_free twice on the same strand pointer without fear.
 *
 */
static inline void spline_strand_free(/*@null@*/ plc_spline_strand *Pl) {

  if (Pl == NULL) {
    return;
  }

  Pl->ns = 0;

  if (Pl->svals != NULL) {
    Pl->svals--; /* undo the original svals++ */
    free(Pl->svals);
    Pl->svals = NULL;
  }

  if (Pl->vt != NULL) {
    Pl->vt--; /* undo our original vt++ (for wraparound) */
    free(Pl->vt);
    Pl->vt = NULL;
  }

  if (Pl->vt2 != NULL) {
    Pl->vt2--; /* undo our original vt++ (for wraparound) */
    free(Pl->vt2);
    Pl->vt2 = NULL;
  }

  Pl->cc = 0;
  if (Pl->clr != NULL) {
    free(Pl->clr);
    Pl->clr = NULL;
  }

  /*@-nullstate@*/ /* We expect to return null fields, it's ok :-) */
  return;
  /*@=nullstate@*/
} /* strand_free */

/*
 * Free the memory associated with a given spline.  We then set all the values
 * in the spline data structure to reflect the fact that the memory has been
 * freed.  We can call plc_spline_free twice on the same spline pointer
 * without fear.
 *
 */
void plc_spline_free(/*@only@*/ /*@null@*/ plc_spline *L) {
  int i;

  /* First, we check the input. */
  if (L == NULL) {
    return; /* Move along, nothing to see here */
  }

  /* Now we can get to work. */
  /*@+loopexec@*/
  for (i=0; i<L->nc; i++) {
    spline_strand_free(&L->cp[i]);
  }
  /*@=loopexec@*/

  /*@-compdestroy@*/
  free(L->cp);
  /*@=compdestroy@*/
  L->cp = NULL;
  L->nc = 0;

  /* Now we destroy the linked list of constraints. */

  plc_constraint *cst;

  for(;L->cst != NULL;) {

    cst = L->cst;
    L->cst = L->cst->next;
    free(cst);

  }

  free(L);
  L = NULL;
} /* plcurve_spline_free */


/*
 * Converts a regular plCurve to spline form. The spline code is adapted from
 * the "Numerical recipes" spline code. Each line from NR is commented next to
 * the (hopefully) equivalent vectorized version below.
 */
plc_spline *plc_convert_to_spline(plCurve * const L,bool *ok)
{
  int   i;
  int  *ns;
  int  *cc;
  bool *open;
  int   comp;

  int         k, I;
  plc_vector p, un, *u;
  double      sig, qn;
  plc_vector yp1,ypn;
  double      scrx;
  plc_vector scrV,scrW;

  plc_spline *spline_L;

  /* First, we check for sanity */

  assert(L != NULL);
  assert(L->nc >= 1);

  plc_fix_wrap(L);

  ns = (int*)calloc((size_t)L->nc,sizeof(int));
  assert(ns != NULL);
  open = (bool *)calloc((size_t)L->nc,sizeof(bool));
  assert(open != NULL);
  cc = (int *)calloc((size_t)L->nc,sizeof(int));
  assert(cc != NULL);

  /*@+loopexec@*/
  for(i=0;i<L->nc;i++) {
    ns[i] = L->cp[i].nv;
    cc[i] = L->cp[i].cc;
    open[i] = L->cp[i].open;
  }
  /*@=loopexec@*/

  /* Now we allocate the new spline. */

  spline_L = plc_spline_new(L->nc,ns,open,cc);

  /* Done with this space now, no matter what happened */
  free(ns); free(open); free(cc);
  ns = NULL; open = NULL; cc = NULL;

  assert(spline_L != NULL);

  /* We now copy the constraint records into the new spline */

  plc_constraint *cst,*newcst;

  for(cst = L->cst;cst != NULL;cst = cst->next) {

    /* Make a copy of the constraint */

    newcst = (plc_constraint *)(calloc(1,sizeof(plc_constraint)));
    *newcst = *cst;

    /* Now insert it at the head of the list of constraints for the spline. */

    newcst->next = spline_L->cst;
    spline_L->cst = newcst;

  }

  /* We now go component-by-component. */

  for (comp = 0;comp < L->nc;comp++) {

    /* Our first task is to assemble svals: the arclength of each vertex. 
     * On an n-vertex open strand, the first point should be at 0 and the nth
     * at the length of the strand.  On an n-vertex closed strand, the first
     * point should be at 0 and the n+1st at the length of the strand. */

    spline_L->cp[comp].svals[-1] = 
      -plc_M_distance(L->cp[comp].vt[0],L->cp[comp].vt[-1]);
    spline_L->cp[comp].svals[0] = 0.0;

    for (i = 1; i <= L->cp[comp].nv; i++) {
      spline_L->cp[comp].svals[i] = spline_L->cp[comp].svals[i-1] +
        plc_M_distance(L->cp[comp].vt[i-1], L->cp[comp].vt[i]);
    }

    /* We have now assembled the "s" values for the spline. */
    /* We go ahead and copy the corresponding positions into place. */

    memcpy(&(spline_L->cp[comp].vt[-1]),&(L->cp[comp].vt[-1]),
        (L->cp[comp].nv+2)*sizeof(plc_vector));

    /* Last, we go ahead and copy color values in. */

    memcpy(spline_L->cp[comp].clr,L->cp[comp].clr,
        (L->cp[comp].cc)*sizeof(plc_color));

    /* We are now prepared to compute vt2 (the second derivatives), which will
       determine the actual form of the splined polyline. */

    /* double p, qn, sig, un, *u */

    /* Definitions moved to the top of the routine */

    /* u = malloc(n,sizeof(double)); */

    u = malloc((L->cp[comp].nv)*sizeof(plc_vector));
    assert(u != NULL);

    /* ...together with values yp1 and ypn for the first derivative at
       the first and last samples... */

    yp1 = plc_mean_tangent(L,comp,0,ok);

    if (L->cp[comp].open) {

      ypn = plc_mean_tangent(L,comp,L->cp[comp].nv-1,ok);

    } else {

      ypn = yp1;

      /* The idea here is that if we're splining a closed curve, we
         want to be able to smooth the last corner as well. */

    }

    spline_L->cp[comp].vt2[-1] = plc_build_vect(0.0,0.0,0.0);

    /* y2[1] = -0.5; */

    spline_L->cp[comp].vt2[0] = plc_build_vect(-0.5,-0.5,-0.5);

    /* u[1] = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1); */

    scrV = plc_vect_diff(spline_L->cp[comp].vt[1],spline_L->cp[comp].vt[0]);
    scrx = spline_L->cp[comp].svals[1] - spline_L->cp[comp].svals[0];
    plc_M_vlincomb(u[0],3.0/(scrx*scrx),scrV,-3.0/scrx,yp1);

    /* for (i=2;i<=n-1;i++) { */

    I = L->cp[comp].open ? L->cp[comp].nv-1 : L->cp[comp].nv;

    /*@+loopexec@*/
    for (i=1; i < I; i++) {

      /* sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]); */

      sig = (spline_L->cp[comp].svals[i] - spline_L->cp[comp].svals[i-1])/
        (spline_L->cp[comp].svals[i+1] - spline_L->cp[comp].svals[i-1]);

      /* p = sig*y2[i-1] + 2.0; */

      plc_M_vlincomb(p,sig,spline_L->cp[comp].vt2[i-1],1.0,
          plc_build_vect(2.0,2.0,2.0));

      /* y2[i] = (sig-1.0)/p; */

      spline_L->cp[comp].vt2[i] = plc_component_div(
          plc_build_vect(sig - 1.0,sig - 1.0,sig - 1.0),p,NULL);

      /* u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]); */

      scrV = plc_vect_diff(spline_L->cp[comp].vt[i+1],
                            spline_L->cp[comp].vt[i]);
      scrW = plc_vect_diff(spline_L->cp[comp].vt[i],
                            spline_L->cp[comp].vt[i-1]);

      plc_M_vlincomb(u[i],
        1.0/(spline_L->cp[comp].svals[i+1] - spline_L->cp[comp].svals[i]),scrV,
       -1.0/(spline_L->cp[comp].svals[i] - spline_L->cp[comp].svals[i-1]),scrW);

      /* u[i]=(6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p; */

      scrx = spline_L->cp[comp].svals[i+1] - spline_L->cp[comp].svals[i-1];
      plc_M_scale_vect(6.0/scrx,u[i]);
      plc_M_sub_vect(u[i],plc_scale_vect(sig,u[i-1]));
      plc_M_component_div(u[i],p);

    }
    /*@=loopexec@*/

    /* qn = 0.5; */

    qn = 0.5;

    /* un = (3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1])); */

    scrx = spline_L->cp[comp].svals[I] - spline_L->cp[comp].svals[I-1];
    scrV = plc_vect_diff(spline_L->cp[comp].vt[I],spline_L->cp[comp].vt[I-1]);

    plc_M_vlincomb(un,3.0/scrx,ypn,-3.0/(scrx*scrx),scrV);

    /* y2[n] = (un - qn*u[n-1])/(qn*y2[n-1]+1.0); */

    scrV.c[0] = scrV.c[1] = scrV.c[2] = 1.0;
    plc_M_vlincomb(scrW,qn,spline_L->cp[comp].vt2[I-1],1.0,scrV);
    /*@-usedef@*/
    plc_M_vlincomb(scrV,1.0,un,-qn,u[I-1]);
    /*@=usedef@*/

    spline_L->cp[comp].vt2[I] = plc_component_div(scrV,scrW,NULL);

    /* for (k=n-1;k>=1;k--) { */

    for (k = I - 1; k >= 0; k--) {

      /* y2[k] = y2[k]*y2[k+1] + u[k]; */

      plc_M_component_mult(spline_L->cp[comp].vt2[k],
                            spline_L->cp[comp].vt2[k+1]);
      /*@-usedef@*/
      plc_M_add_vect(spline_L->cp[comp].vt2[k],u[k]);
      /*@=usedef@*/

    }

    /* free(u); */

    free(u);
    u = NULL;

  }

  return spline_L;
}

inline static plc_vector evaluate_poly(const plc_spline * const spL,
                                        const int comp, 
                                        const int shi, 
                                        const int slo, 
                                        const double s,
                                    /*@out@*/ double *a,
                                    /*@out@*/ double *b) {
  double h;
  plc_vector scrV, scrW;

  h = spL->cp[comp].svals[shi] - spL->cp[comp].svals[slo];
  assert(h > DBL_EPSILON); /* They'd best not be the same point */

  /* FInd the appropriate weights */
  *a = (spL->cp[comp].svals[shi] - s)/h;
  *b = (s - spL->cp[comp].svals[slo])/h;

  plc_M_vlincomb(scrV,*a,spL->cp[comp].vt[slo],*b,spL->cp[comp].vt[shi]);

  plc_M_vlincomb(scrW, 
      ((*a)*(*a)*(*a) - (*a))*(h*h)/6.0, spL->cp[comp].vt2[slo],
      ((*b)*(*b)*(*b) - (*b))*(h*h)/6.0, spL->cp[comp].vt2[shi]);

  return plc_vect_sum(scrV,scrW);
}

inline static plc_vector evaluate_poly_derivative(const plc_spline * const spL,
						  const int comp, 
						  const int shi, 
						  const int slo, 
						  const double s,
						  /*@out@*/ double *a,
						  /*@out@*/ double *b) {
  double h;
  plc_vector dscrV, dscrW;

  h = spL->cp[comp].svals[shi] - spL->cp[comp].svals[slo];
  assert(h > DBL_EPSILON); /* They'd best not be the same point */

  /* Find the appropriate weights */
  *a = (spL->cp[comp].svals[shi] - s)/h;
  *b = (s - spL->cp[comp].svals[slo])/h;

  /* Now get ready to differentiate */

  double da, db;

  da = -1.0/h;
  db = 1.0/h;

  /* We give the old (vector) function as a comment, followed by the s derivative */

  /* scrV function: plc_M_vlincomb(scrV,*a,spL->cp[comp].vt[slo],*b,spL->cp[comp].vt[shi]); */

  plc_M_vlincomb(dscrV,da,spL->cp[comp].vt[slo],db,spL->cp[comp].vt[shi]);

  /* scrW function: 
     
     plc_M_vlincomb(scrW, 
     ((*a)*(*a)*(*a) - (*a))*(h*h)/6.0, spL->cp[comp].vt2[slo],
     ((*b)*(*b)*(*b) - (*b))*(h*h)/6.0, spL->cp[comp].vt2[shi]);

  */
     
  plc_M_vlincomb(dscrW, 
      (3*(*a)*(*a)*(da) - (da))*(h*h)/6.0, spL->cp[comp].vt2[slo],
      (3*(*b)*(*b)*(db) - (db))*(h*h)/6.0, spL->cp[comp].vt2[shi]);

  return plc_vect_sum(dscrV,dscrW);

}

/*
 * Converts spline back to regular plCurve, but changes the number of verts
 * to those in nv. We require that the number of components in nv match those
 * in spL.
 *
 * A new plCurve is allocated.
 */
plCurve *plc_convert_from_spline(const plc_spline * const spL,
                                     const int * const nv)
{
  int *cc;
  bool *open;
  int comp, i;
  double s, s_max;
  int shi, slo;
  double b, a;
  bool do_colors;

  plCurve *L;

  /* First, we do some elementary sanity checking. */

  assert(spL != NULL);
  assert(spL->nc >= 1);

  /* Now we allocate the new plCurve. */

  open = (bool *)malloc(spL->nc*sizeof(bool));
  assert(open != NULL);
  cc = (int *)malloc(spL->nc*sizeof(int));
  assert(cc != NULL);

  /*@+loopexec@*/
  for (i=0; i<spL->nc; i++) {
    open[i] = spL->cp[i].open;
    cc[i] = (spL->cp[i].cc <= 1) ? spL->cp[i].cc : nv[i];
  }
  /*@=loopexec@*/

  L = plc_new(spL->nc,nv,open,cc);
  assert(L != NULL);

  free(open); free(cc);
  open = NULL; cc = NULL;

  /* We now copy the constraint records into the new polyline */

  plc_constraint *cst,*newcst;
  
  for(cst = spL->cst;cst != NULL;cst = cst->next) {

    /* Make a copy of the constraint */

    newcst = (plc_constraint *)(calloc(1,sizeof(plc_constraint)));
    *newcst = *cst;

    /* Now insert it at the head of the list of constraints for the spline. */

    newcst->next = L->cst;
    L->cst = newcst;

  }

  /* Now we can start filling in vertex positions, component by component. */

  for (comp=0; comp < L->nc; comp++) {

    do_colors = (L->cp[comp].cc > 1);
    
    if (L->cp[comp].cc == 1) {
      L->cp[comp].clr[0] = spL->cp[comp].clr[0];
    } 

    if (spL->cp[comp].open) {
      s_max = spL->cp[comp].svals[spL->cp[comp].ns-1] * 
        ((double)L->cp[comp].nv / (double)(L->cp[comp].nv - 1));
    } else {
      s_max = spL->cp[comp].svals[spL->cp[comp].ns];
    }

    for (i=0, slo=0, shi=1; i < L->cp[comp].nv; i++) {

      s = s_max*i/L->cp[comp].nv;

      /* Search to make sure that slo and shi bracket our current s.*/

      while (spL->cp[comp].svals[shi] <= s && 
          shi < ((L->cp[comp].open) ? spL->cp[comp].ns-1 : spL->cp[comp].ns)) {
        shi++;
        slo++;
      }
      assert(spL->cp[comp].svals[shi] >= s);
      assert(spL->cp[comp].svals[slo] <= s);

      /* Now we can evaluate the spline polynomial. */

      L->cp[comp].vt[i] = evaluate_poly(spL,comp,shi,slo,s,&a,&b);

      /* And interpolate the colors */
      if (do_colors) {
        L->cp[comp].clr[i] = plc_build_color(
            a*spL->cp[comp].clr[slo].r     + b*spL->cp[comp].clr[shi].r,
            a*spL->cp[comp].clr[slo].g     + b*spL->cp[comp].clr[shi].g,
            a*spL->cp[comp].clr[slo].b     + b*spL->cp[comp].clr[shi].b,
            a*spL->cp[comp].clr[slo].alpha + b*spL->cp[comp].clr[shi].alpha);
      }
    }

    /* We now search the constraints to figure out which ones apply to this component */

    for(cst=L->cst;cst != NULL;cst=cst->next) {

      if (cst->cmp == comp) { /* We need to translate this component */

	double cst_slo,cst_shi;
	  
	cst_slo = spL->cp[comp].svals[cst->vert];
	cst_shi = spL->cp[comp].svals[cst->vert+cst->num_verts];
	
	/* We first set the start vert  */
	
	if (cst->vert == spL->cp[comp].ns-1) { /* Last vert */

	  cst->vert = L->cp[comp].nv-1;
	  
	} else if (cst->vert == 0) { /* First vert */

	  cst->vert = 0;

	} else { /* Somewhere in between-- work it out from the s values */

	   cst->vert = floor((cst_slo/s_max)*L->cp[comp].nv);

	}

	/* Now we set the number of verts */

	if (cst->kind == fixed || cst->num_verts == 1) { /* One vert */

	  cst->num_verts = 1;

	} else if (cst->num_verts == spL->cp[comp].ns) { /* All verts */

	  cst->num_verts = L->cp[comp].nv;

	} else { /* Scale the length from the old resolution */

	  cst->num_verts = ceil(((cst_shi-cst_slo)/s_max)*L->cp[comp].nv);
	  
	}	
	
      }

    }

  }

  plc_fix_wrap(L);
  return L;
}


/* Procedure samples the spline at a particular s value, returning a
 * spatial position. */
plc_vector plc_sample_spline(const plc_spline * const spL,
                                  const int cmp,
                                  double s)
{
  int klo,khi,k;
  double cmpLen;
  double b, a;

  /* We begin with a bit of checking to make sure that cmp and s seem
     compatible with the given spL. */

  assert(cmp >= 0);
  assert(cmp < spL->nc);

  /* Now fix any wraparound, so that the s value given is in [0,cmpLen). */

  cmpLen = spL->cp[cmp].svals[spL->cp[cmp].ns];
  while (s < DBL_EPSILON) {    s += cmpLen; }
  while (s >= cmpLen) { s -= cmpLen; }

  /* We now search for the correct s value. */

  klo = 0;
  khi = spL->cp[cmp].ns;

  while (khi-klo > 1) {

    /*@-shiftimplementation@*/
    k = (khi+klo) >> 1;
    /*@=shiftimplementation@*/
    if (spL->cp[cmp].svals[k] >= s) khi=k;
    else klo = k;

  }

  return evaluate_poly(spL,cmp,khi,klo,s,&a,&b);
}

/* Procedure samples the spline at a particular s value, returning a
 * spatial position. */
plc_vector plc_spline_tangent(const plc_spline * const spL,
                                  const int cmp,
                                  double s)
{
  int klo,khi,k;
  double cmpLen;
  double b, a;

  /* We begin with a bit of checking to make sure that cmp and s seem
     compatible with the given spL. */

  assert(cmp >= 0);
  assert(cmp < spL->nc);

  /* Now fix any wraparound, so that the s value given is in [0,cmpLen). */

  cmpLen = spL->cp[cmp].svals[spL->cp[cmp].ns];
  while (s < DBL_EPSILON) {    s += cmpLen; }
  while (s >= cmpLen) { s -= cmpLen; }

  /* We now search for the correct s value. */

  klo = 0;
  khi = spL->cp[cmp].ns;

  while (khi-klo > 1) {

    /*@-shiftimplementation@*/
    k = (khi+klo) >> 1;
    /*@=shiftimplementation@*/
    if (spL->cp[cmp].svals[k] >= s) khi=k;
    else klo = k;

  }

  return evaluate_poly_derivative(spL,cmp,khi,klo,s,&a,&b);
}

