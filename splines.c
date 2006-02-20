/*
 * Routines to create, destroy, and convert spline equivalents of plCurves
 *
 * $Id: splines.c,v 1.17 2006-02-20 22:23:39 ashted Exp $
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

#include "plCurve.h"

#ifdef HAVE_MATH_H
  #include <math.h>
#endif
#ifdef HAVE_MALLOC_H
  #include <malloc.h>
#endif
#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
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

static inline void spline_strand_new(plCurve_spline_strand *Pl,
                                                     int  ns, 
                                                     int  open, 
                                                     int  cc) {
 
  assert(ns >= 1);

  Pl->open = open;
  Pl->ns = ns;

  Pl->svals = (double *)malloc((ns+2)*sizeof(double));
  assert(Pl->svals != NULL);
  Pl->svals++; /* so that Pl->svals[-1] is valid. */

  Pl->vt = (plcl_vector *)malloc((ns+2)*sizeof(plcl_vector));
  assert(Pl->vt != NULL);
  Pl->vt++; /* so that Pl->vt[-1] is a valid space */

  Pl->vt2 = (plcl_vector *)malloc((ns+2)*sizeof(plcl_vector));
  assert(Pl->vt2 != NULL);
  Pl->vt2++; /* so that Pl->vt2[-1] is a valid space */

  Pl->cc = cc;
  Pl->clr = (plCurve_color *)malloc(cc*sizeof(plCurve_color));
  assert(Pl->clr != NULL);
  Pl->clr2 = (plCurve_color *)malloc(cc*sizeof(plCurve_color));
  assert(Pl->clr2 != NULL);
}

/*
 * Procedure allocates memory for a new spline. The number of
 * components is given by components. The number of data samples in each
 * component shows up in the buffer pointed to by ns.  The closed/open
 * nature of each strand is given in the array pointed to by open.
 *
 */
plCurve_spline *plCurve_spline_new(const int components, 
                                   const int * const ns, 
                                   const int * const open, 
                                   const int * const cc) 
{
  plCurve_spline *L;
  int i;

  /* First, we check to see that the input values are reasonable. */

  assert(components >= 1);
  assert(ns != NULL);
  assert(open != NULL);
  assert(cc != NULL);

  /* Now we attempt to allocate space for these components. */
  
  L = malloc(sizeof(plCurve_spline));
  assert(L != NULL);
  L->nc = components;
  L->cp = malloc(L->nc*sizeof(plCurve_spline_strand));
  assert(L->cp != NULL);

  for (i = 0; i < L->nc; i++) {
    spline_strand_new(&L->cp[i],ns[i],open[i],cc[i]);
  }

  return L;
} /* plCurve_spline_new */

/*
 * Free the memory used to hold vertices in a given strand (not the memory of
 * the strand itself).  We then set all the values in the strand data structure
 * to reflect the fact that the memory has been freed.  We can call
 * spline_strand_free twice on the same strand pointer without fear. 
 *
 */ 
static inline void spline_strand_free(plCurve_spline_strand *Pl) {
  
  if (Pl == NULL) {
    return;
  }

  Pl->ns = 0;
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
  }
  if (Pl->clr2 != NULL) {
    free(Pl->clr2);
  }

} /* strand_free */

/*
 * Free the memory associated with a given spline.  We then set all the values
 * in the spline data structure to reflect the fact that the memory has been
 * freed.  We can call plCurve_spline_free twice on the same spline pointer
 * without fear. 
 *
 */ 
void plCurve_spline_free(plCurve_spline *L) {
  int i;

  plcl_error_num = 0;

  /* First, we check the input. */
  if (L == NULL) {
    return; /* Move along, nothing to see here */
  }

  /* Now we can get to work. */
  for (i=0; i<L->nc; i++) {
    spline_strand_free(&L->cp[i]);
  }

  free(L->cp);
  L->nc = 0;

  free(L);
  L = NULL;
} /* plcurve_spline_free */


/*
 * Converts a regular plCurve to spline form. The spline code is adapted from
 * the "Numerical recipes" spline code. Each line from NR is commented next to
 * the (hopefully) equivalent vectorized version below. 
 */
plCurve_spline *plCurve_convert_to_spline(plCurve * const L)
{
  int    i;
  int    *ns,*cc,*open;
  int    comp;

  plCurve_spline *spline_L;

  /* First, we check for sanity */

  plcl_error_num = 0;

  assert(L != NULL);
  assert(L->nc >= 1);

  plCurve_fix_wrap(L);

  ns = calloc(L->nc,sizeof(int));
  assert(ns != NULL);
  open = calloc(L->nc,sizeof(int));
  assert(open != NULL);
  cc = calloc(L->nc,sizeof(int));
  assert(cc != NULL);

  /*@+loopexec@*/
  for(i=0;i<L->nc;i++) {
    ns[i] = L->cp[i].nv;
    cc[i] = L->cp[i].cc;
    open[i] = L->cp[i].open;
  }
  /*@=loopexec@*/

  /* Now we allocate the new spline. */

  spline_L = plCurve_spline_new(L->nc,ns,open,cc);

  /* Done with this space now, no matter what happened */
  free(ns); free(open); free(cc);

  if (plcl_error_num != 0 || spline_L == NULL) {
    return NULL;
  }
          
  /* We now go component-by-component. */

  for (comp = 0;comp < L->nc;comp++) {

    /* Our first task is to assemble svals: the arclength of each vertex. */

    spline_L->cp[comp].svals[0] = 0.0;

    for (i=1;i<L->cp[comp].nv;i++) {

      spline_L->cp[comp].svals[i] = 
        spline_L->cp[comp].svals[i-1] + plcl_M_distance(L->cp[comp].vt[i-1],
                                            L->cp[comp].vt[i]);
    }
    
    if (!L->cp[comp].open) { 

      spline_L->cp[comp].svals[i] = 
        spline_L->cp[comp].svals[i-1] + plcl_M_distance(L->cp[comp].vt[i-1],
                                                      L->cp[comp].vt[i]);     

      spline_L->cp[comp].svals[-1] = -plcl_M_distance(L->cp[comp].vt[0],
                                                    L->cp[comp].vt[-1]);

    } else {

      spline_L->cp[comp].svals[i] = spline_L->cp[comp].svals[i-1];
      spline_L->cp[comp].svals[-1] = spline_L->cp[comp].svals[0];

    }
    
    /* We have now assembled the "s" values for the spline. */
    /* We go ahead and copy the corresponding positions into place. */

    for(i=-1;i<=L->cp[comp].nv;i++) {

      spline_L->cp[comp].vt[i] = L->cp[comp].vt[i];

    }

    /* Last, we go ahead and copy color values in. */

    for(i=0;i<L->cp[comp].cc;i++) {

      spline_L->cp[comp].clr[i] = L->cp[comp].clr[i];

    }
    
    /* We are now prepared to compute vt2 and clr2, which will
       determine the actual form of the splined polyline. */

    int i,k,I;

    /* double p, qn, sig, un, *u */

    plcl_vector p,un,*u;
    double sig,qn;

    /* u = malloc(n,sizeof(double)); */

    u = malloc((spline_L->cp[comp].ns+2)*sizeof(plcl_vector));
    assert(u != NULL);

    /* ...together with values yp1 and ypn for the first derivative at 
       the first and last samples... */

    plcl_vector yp1,ypn;

    yp1 = plCurve_tangent_vector(L,comp,0);

    if (L->cp[comp].open) { 

      ypn = plCurve_tangent_vector(L,comp,L->cp[comp].nv-1);

    } else {

      ypn = yp1;

      /* The idea here is that if we're splining a closed curve, we
         want to be able to smooth the last corner as well. */

    }

    /* y2[1] = -0.5; */

    spline_L->cp[comp].vt2[0].c[0] = spline_L->cp[comp].vt2[0].c[1] 
      = spline_L->cp[comp].vt2[0].c[2] = -0.5;

    /* u[1] = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1); */

    double scrx;
    plcl_vector scrV,scrW;

    scrV = plcl_vect_diff(spline_L->cp[comp].vt[1],spline_L->cp[comp].vt[0]);
    scrx = spline_L->cp[comp].svals[1] - spline_L->cp[comp].svals[0];
    plcl_M_vlincomb(u[0],3.0/pow(scrx,2),scrV,-3.0/scrx,yp1);
    
    /* for (i=2;i<=n-1;i++) { */

    I = L->cp[comp].open ? L->cp[comp].nv-2 : L->cp[comp].nv-1;

    for(i=1;i<=I;i++) {

      /* sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]); */

      sig = (spline_L->cp[comp].svals[i] - spline_L->cp[comp].svals[i-1])/
        (spline_L->cp[comp].svals[i+1] - spline_L->cp[comp].svals[i-1]);

      /* p = sig*y2[i-1] + 2.0; */

      scrV.c[0] = scrV.c[1] = scrV.c[2] = 2.0;
      plcl_M_vlincomb(p,sig,spline_L->cp[comp].vt2[i-1],1.0,scrV);

      /* y2[i] = (sig-1.0)/p; */

      scrV.c[0] = scrV.c[1] = scrV.c[2] = sig - 1.0;
      spline_L->cp[comp].vt2[i] = plcl_component_div(scrV,p);

      /* u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]); */

      scrV = plcl_vect_diff(spline_L->cp[comp].vt[i+1],
                            spline_L->cp[comp].vt[i]);
      scrW = plcl_vect_diff(spline_L->cp[comp].vt[i],
                            spline_L->cp[comp].vt[i-1]);

      plcl_M_vlincomb(u[i],
        1.0/(spline_L->cp[comp].svals[i+1] - spline_L->cp[comp].svals[i]),scrV,
       -1.0/(spline_L->cp[comp].svals[i] - spline_L->cp[comp].svals[i-1]),scrW);

      /* u[i]=(6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p; */

      scrx = spline_L->cp[comp].svals[i+1] - spline_L->cp[comp].svals[i-1];
      plcl_M_scale_vect(6.0/scrx,u[i]);
      scrV = plcl_scale_vect(sig,u[i-1]);
      plcl_M_sub_vect(u[i],scrV);
      plcl_M_component_div(u[i],p);
 
    }

    /* qn = 0.5; */

    qn = 0.5;

    /* un = (3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1])); */

    scrx = spline_L->cp[comp].svals[i] - spline_L->cp[comp].svals[i-1];
    scrV = plcl_vect_diff(spline_L->cp[comp].vt[i],spline_L->cp[comp].vt[i-1]);

    plcl_M_vlincomb(un,3.0/scrx,ypn,-3.0/pow(scrx,2),scrV);

    /* y2[n] = (un - qn*u[n-1])/(qn*y2[n-1]+1.0); */

    scrV.c[0] = scrV.c[1] = scrV.c[2] = 1.0;
    plcl_M_vlincomb(scrW,qn,spline_L->cp[comp].vt2[i-1],1.0,scrV);
    plcl_M_vlincomb(scrV,1.0,un,-qn,u[i-1]);

    spline_L->cp[comp].vt2[i] = plcl_component_div(scrV,scrW);
                     
    /* for (k=n-1;k>=1;k--) { */

    for (k=i-1;k>=0;k--) {
      
      /* y2[k] = y2[k]*y2[k+1] + u[k]; */

      plcl_M_component_mult(spline_L->cp[comp].vt2[k],
                            spline_L->cp[comp].vt2[k+1]);
      plcl_M_add_vect(spline_L->cp[comp].vt2[k],u[k]);

    }

    /* free(u); */

    free(u);

  }

  return spline_L;
  /* Note to self: color splining is not yet implemented! */
}

/*
 * Converts spline back to regular plCurve, but changes the number of verts
 * to those in nv. We require that the number of components in nv match those
 * in spL. 
 *
 * A new plCurve is allocated. 
 */
plCurve *plCurve_convert_from_spline(const plCurve_spline * const spL,
                                     const int * const nv)
{
  int *cc; 
  bool *open;
  int comp, i;
  double s,s_step;
  int shi, slo;
  double h, b, a;

  plcl_vector scrV, scrW;
  plCurve *L;

  /* First, we do some elementary sanity checking. */

  assert(spL != NULL);
  assert(spL->nc >= 1);

  /* Now we allocate the new plCurve. */

  cc = (int *)malloc(spL->nc*sizeof(int));
  assert(cc != NULL);
  open = (bool *)malloc(spL->nc*sizeof(bool));
  assert(open != NULL);

  /*@+loopexec@*/
  for (i=0; i<spL->nc; i++) { 
    open[i] = spL->cp[i].open; 
    cc[i] = spL->cp[i].cc; 
  }
  /*@=loopexec@*/

  L = plCurve_new(spL->nc,nv,open,cc);
  
  if (L == NULL || plcl_error_num != 0) { return NULL; }

  free(open); free(cc);

  /* Now we can start filling in vertex positions, component by component. */

  for(comp=0;comp < L->nc;comp++) {

    s_step = spL->cp[comp].svals[spL->cp[comp].ns] / 
      (L->cp[comp].open ? L->cp[comp].nv-1 : L->cp[comp].nv);

    for(i=0,s=0,slo=0,shi=1;i<L->cp[comp].nv;i++,s+=s_step) {

      /* First, search to make sure that slo and shi bracket our current s.*/

      s = (s > spL->cp[comp].svals[spL->cp[comp].ns] ? 
           spL->cp[comp].svals[spL->cp[comp].ns] : s);

      while(spL->cp[comp].svals[shi] < s && slo < spL->cp[comp].ns) 
        { shi++; slo++; }

      /* Now we can evaluate the spline polynomial. */

      h = spL->cp[comp].svals[shi] - spL->cp[comp].svals[slo];
      assert(h > DBL_EPSILON);

      a = (spL->cp[comp].svals[shi] - s)/h;
      b = (s - spL->cp[comp].svals[slo])/h;

      plcl_M_vlincomb(scrV,a,spL->cp[comp].vt[slo],b,spL->cp[comp].vt[shi]);

      plcl_M_vlincomb(scrW, (a*a*a - a)*(h*h)/6.0, spL->cp[comp].vt2[slo],
                            (b*b*b - b)*(h*h)/6.0, spL->cp[comp].vt2[shi]);

      L->cp[comp].vt[i] = plcl_vect_sum(scrV,scrW);
    }

    /* Now we copy color data back to the buffer. */
    /* This is probably meaningless, unless cc == 1. */

    for (i=0;i<L->cp[comp].cc;i++) {
      L->cp[comp].clr[i] = spL->cp[comp].clr[i];
    }
  }
  return L;
}


/* Procedure samples the spline at a particular s value, returning a
 * spatial position. */
plcl_vector plCurve_sample_spline(const plCurve_spline * const spL,
                                 const int cmp,
                                 double s)
{
  int klo,khi,k;
  double cmpLen;
  double h, b, a;

  plcl_vector scrV, scrW;
  plcl_vector retV;

  /* We begin with a bit of checking to make sure that cmp and s seem
     compatible with the given spL. */

  assert(cmp >= 0);
  assert(cmp < spL->nc);

  /* Now fix any wraparound, so that the s value given is in [0,cmpLen). */

  cmpLen = spL->cp[cmp].svals[spL->cp[cmp].ns];
  while (s < 0) {    s += cmpLen; }
  while (s >= cmpLen) { s -= cmpLen; }

  /* We now search for the correct s value. */

  klo = 0;
  khi = spL->cp[cmp].ns;

  while (khi-klo > 1) {

    k = (khi+klo) >> 1;
    if (spL->cp[cmp].svals[k] > s) khi=k;
    else klo = k;

  }

  /* Now we can evaluate the polynomial. */
  
  h = spL->cp[cmp].svals[khi] - spL->cp[cmp].svals[klo];
  assert(h > DBL_EPSILON);

  a = (spL->cp[cmp].svals[khi] - s)/h;
  b = (s - spL->cp[cmp].svals[klo])/h;
  
  plcl_M_vlincomb(scrV,a,spL->cp[cmp].vt[klo],b,spL->cp[cmp].vt[khi]);
  
  plcl_M_vlincomb(scrW,
    (a*a*a - a)*(h*h)/6.0,spL->cp[cmp].vt2[klo],
    (b*b*b - b)*(h*h)/6.0,spL->cp[cmp].vt2[khi]);
  
  retV = plcl_vect_sum(scrV,scrW);  
  return retV;
}

/*
 * Take a plCurve, make sure it satisfies its constraints, convert it to a
 * spline and back, returning a plCurve with the desired number of points (or
 * as close to it as can be, given the constraints.
 *
 */
plCurve *meliorate(plCurve *L) {
  return L;
} 
