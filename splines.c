/*
 * Routines to create, destroy, and convert spline equivalents of plCurves
 *
 * $Id: splines.c,v 1.13 2006-02-15 22:39:19 ashted Exp $
 *
 * This code generates refinements of links, component by component, using the
 * Numerical Recipes spline code for interpolation. 
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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

#include "plCurve.h"

static inline void spline_pline_new(linklib_spline_pline *Pl,
                                                     int  ns, 
                                                     int  open, 
                                                     int  cc) {
 
  if (ns < 1) {
    plcl_error_num = PLCL_E_TOO_FEW_SAMPS;
    sprintf(plcl_error_str,
      "spline_pline_new: Can't create a spline_pline with %d samples.\n",ns);
    return;
  }

  Pl->open = open;
  Pl->ns = ns;

  if ((Pl->svals = (double *)malloc((ns+2)*sizeof(double))) == NULL) {
    plcl_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plcl_error_str,"spline_pline_new: Can't allocation space for %d"
            " samples in spline_pline_new.\n",ns);
    return;
  }
  Pl->svals++; /* so that Pl->svals[-1] is valid. */

  if ((Pl->vt = 
       (plcl_vector *)malloc((ns+2)*sizeof(plcl_vector))) == NULL) {
    plcl_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plcl_error_str,"spline_pline_new: Can't allocate space for %d "
      "samples in spline_pline_new.\n",ns);
    return;
  }
  Pl->vt++; /* so that Pl->vt[-1] is a valid space */

  if ((Pl->vt2 = 
       (plcl_vector *)malloc((ns+2)*sizeof(plcl_vector))) == NULL) {
    plcl_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plcl_error_str,"spline_pline_new: Can't allocate space for %d "
      "samples in spline_pline_new.\n",ns);
    return;
  }
  Pl->vt2++; /* so that Pl->vt2[-1] is a valid space */

  Pl->cc = cc;
  if ((Pl->clr = (plCurve_color *)malloc(cc*sizeof(plCurve_color))) == NULL) {
    plcl_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plcl_error_str,"spline_pline_new: Can't allocate space for %d "
      "colors in pline_new.\n",cc);
    return;
  }
  if ((Pl->clr2 = (plCurve_color *)malloc(cc*sizeof(plCurve_color))) == NULL) {
    plcl_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plcl_error_str,"spline_pline_new: Can't allocate space for %d "
      "colors in pline_new.\n",cc);
    return;
  }
}

/*
 * Procedure allocates memory for a new spline_link. The number of
 * components is given by components. The number of data samples in each
 * component shows up in the buffer pointed to by ns.  The closed/open
 * nature of each pline is given in the array pointed to by open.
 *
 */
plCurve_spline *linklib_spline_link_new(const int components, 
                                        const int * const ns, 
                                        const int * const open, 
                                        const int * const cc) 
{
  plCurve_spline *L;
  int i;

  /* First, we check to see that the input values are reasonable. */

  plcl_error_num = plcl_error_str[0] = 0;

  if (components < 1) {
    plcl_error_num = PLCL_E_TOO_FEW_COMPS;
    sprintf(plcl_error_str,"linklib_spline_link_new: Can't create a link "
      "with %d components.",components);
    return NULL;
  }

  if (ns == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
#ifdef HAVE_STRLCPY
    strlcpy(plcl_error_str,"plCurve_new: ns is NULL.\n",
      sizeof(plcl_error_str));
#else
    strlncpy(plcl_error_str,"plCurve_new: ns is NULL.\n",
      sizeof(plcl_error_str)-1);
    plcl_error_str[sizeof(plcl_error_str)-1] = '\0';
#endif
    return NULL;
  }
  if (open == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
#ifdef HAVE_STRLCPY
    strlcpy(plcl_error_str,"plCurve_new: open is NULL.\n",
      sizeof(plcl_error_str));
#else
    strlncpy(plcl_error_str,"plCurve_new: open is NULL.\n",
      sizeof(plcl_error_str)-1);
    plcl_error_str[sizeof(plcl_error_str)-1] = '\0';
#endif
    return NULL;
  }
  if (cc == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
#ifdef HAVE_STRLCPY
    strlcpy(plcl_error_str,"plCurve_new: cc is NULL.\n",
      sizeof(plcl_error_str));
#else
    strlncpy(plcl_error_str,"plCurve_new: cc is NULL.\n",
      sizeof(plcl_error_str)-1);
    plcl_error_str[sizeof(plcl_error_str)-1] = '\0';
#endif
    return NULL;
  }

  /* Now we attempt to allocate space for these components. */
  
  if ((L = malloc(sizeof(plCurve_spline))) == NULL) {
    plcl_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plcl_error_str,"linklib_spline_link_new: Could not allocate "
      "space for link in link_new.\n");
    return NULL;
  }
  L->nc = components;
  if ((L->cp = malloc(L->nc*sizeof(linklib_spline_pline))) == NULL) {
    plcl_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plcl_error_str,"Can't allocate array of spline_pline ptrs "
      "in link_new.\n");
    return NULL;
  }

  for (i = 0; i < L->nc; i++) {
    spline_pline_new(&L->cp[i],ns[i],open[i],cc[i]);
  }

  return L;
}

/*
 * Free the memory used to hold vertices in a given pline (not the memory of
 * the pline itself).  We then set all the values in the link data structure to
 * reflect the fact that the memory has been freed.  We can call pline_free
 * twice on the same pline without fear. 
 *
 */ 
static inline void spline_pline_free(linklib_spline_pline *Pl) {
  
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

} /* pline_free */

/*
 * Free the memory associated with a given link.  We then set all the values in
 * the link data structure to reflect the fact that the memory has been freed.
 * We can call link_free twice on the same link without fear. 
 *
 */ 
void linklib_spline_link_free(plCurve_spline *L) {
  int i;

  plcl_error_num = plcl_error_str[0] = 0;

  /* First, we check the input. */
  if (L == NULL) {
    return; /* Move along, nothing to see here */
  }

  /* Now we can get to work. */
  for (i=0; i<L->nc; i++) {
    spline_pline_free(&L->cp[i]);
  }

  free(L->cp);
  L->nc = 0;

  free(L);
  L = NULL;
} /* linklib_spline_link_free */


/*
 * Converts a regular link to spline link form. The spline code is adapted from
 * the "Numerical recipes" spline code. Each line from NR is commented next to
 * the (hopefully) equivalent vectorized version below. 
 */
plCurve_spline *convert_to_spline_link(plCurve * const L)
{
  int    i;
  int    *ns,*cc,*open;
  int    comp;

  plCurve_spline *spline_L;

  /* First, we check for sanity */

  plcl_error_num = plcl_error_str[0] = 0;

  if (L == NULL) {
    plcl_error_num = PLCL_E_NULL_PTR;
    sprintf(plcl_error_str,"convert_to_spline_link: Passed NULL "
            "pointer.\n");
    return NULL;
  }

  if (L->nc < 1) {
    plcl_error_num = PLCL_E_TOO_FEW_COMPS;
    sprintf(plcl_error_str,"convert_to_spline_link: Link pointer "
            "appears corrupted (%d components).\n",L->nc);
    return NULL;
  }

  plCurve_fix_wrap(L);

  ns = calloc(L->nc,sizeof(int));
  open = calloc(L->nc,sizeof(int));
  cc = calloc(L->nc,sizeof(int));

  if (ns == NULL || open == NULL || cc == NULL) {
    plcl_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plcl_error_str,"convert_to_spline_link: Couldn't allocate"
            " memory for buffers to make %d component spline_link.",L->nc);
    return NULL;
  } 

  for(i=0;i<L->nc;i++) {
    ns[i] = L->cp[i].nv;
    cc[i] = L->cp[i].cc;
    open[i] = L->cp[i].open;
  }

  /* Now we allocate the new spline link. */

  spline_L = linklib_spline_link_new(L->nc,ns,open,cc);

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

    if (u == NULL) {

      plcl_error_num = PLCL_E_CANT_ALLOC;
      sprintf(plcl_error_str,"convert_to_spline_link: Unable to allocate"
              " splining buffer for %d verts.\n",spline_L->cp[comp].ns+2);
      return NULL;

    }

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
 * Converts spline_link back to regular link, but changes the number of verts
 * to those in nv. We require that the number of components in nv match those
 * in spL. 
 *
 * A new plCurve is allocated. 
 */
plCurve *convert_spline_to_link(const plCurve_spline * const spL,
                                const int * const nv)
{
  int *cc, *open;
  int comp, i;
  double s,s_step;
  int shi, slo;
  double h, b, a;

  plcl_vector scrV, scrW;
  plCurve *L;

  /* First, we do some elementary sanity checking. */

  plcl_error_num = plcl_error_str[0] = 0;

  if (spL == NULL) {

    plcl_error_num = PLCL_E_NULL_PTR;
    sprintf(plcl_error_str,"convert_spline_to_link: Passed NULL pointer.\n");
    return NULL;

  }

  if (spL->nc < 1) {
    plcl_error_num = PLCL_E_TOO_FEW_COMPS;
    sprintf(plcl_error_str,
            "convert_spline_to_link: spline_link appears corrupt (nc = %d).\n",
            spL->nc);
    return NULL;
  }

  /* Now we allocate the new (conventional) link. */

  cc = malloc(spL->nc*sizeof(int));
  open = malloc(spL->nc*sizeof(int));

  if (cc == NULL || open == NULL) {

    plcl_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plcl_error_str,
            "convert_spline_to_link: Couldn't allocate cc, open buffer"
            " of size %d.\n",spL->nc);
    return NULL;

  }

  for(i=0;i<spL->nc;i++) { open[i] = spL->cp[i].open; cc[i] = spL->cp[i].cc; }

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

      if (h < 10e-14) { 
        
        plcl_error_num = PLCL_E_SMP_TOO_CLOSE;
        sprintf(plcl_error_str,"convert_spline_to_link: svals of samples %d "
                "and %d are"
                " too close for spline algorithm (svals are %g and %g).\n",
                shi,slo,spL->cp[comp].svals[shi],spL->cp[comp].svals[slo]);
        return NULL;

      }

      a = (spL->cp[comp].svals[shi] - s)/h;
      b = (s - spL->cp[comp].svals[slo])/h;

      plcl_M_vlincomb(scrV,a,spL->cp[comp].vt[slo],b,spL->cp[comp].vt[shi]);

      plcl_M_vlincomb(scrW,
        (a*a*a - a)*(h*h)/6.0, spL->cp[comp].vt2[slo],
        (b*b*b - b)*(h*h)/6.0, spL->cp[comp].vt2[shi]);

      L->cp[comp].vt[i] = plcl_vect_sum(scrV,scrW);

    }

    /* Now we copy color data back to the buffer. */
    /* This is probably meaningless, unless cc == 1. */

    for(i=0;i<L->cp[comp].cc;i++) {

      L->cp[comp].clr[i] = spL->cp[comp].clr[i];

    }

  }

  return L;

}


/* Procedure evaluates the spline link at a particular s value, 
   returning a spatial position. */
plcl_vector evaluate_spline_link(const plCurve_spline * const spL,
                                 const int cmp,
                                 double s)
{
  int klo,khi,k;
  double cmpLen;
  double h, b, a;

  plcl_vector scrV, scrW;
  plcl_vector retV;
  plcl_vector zeroVec = {{0,0,0}};

  /* We begin with a bit of checking to make sure that cmp and s seem
     compatible with the given spL. */

  plcl_error_num = plcl_error_str[0] = 0;

  if (cmp < 0 || cmp > spL->nc) {
    plcl_error_num = PLCL_E_CANT_FIND_POS;
    sprintf(plcl_error_str,"evaluate_spline_link: Can't find position %g"
            " on component %d of the %d component link spL.\n",
            s,cmp,spL->nc);
    return zeroVec;
  }

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

  if (h < 10e-14) { 
    
    plcl_error_num = PLCL_E_SMP_TOO_CLOSE;
    sprintf(plcl_error_str,"evaluate_spline_link: svals of samples %d "
            "and %d are"
            " too close for spline algorithm (svals are %g and %g).\n",
            khi,klo,spL->cp[cmp].svals[khi],spL->cp[cmp].svals[klo]);
    return zeroVec;

  }
  
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
