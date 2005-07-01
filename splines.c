/*
 *  Routines to create, destroy, and convert spline_links (and spline_plines)
 * 
 *  $Id: splines.c,v 1.2 2005-07-01 01:56:33 cantarel Exp $
 *
 */

/* Copyright 2004 The University of Georgia. */

/* This file is part of vecttools.
   
vecttools is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

vecttools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vecttools; if not, write to the Free Software
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

#include <link.h>
#include <spline_links.h>

extern int  linklib_error_num;
extern char linklib_error_str[80];

static void spline_pline_new(linklib_spline_pline *Pl,int ns, int acyclic, int cc) {
 
  if (ns < 1) {
    linklib_error_num = 71;
    sprintf(linklib_error_str,"spline_pline_new: Can't create a spline_pline with %d samples.\n",ns);
    return;
  }

  Pl->acyclic = acyclic;
  Pl->ns = ns;

  if ((Pl->svals = (double *)malloc((ns+2)*sizeof(double))) == NULL) {

    linklib_error_num = 317;
    sprintf(linklib_error_str,"spline_pline_new: Can't allocation space for %d"
	    " samples in spline_pline_new.\n",ns);
    return;
  }
  Pl->svals++; /* so that Pl->svals[-1] is valid. */

  if ((Pl->vt = 
       (linklib_vector *)malloc((ns+2)*sizeof(linklib_vector))) == NULL) {
    linklib_error_num = 72;
    sprintf(linklib_error_str,"spline_pline_new: Can't allocate space for %d samples in spline_pline_new.\n",ns);
    return;
  }
  Pl->vt++; /* so that Pl->vt[-1] is a valid space */

  if ((Pl->vt2 = 
       (linklib_vector *)malloc((ns+2)*sizeof(linklib_vector))) == NULL) {
    linklib_error_num = 72;
    sprintf(linklib_error_str,"spline_pline_new: Can't allocate space for %d samples in spline_pline_new.\n",ns);
    return;
  }
  Pl->vt2++; /* so that Pl->vt2[-1] is a valid space */

  Pl->cc = cc;
  if ((Pl->clr = (linklib_color *)malloc(cc*sizeof(linklib_color))) == NULL) {
    linklib_error_num = 73;
    sprintf(linklib_error_str,"spline_pline_new: Can't allocate space for %d colors in pline_new.\n",cc);
    return;
  }
  if ((Pl->clr2 = (linklib_color *)malloc(cc*sizeof(linklib_color))) == NULL) {
    linklib_error_num = 73;
    sprintf(linklib_error_str,"spline_pline_new: Can't allocate space for %d colors in pline_new.\n",cc);
    return;
  }
}

/*
 * Procedure allocates memory for a new spline_link. The number of
 * components is given by components. The number of data samples in each
 * component shows up in the buffer pointed to by ns.  The closed/open
 * nature of each pline is given in the array pointed to by acyclic.
 *
 */
linklib_spline_link *linklib_spline_link_new(int components, 
					     const int *ns, 
					     const int *acyclic, 
					     const int *cc) 
{
  linklib_spline_link *L;
  int i;

  /* First, we check to see that the input values are reasonable. */

  if (components < 1) {
    linklib_error_num = 81;
    sprintf(linklib_error_str,"linklib_spline_link_new: Can't create a link with %d components.",components);
    return NULL;
  }

  if (ns == NULL || acyclic == NULL || cc == NULL) {
    linklib_error_num = 82;
    sprintf(linklib_error_str,"linklib_spline_link_new: "
	    "ns, acyclic or cc is NULL.");
    return NULL;
  }

  /* Now we attempt to allocate space for these components. */
  
  if ((L = (linklib_spline_link *)malloc(sizeof(linklib_spline_link))) == NULL) {
    linklib_error_num = 83;
    sprintf(linklib_error_str,"linklib_spline_link_new: Could not allocate space for link in link_new.\n");
    return NULL;
  }
  L->nc = components;
  if ((L->cp = (linklib_spline_pline *)malloc(L->nc*sizeof(linklib_spline_pline))) == NULL) {
    linklib_error_num = 84;
    sprintf(linklib_error_str,"Can't allocate array of spline_pline ptrs in link_new.\n");
    return NULL;
  }

  for (i = 0; i < L->nc; i++) {
    spline_pline_new(&L->cp[i],ns[i],acyclic[i],cc[i]);
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
void spline_pline_free(linklib_spline_pline *Pl) {
  
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
void linklib_spline_link_free(linklib_spline_link *L) {
  int i;

  /* First, we check the input. */
  if (L == NULL) {
    return; /* Move along, nothing to see here */
  }

  if (L->nc < 0) {
    linklib_error_num = 91;
    sprintf(linklib_error_str,"linklib_spline_link_free: Link appears corrupted. L.nc = %d.",L->nc);
    return;
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


linklib_spline_link *convert_to_spline_link(linklib_link *L)

     /* Converts a regular link to spline link form. The spline */
     /* code is adapted from the "Numerical recipes" spline code. */
     /* Each line from NR is commented next to the (hopefully) */
     /* equivalent vectorized version below. */
{
  int    i;
  int    *ns,*cc,*acyclic;
  int    comp;

  linklib_spline_link *spline_L;

  /* First, we check for sanity */

  if (L == NULL) {

    linklib_error_num = 519;
    sprintf(linklib_error_str,"convert_to_spline_link: Passed NULL "
	    "pointer.\n");
    return NULL;

  }

  if (L->nc < 1) {

    linklib_error_num = 520;
    sprintf(linklib_error_str,"convert_to_spline_link: Link pointer "
	    "appears corrupted (%d components).\n",L->nc);
    return NULL;

  }

  linklib_link_fix_wrap(L);

  ns = calloc(L->nc,sizeof(int));
  acyclic = calloc(L->nc,sizeof(int));
  cc = calloc(L->nc,sizeof(int));

  if (ns == NULL || acyclic == NULL || cc == NULL) {

    linklib_error_num = 521;
    sprintf(linklib_error_str,"convert_to_spline_link: Couldn't allocate"
	    " memory for buffers to make %d component spline_link.",L->nc);
    return NULL;

  } 

  for(i=0;i<L->nc;i++) {

    ns[i] = L->cp[i].nv;
    cc[i] = L->cp[i].cc;
    acyclic[i] = L->cp[i].acyclic;

  }

  /* Now we allocate the new spline link. */

  spline_L = linklib_spline_link_new(L->nc,ns,acyclic,cc);

  if (linklib_error_num != 0 || spline_L == NULL) {

    return NULL;

  }

  free(ns); free(acyclic); free(cc);
	  
  /* We now go component-by-component. */

  for(comp = 0;comp < L->nc;comp++) {

    /* Our first task is to assemble svals: the arclength of each vertex. */

    spline_L->cp[comp].svals[0] = 0.0;

    for(i=1;i<L->cp[comp].nv;i++) {

      spline_L->cp[comp].svals[i] = 
	spline_L->cp[comp].svals[i-1] + linklib_vdist(L->cp[comp].vt[i-1],
					    L->cp[comp].vt[i]);
    }
    
    if (!L->cp[comp].acyclic) { 

      spline_L->cp[comp].svals[i] = 
	spline_L->cp[comp].svals[i-1] + linklib_vdist(L->cp[comp].vt[i-1],
						      L->cp[comp].vt[i]);     

      spline_L->cp[comp].svals[-1] = -linklib_vdist(L->cp[comp].vt[0],
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

    linklib_vector p,un,*u;
    double sig,qn;

    /* u = malloc(n,sizeof(double)); */

    u = malloc((spline_L->cp[comp].ns+2)*sizeof(linklib_vector));

    if (u == NULL) {

      linklib_error_num = 522;
      sprintf(linklib_error_str,"convert_to_spline_link: Unable to allocate"
	      " splining buffer for %d verts.\n",spline_L->cp[comp].ns+2);
      return NULL;

    }

    /* ...together with values yp1 and ypn for the first derivative at 
       the first and last samples... */

    linklib_vector yp1,ypn;

    yp1 = linklib_link_tangent_vector(L,comp,0);

    if (L->cp[comp].acyclic) { 

      ypn = linklib_link_tangent_vector(L,comp,L->cp[comp].nv-1);

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
    linklib_vector scrV,scrW;

    scrV = linklib_vminus(spline_L->cp[comp].vt[1],spline_L->cp[comp].vt[0]);
    scrx = spline_L->cp[comp].svals[1] - spline_L->cp[comp].svals[0];
    linklib_vlincombine(scrV,3.0/pow(scrx,2),yp1,-3.0/scrx,u[0]);
    
    /* for (i=2;i<=n-1;i++) { */

    I = L->cp[comp].acyclic ? L->cp[comp].nv-2 : L->cp[comp].nv-1;

    for(i=1;i<=I;i++) {

      /* sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]); */

      sig = (spline_L->cp[comp].svals[i] - spline_L->cp[comp].svals[i-1])/
	(spline_L->cp[comp].svals[i+1] - spline_L->cp[comp].svals[i-1]);

      /* p = sig*y2[i-1] + 2.0; */

      scrV.c[0] = scrV.c[1] = scrV.c[2] = 2.0;
      linklib_vlincombine(spline_L->cp[comp].vt2[i-1],sig,scrV,1.0,p);

      /* y2[i] = (sig-1.0)/p; */

      scrV.c[0] = scrV.c[1] = scrV.c[2] = sig - 1.0;
      spline_L->cp[comp].vt2[i] = linklib_vdivide(scrV,p);

      /* u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]); */

      scrV = linklib_vminus(spline_L->cp[comp].vt[i+1],spline_L->cp[comp].vt[i]);
      scrW = linklib_vminus(spline_L->cp[comp].vt[i],spline_L->cp[comp].vt[i-1]);

      linklib_vlincombine(scrV,
			  1.0/(spline_L->cp[comp].svals[i+1] - spline_L->cp[comp].svals[i]),
			  scrW,
			  -1.0/(spline_L->cp[comp].svals[i] - spline_L->cp[comp].svals[i-1]),
			  u[i]);

      /* u[i]=(6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p; */

      scrx = spline_L->cp[comp].svals[i+1] - spline_L->cp[comp].svals[i-1];
      linklib_vsmult(6.0/scrx,u[i]);
      scrV = linklib_scalarmult(sig,u[i-1]);
      linklib_vsub(u[i],scrV);
      linklib_vdiv(u[i],p);
 
    }

    /* qn = 0.5; */

    qn = 0.5;

    /* un = (3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1])); */

    scrx = spline_L->cp[comp].svals[i] - spline_L->cp[comp].svals[i-1];
    scrV = linklib_vminus(spline_L->cp[comp].vt[i],spline_L->cp[comp].vt[i-1]);

    linklib_vlincombine(ypn,3.0/scrx,scrV,-3.0/pow(scrx,2),un);

    /* y2[n] = (un - qn*u[n-1])/(qn*y2[n-1]+1.0); */

    scrV.c[0] = scrV.c[1] = scrV.c[2] = 1.0;
    linklib_vlincombine(spline_L->cp[comp].vt2[i-1],qn,scrV,1.0,scrW);
    linklib_vlincombine(un,1.0,u[i-1],-qn,scrV);

    spline_L->cp[comp].vt2[i] = linklib_vdivide(scrV,scrW);
		     
    /* for (k=n-1;k>=1;k--) { */

    for (k=i-1;k>=0;k--) {
      
      /* y2[k] = y2[k]*y2[k+1] + u[k]; */

      linklib_vmul(spline_L->cp[comp].vt2[k],spline_L->cp[comp].vt2[k+1]);
      linklib_vadd(spline_L->cp[comp].vt2[k],u[k]);

    }

    /* free(u); */

    free(u);

  }

  return spline_L;
  
  /* Note to self: color splining is not yet implemented! */

}

linklib_link *convert_spline_to_link(linklib_spline_link *spL,int *nv)

     /* Converts spline_link back to regular link, but changes the number 
	of verts to those in nv. We require that the number of components 
	in nv match those in spL. 

	A new linklib_link is allocated. */

{
  int *cc, *acyclic;
  int comp, i;
  double s,s_step;
  int shi, slo;
  double h, b, a;

  linklib_vector scrV, scrW;
  linklib_link *L;

  /* First, we do some elementary sanity checking. */

  if (spL == NULL) {

    linklib_error_num = 709;
    sprintf(linklib_error_str,"convert_spline_to_link: Passed NULL pointer.\n");
    return NULL;

  }

  if (spL->nc < 1) {

    linklib_error_num = 710;
    sprintf(linklib_error_str,
	    "convert_spline_to_link: spline_link appears corrupt (nc = %d).\n",spL->nc);
    return NULL;
  }

  /* Now we allocate the new (conventional) link. */

  cc = malloc(spL->nc*sizeof(int));
  acyclic = malloc(spL->nc*sizeof(int));

  if (cc == NULL || acyclic == NULL) {

    linklib_error_num = 711;
    sprintf(linklib_error_str,"convert_spline_to_link: Couldn't allocate cc, acyclic buffer"
	    " of size %d.\n",spL->nc);
    return NULL;

  }

  for(i=0;i<spL->nc;i++) { acyclic[i] = spL->cp[i].acyclic; cc[i] = spL->cp[i].cc; }

  L = linklib_link_new(spL->nc,nv,acyclic,cc);
  
  if (L == NULL || linklib_error_num != 0) { return NULL; }

  free(acyclic); free(cc);

  /* Now we can start filling in vertex positions, component by component. */

  for(comp=0;comp < L->nc;comp++) {

    s_step = spL->cp[comp].svals[spL->cp[comp].ns] / 
      (L->cp[comp].acyclic ? L->cp[comp].nv-1 : L->cp[comp].nv);

    for(i=0,s=0,slo=0,shi=1;i<L->cp[comp].nv;i++,s+=s_step) {

      /* First, search to make sure that slo and shi bracket our current s.*/

      s = (s > spL->cp[comp].svals[spL->cp[comp].ns] ? 
	   spL->cp[comp].svals[spL->cp[comp].ns] : s);

      while(spL->cp[comp].svals[shi] < s && slo < spL->cp[comp].ns) 
	{ shi++; slo++; }

      /* Now we can evaluate the spline polynomial. */

      h = spL->cp[comp].svals[shi] - spL->cp[comp].svals[slo];

      if (h < 10e-14) { 
	
	linklib_error_num = 713;
	sprintf(linklib_error_str,"convert_spline_to_link: svals of samples %d "
		"and %d are"
		" too close for spline algorithm (svals are %g and %g).\n",
		shi,slo,spL->cp[comp].svals[shi],spL->cp[comp].svals[slo]);
	return NULL;

      }

      a = (spL->cp[comp].svals[shi] - s)/h;
      b = (s - spL->cp[comp].svals[slo])/h;

      linklib_vlincombine(spL->cp[comp].vt[slo],a,spL->cp[comp].vt[shi],b,
			  scrV);

      linklib_vlincombine(spL->cp[comp].vt2[slo],(a*a*a - a)*(h*h)/6.0,
			  spL->cp[comp].vt2[shi],(b*b*b - b)*(h*h)/6.0,
			  scrW);

      L->cp[comp].vt[i] = linklib_vplus(scrV,scrW);

    }

    /* Now we copy color data back to the buffer. */
    /* This is probably meaningless, unless cc == 1. */

    for(i=0;i<L->cp[comp].cc;i++) {

      L->cp[comp].clr[i] = spL->cp[comp].clr[i];

    }

  }

  return L;

}

linklib_vector evaluate_spline_link(linklib_spline_link *spL,int cmp,double s)

/* Procedure evaluates the spline link at a particular s value, 
   returning a spatial position. */

{
  int klo,khi,k;
  double cmpLen;
  double h, b, a;

  linklib_vector scrV, scrW;
  linklib_vector retV;
  linklib_vector zeroVec = {{0,0,0}};

  /* We begin with a bit of checking to make sure that cmp and s seem
     compatible with the given spL. */

  if (cmp < 0 || cmp > spL->nc) {

    linklib_error_num = 1813;
    sprintf(linklib_error_str,"evaluate_spline_link: Can't find position %g"
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
    
    linklib_error_num = 3030;
    sprintf(linklib_error_str,"evaluate_spline_link: svals of samples %d "
	    "and %d are"
	    " too close for spline algorithm (svals are %g and %g).\n",
	    khi,klo,spL->cp[cmp].svals[khi],spL->cp[cmp].svals[klo]);
    return zeroVec;

  }
  
  a = (spL->cp[cmp].svals[khi] - s)/h;
  b = (s - spL->cp[cmp].svals[klo])/h;
  
  linklib_vlincombine(spL->cp[cmp].vt[klo],a,spL->cp[cmp].vt[khi],b,
		      scrV);
  
  linklib_vlincombine(spL->cp[cmp].vt2[klo],(a*a*a - a)*(h*h)/6.0,
		      spL->cp[cmp].vt2[khi],(b*b*b - b)*(h*h)/6.0,
		      scrW);
  
  retV = linklib_vplus(scrV,scrW);  
  return retV;

}
