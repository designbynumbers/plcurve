/*
 *
 * Data structures and prototypes for spline_links
 *
 *  $Id: spline_links.h,v 1.3 2006-02-03 13:10:20 ashted Exp $
 *
 */

/* Copyright 2004 The University of Georgia */

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

#ifndef __SPLINE_LINK_H
#define __SPLINE_LINK_H

#if (__cplusplus || c_plusplus)
extern "C" {
#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
  
#ifdef HAVE_MATH_H
#include <math.h>
#endif
  
#include "plCurve.h"
  
  /* This code generates refinements of links, component by component,
     using the Numerical Recipes spline code for interpolation. */
  
  typedef struct linklib_spline_pline_type {
    int          open;   /* This is an "open" pline (with distinct ends) */
    int          ns;     /* Number of samples used to build spline. */
    double      *svals;  /* s values at samples */
    plcl_vector *vt;     /* positions at these s values */
    plcl_vector *vt2;    /* second derivatives at these s values */
    
    /***** Need a way to specify constraints on endpoints here *****/
    
    int             cc;
    plCurve_color  *clr;       /* color values at samples */
    plCurve_color  *clr2;      /* second derivatives at these s values */
    
  } linklib_spline_pline;
  
  typedef struct linklib_spline_link_type {	
    int nc;			/* Number of components */
    linklib_spline_pline *cp;     /* Components */
  } linklib_spline_link;
  

  /* Allocate new spline_link. */
  linklib_spline_link *linklib_spline_link_new(int components, 
					       const int *ns, 
					       const int *open, 
					       const int *cc);
  
  /* Free memory for spline_link. */
  void linklib_spline_link_free(linklib_spline_link *L);

  /* Convert conventional link to spline_link. */
  linklib_spline_link *convert_to_spline_link(plCurve *L);

  /* Convert spline_link to conventional link (with resampling). */
  plCurve *convert_spline_to_link(linklib_spline_link *spL,int *nv);

  /* Evaluate a spline_link at a particular s value. */
  plcl_vector evaluate_spline_link(linklib_spline_link *spL,int cmp,double s);

#if (__cplusplus || c_plusplus)
};
#endif
#endif

