/*
 *  plwhitten.c : Perform a whitten group operation on a plCurve. 
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

#include <stdlib.h>
#include <math.h>

#ifdef HAVE_MALLOC_H
  #include <malloc.h>
#endif

#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <float.h>
 
/* Performs a mirror/reverse/permute operation from the Whitten group on L. */
/* We expect input in the form:
   
   mirror = +-1,
   epsilon = array of L->nc +-1
   perm = array of 2*L->nc 0-denominated indices
   
   so that if perm[i] is the first occurrence of 
   permutation j then perm[i+1] is p(j). 

   The transformation rule is that for a general element in the form

   e_0 x (e_1 x ... x e_n semidirect p)

   we take

   (K_1,...K_n) -> e_0 * (e_1 K_p(1),...,e_n K_p(n))

   where -K_i is K_i with orientations reversed, and -1 * K_i is the 
   mirror image of K_i.

   The operation works in place on L.

*/

int  perm(int permsize,int *perm,int ind)
{
  int i;

  for(i=0;i<permsize;i++) {
    
    if (perm[i] == ind) { return perm[i+1]; }

  }

  return ind;
}

int  iperm(int permsize,int *perm,int ind)
{
  int i;

  for(i=permsize;i>0;i--) {
    
    if (perm[i] == ind) { return perm[i-1]; }

  }

  return ind;
}

void plc_whitten(plCurve *L, int mirror, int *eps, int *thisperm)

{
  plc_strand *newpos;
  int i;

  /* First, create a new, permuted version of the strands in L */

  newpos = calloc(L->nc,sizeof(plc_strand));

  for(i=0;i<L->nc;i++) {

    newpos[i] = L->cp[perm(2*L->nc,thisperm,i)];

  }

  /* Now recopy into L's cp array */

  for (i=0;i<L->nc;i++) {
    
    L->cp[i] = newpos[i];

  }

  /* and get rid of newpos */

  free(newpos);

  /* Of course, we need to change the cmp indices in constraints, too.     */
  /* Here we observe that a constraint which applied to component p(i) now */
  /* applies to component i, so we apply the _inverse_ permutation to tc->cmp */

  plc_constraint *tc;

  for(tc=L->cst;tc != NULL;tc=tc->next) {

    tc->cmp = iperm(2*L->nc,thisperm,tc->cmp);

  }

  /* Now perform the mirror and reverse operations */

  if (mirror == -1) {

    int cp, vt;

    /* Mirror over the xy plane */

    for(cp=0;cp<L->nc;cp++) {

      for(vt=0;vt<L->cp[cp].nv;vt++) {

	L->cp[cp].vt[vt].c[2] *= -1;

      }

    }

    /* Now mirror the constraints */

    plc_constraint *tc;

    for(tc=L->cst;tc != NULL;tc = tc->next) {

      tc->vect[0].c[2] *= -1;  /* For each type, we reflect this... */

      if (tc->kind == line) {

	tc->vect[1].c[2] *= -1; /* This is a point on the line... */

      }

      /* For a plane constraint, the third vector is a distance, which 
	 is preserved under reflection, so there's nothing to do. */

    }
    
  }

  plc_fix_wrap(L); // We changed vertices, so we fix the wrap.

  /* Having mirrored, we now reverse the components that need it. */

  int cp,vt;

  for(cp=0;cp<L->nc;cp++) {

    if (eps[cp] == -1) { /* We should reverse this guy... */

      plc_vector swapv;
      plc_color  swapc;
      int nv;

      nv = L->cp[cp].nv;

      for(vt=0;vt<nv/2;vt++) {

	/* We always swap the vertices */

	swapv = L->cp[cp].vt[vt];
	L->cp[cp].vt[vt] = L->cp[cp].vt[nv-vt-1];
	L->cp[cp].vt[nv-vt-1] = swapv;

	
      }

      if (L->cp[cp].cc == nv) { 
	/* We have per-vertex colors and need to swap 'em */
	
	for(vt=0;vt<nv/2;vt++) {

	  swapc = L->cp[cp].clr[vt];
	  L->cp[cp].clr[vt] = L->cp[cp].clr[nv-vt-1];
	  L->cp[cp].clr[nv-vt-1] = swapc;
	  	
	}

      } /* The other possibilities for cc are 0 or 1, 
	   which don't need swapping. */

      /* Now we handle the changed vertex indices in constraints. */

      for(tc=L->cst;tc != NULL;tc=tc->next) {

	tc->vert = (nv-1) - (tc->vert + tc->num_verts); 

      }

    }

  }

  plc_fix_wrap(L);

  /* We have now applied everything. */

}
