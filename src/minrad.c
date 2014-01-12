/* 

     minrad.c : This file contains a single procedure that computes the
                "minrad" of a link L (in the sense of Eric Rawdon) and 
                returns it as a double.    

                We may later fold this into another source file.
                Note that we use some of the linear algebra stuff
                from vector.c. 

  $Id: minrad.c,v 1.18 2006-04-18 19:14:36 ashted Exp $

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
#ifdef HAVE_MATH_H
#include <math.h>
#endif
#ifdef HAVE_FLOAT_H
#include <float.h>
#endif
#include <plCurve.h>
#include "octrope.h"

#define mr_up_bound ((cutoff > 0) ? cutoff : mr + epsilon)

double octrope_minrad(plCurve *L, 
                      const double cutoff, const double epsilon,
                      octrope_mrloc *min_rad_locs,
                      int mr_size, int *num_min_rad_locs)
{
  int i,j,k,l;
  double mr = {DBL_MAX}, rad;
  plc_vector in,out;
  double normin, normout;
  double mrplus, mrminus, this_mr;
  double dot_prod,cross_prod_norm;
  int start,end;

  /* We start with some initializations. */
  
  octrope_error_num = octrope_error_str[0] = 0;
  plc_fix_wrap(L);
  if (num_min_rad_locs != NULL) *num_min_rad_locs = 0;
  
  /* Now we work. */

  for(i=0;i<L->nc;i++) {                  /* Loop over components. */

    if (L->cp[i].open) {
      start = 1;
      end = L->cp[i].nv - 1;
      out = plc_vect_diff(L->cp[i].vt[1],L->cp[i].vt[0]);      
      normout = plc_norm(out); 
    } else {
      start = 0;
      end = L->cp[i].nv;
      out = plc_vect_diff(L->cp[i].vt[0],L->cp[i].vt[-1]);      
      normout = plc_norm(out); 
    }

    /* Now we handle the main loop. */

    for (j = start; j < end; j++) {
      
      in     = out;                   /* Pass back "outgoing" data. */
      normin = normout;               
      
      out = plc_vect_diff(L->cp[i].vt[j+1],L->cp[i].vt[j]);
      normout = plc_norm(out);
      
      dot_prod = plc_M_dot(in,out);
      cross_prod_norm = plc_norm(plc_cross_prod(in,out));
      
      rad = (normin*normout + dot_prod)/(2*cross_prod_norm);

      if (cross_prod_norm < 1e-12) { /* Presumably, cross_prod_norm == 0 */

	continue;   /* We don't need to record this value, or indeed do 
		       anything else with it, since rad = Inf or NaN. 
		       So skip to the next iteration of the loop. */

      }
      
      mrminus = normin*rad;
      mrplus  = normout*rad;

      this_mr = (mrminus < mrplus) ? mrminus : mrplus;
      mr = (mr < this_mr) ? mr : this_mr;

      /* Now we need to enter one or both of mrminus and mrplus into the list. */

      if (num_min_rad_locs != NULL) { /* We are recording locations. */

	/* We start by trying to add mrminus */

        if (mrminus <= mr_up_bound) { /* This one is of interest. */
          if (*num_min_rad_locs < mr_size) { /* And we have room. */
            min_rad_locs[*num_min_rad_locs].component = i;
            min_rad_locs[*num_min_rad_locs].vert = j;
	    min_rad_locs[*num_min_rad_locs].svert = j-1;
            min_rad_locs[*num_min_rad_locs].mr = mrminus;
	    min_rad_locs[*num_min_rad_locs].compression = 0;
  
            (*num_min_rad_locs)++;
          } else { /* Or we can create room. */ 
            k = *num_min_rad_locs - 1;
            while ((k >= 0) && (min_rad_locs[k].mr > mr_up_bound)) { k--; }
            for (l = k; l >= 0; l--) {
              if (min_rad_locs[l].mr > mr_up_bound) {
                min_rad_locs[l] = min_rad_locs[k--];
              }
            }
            *num_min_rad_locs = k + 1;
            if (*num_min_rad_locs < mr_size) { /* We now have room. */
              min_rad_locs[*num_min_rad_locs].component = i;
              min_rad_locs[*num_min_rad_locs].vert = j;
	      min_rad_locs[*num_min_rad_locs].svert = j-1;
              min_rad_locs[*num_min_rad_locs].mr = mrminus;
	      min_rad_locs[*num_min_rad_locs].compression = 0;
	    
              (*num_min_rad_locs)++;
            }
          }
        }

	/* Now we try to add mrplus */

	if (mrplus <= mr_up_bound) { /* This one is of interest. */
          if (*num_min_rad_locs < mr_size) { /* And we have room. */
            min_rad_locs[*num_min_rad_locs].component = i;
            min_rad_locs[*num_min_rad_locs].vert = j;
	    min_rad_locs[*num_min_rad_locs].svert = j+1;
            min_rad_locs[*num_min_rad_locs].mr = mrplus;
	    min_rad_locs[*num_min_rad_locs].compression = 0;

            (*num_min_rad_locs)++;
          } else { /* Or we can create room. */ 
            k = *num_min_rad_locs - 1;
            while ((k >= 0) && (min_rad_locs[k].mr > mr_up_bound)) { k--; }
            for (l = k; l >= 0; l--) {
              if (min_rad_locs[l].mr > mr_up_bound) {
                min_rad_locs[l] = min_rad_locs[k--];
              }
            }
            *num_min_rad_locs = k + 1;
            if (*num_min_rad_locs < mr_size) { /* We now have room. */
              min_rad_locs[*num_min_rad_locs].component = i;
              min_rad_locs[*num_min_rad_locs].vert = j;
	      min_rad_locs[*num_min_rad_locs].svert = j+1;
              min_rad_locs[*num_min_rad_locs].mr = mrplus;
	      min_rad_locs[*num_min_rad_locs].compression = 0;
  
              (*num_min_rad_locs)++;
            }
          }
        }

	
      }
      
#ifdef DEBUG
      if (octrope_debug_level() > 8) {
        printf("%d/%d & %d/%d:\n",j,j-1,j+1,j);
        printf("  (%3.3f,%3.3f,%3.3f)-(%3.3f,%3.3f,%3.3f)=(%3.3f,%3.3f,%3.3f) |%3.3f|\n",
          L->cp[i].vt[j].c[0], L->cp[i].vt[j].c[1], L->cp[i].vt[j].c[2],
          L->cp[i].vt[j-1].c[0], L->cp[i].vt[j-1].c[1], L->cp[i].vt[j-1].c[2],
          in.c[0],in.c[1],in.c[2],normin);
        printf("  (%3.3f,%3.3f,%3.3f)-(%3.3f,%3.3f,%3.3f)=(%3.3f,%3.3f,%3.3f) |%3.3f|\n",
          L->cp[i].vt[j+1].c[0], L->cp[i].vt[j+1].c[1], L->cp[i].vt[j+1].c[2],
          L->cp[i].vt[j].c[0], L->cp[i].vt[j].c[1], L->cp[i].vt[j].c[2],
          out.c[0],out.c[1],out.c[2],normout);
        printf("  rad: %3.3f, mrminus: %3.3f, mrplus: %3.3f, mr: %3.3f\n",
          rad,mrminus,mrplus,mr);
      }
#endif /* DEBUG */
    }
    
  }

#ifdef DEBUG
   if (octrope_debug_level() > 5) {
     printf(
       "Paring minrad list in octrope_minrad, using bound: %f\n",mr_up_bound);
   }
#endif /* DEBUG */
  /* Pare down the list to be only the mr's within epsilon of minimum */
  if (num_min_rad_locs != NULL) {
    /* Note that by the way we built the list, all of the minrads 
       with value <= mr+epsilon will be at the end of the list. */
    k = *num_min_rad_locs - 1;
    while ((k >= 0) && (min_rad_locs[k].mr > mr_up_bound)) { k--; }
    for (l = k; l >= 0; l--) {
      if (min_rad_locs[l].mr > mr_up_bound) {
        min_rad_locs[l] = min_rad_locs[k--];
      }
    }
    *num_min_rad_locs = k + 1;
  }
  
  return mr;
}
