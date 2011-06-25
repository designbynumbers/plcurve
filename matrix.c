/*
 * @COPYRIGHT@
 *
 * Routines for working with matrices.
 *
 * $Id: matrix.c,v 1.37 2007-07-12 15:27:52 cantarel Exp $
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
#include <plCurve.h>

#ifdef HAVE_MATH_H
  #include <math.h>
#endif
#ifdef HAVE_STRING_H
  #include <string.h>
#endif
#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif
#ifdef HAVE_FLOAT_H
  #include <float.h>
#endif
#ifdef HAVE_ASSERT_H
  #include <assert.h>
#endif

plc_matrix  *plc_build_matrix(plc_vector row0,plc_vector row1,plc_vector row2)
     /* Builds a matrix from rows. */
{
  plc_matrix *A;
  int i;

  A = calloc(1,sizeof(plc_matrix));
  
  for(i=0;i<3;i++) {(*A)[0][i] = row0.c[i];}
  for(i=0;i<3;i++) {(*A)[1][i] = row1.c[i];}
  for(i=0;i<3;i++) {(*A)[2][i] = row2.c[i];}

  return A;
}

plc_matrix  *lu_decmp(plc_matrix *A,int index[3],double *d,bool *ok)

     /* Procedure uses Crout's method to perform an LU decomposition
	on A (transformed into an n_matrix for simplicity). The pivots mean
	that we decompose a row-wise permutation of A...

	Index should be an _allocated_ pointer to an array of 3 ints. 
	These will be set to the row permutation. 

	d = +/- 1.0, depending on whether the number of 
	transpositions is odd or even. */

     /* This is essentially copied directly from Numerical Recipes */

     /* Since all the following procedures use n_matrices internally, 
	the result is returned as an n_matrix to save some time. */

{
  int      i,imax={-1},j,k,idum;
  double   big,dum,sum,temp;
  double   vv[3];
  plc_matrix *cvt;
  
  double   TINY = {1.0e-20};
  
  /* First, we copy the input matrix. */
  
  cvt = calloc(1,sizeof(plc_matrix));
  for(i=0;i<3;i++) { for(j=0;j<3;j++) { (*cvt)[i][j] = (*A)[i][j]; } }
  
  /* Now we get to work. */
  
  for(i=0;i<3;i++) index[i] = i;         /* Initialize index. */
  
  *d = 1.0;                              /* No changes yet. */
  
  for(i=0;i<3;i++) {                     /* Get scaling for each row. */

    big = 0.0;

    for(j=0;j<3;j++) {

      if ((temp = fabs((*cvt)[i][j])) > big) big = temp;

    }

    if (big == 0.0) {

      *ok = false;
      return NULL;

    }

    vv[i] = 1.0/big;

  }

  for(j=0;j<3;j++) {                       /* Loop over columns */

    for(i=0;i<j;i++) {

      sum = (*cvt)[i][j];
      
      for(k=0;k<i;k++) {

	sum -= (*cvt)[i][k] * (*cvt)[k][j];

      }

      (*cvt)[i][j] = sum;

    }

    big = 0.0;                             /* Now we search for a pivot. */
      
    for(i=j;i<3;i++) {

      sum = (*cvt)[i][j];

      for(k=0;k<j;k++) {
	
	sum -= (*cvt)[i][k] * (*cvt)[k][j];

      }

      (*cvt)[i][j] = sum;

      if ((dum = vv[i] * fabs(sum)) >= big) {   

	big = dum;
	imax = i;

      }

    }

    if (j != imax) {                        /* Must we interchange rows? */

      for(k=0;k<3;k++) {                    /* Perform the interchange.. */

	dum = (*cvt)[imax][k];
	(*cvt)[imax][k] = (*cvt)[j][k];
	(*cvt)[j][k] = dum;

      }

      *d = -(*d);                           /* Change parity of d. */

      dum = vv[imax];                       /* Swap scale factors... */
      vv[imax] = vv[j];                     
      vv[j] = dum;

      idum = index[j];                       /* And swap index. */
      index[j] = index[imax];
      index[imax] = idum;

    }

    if ((*cvt)[j][j] == 0.0) (*cvt)[j][j] = TINY;

    if (j != 2) {

      dum = 1.0/((*cvt)[j][j]);

      for (i=j+1;i<3;i++) {

	(*cvt)[i][j] *= dum;

      }

    }

  }

  /* We are now done. Our mission is to return the data. */

  return cvt;

  }

plc_vector solve_linear(plc_matrix *LU,int index[3],plc_vector b)

     /* Procedure uses backsubstitution to solve the linear system

	                   Ax = b,

	given an LU-decomposition of A, the vector b, and the 
	permutation matrix returned by lu_decmp.

	It is worthwhile to break out the step of finding the 
	LU-decomposition-- this procedure is fast compared to the
	one above, and we may want to find many solutions of this
	equation with the same A, but different B.

     */

{
  plc_vector result;
  double B[3],scr_B[3];

  int i,j;
  double sum;

  /* We start by unscrambling the permutation. */

  scr_B[0] = b.c[0]; scr_B[1] = b.c[1]; scr_B[2] = b.c[2];

  for(i=0;i<3;i++) {

    B[i] = scr_B[index[i]];

  }

  for(i=0;i<3;i++) {   /* Now, we are ready to work.... */

    sum = B[i];

    for (j = 0; j <= i-1; j++) {

	sum -= (*LU)[i][j]*B[j];

    }

    B[i] = sum;

  }

  for(i=2;i>=0;i--) {

    sum = B[i];
    
    for (j=i+1;j<3;j++) { 

      sum -= (*LU)[i][j]*B[j];
     
    }

    B[i] = sum/((*LU)[i][i]);

  }
  
  /* Now prepare and return the result! */

  result = plc_build_vect(B[0],B[1],B[2]);

  return result;

}

plc_vector   plc_matrix_vector_multiply(plc_matrix *A,plc_vector x) 
{
  plc_vector build;
  double Ax[3] = {0,0,0};
  int i,j;

  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      Ax[i] += (*A)[i][j]*x.c[j];
    }
  }

  build = plc_build_vect(Ax[0],Ax[1],Ax[3]);
  return build;
}

plc_matrix   *plc_matrix_matrix_multiply(plc_matrix *A,plc_matrix *B)
{
  plc_matrix *AB;
  int i,j,k;

  AB = calloc(1,sizeof(plc_matrix));

  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      for(k=0;k<3;k++) {

	(*AB)[i][j] += (*A)[i][k] * (*B)[k][j];

      }
    }
  }

  return AB;
}

plc_matrix *plc_matrix_copy(plc_matrix *A)
{
  if (A == NULL) { return NULL; }
  
  plc_matrix *build;
  build = calloc(1,sizeof(plc_matrix));
  
  int i,j;
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      (*build)[i][j] = (*A)[i][j];
    }
  }

  return build;
}
  
