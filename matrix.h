/*

  matrix.h : An internal (crude) set of linear algebra functions for plCurve. These are BLAS-independent 
  for portability reasons (and also because for matrices this small, there is little advantage).

*/

#ifndef PLC_MATRIX_H
#define PLC_MATRIX_H

#if (__cplusplus || c_plusplus)
extern "C" {
#endif

typedef double (plc_matrix)[3][3];
  
plc_matrix  *plc_build_matrix(plc_vector row0,plc_vector row1,plc_vector row2);
plc_matrix  *lu_decmp(plc_matrix *A,int index[3],double *d,bool *ok);
plc_vector   solve_linear(plc_matrix *LU,int index[3],plc_vector b);
plc_vector   plc_matrix_vector_multiply(plc_matrix *A,plc_vector x);
plc_matrix   *plc_matrix_matrix_multiply(plc_matrix *A,plc_matrix *B);
plc_matrix   *plc_matrix_copy(plc_matrix *A);

#if (__cplusplus || c_plusplus)
};
#endif
#endif
