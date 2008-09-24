typedef double (plc_matrix)[3][3];
  
plc_matrix  *plc_build_matrix(plc_vector row0,plc_vector row1,plc_vector row2);
plc_matrix  *lu_decmp(plc_matrix *A,int index[3],double *d,bool *ok);
plc_vector  solve_linear(plc_matrix *LU,int index[3],plc_vector b);
