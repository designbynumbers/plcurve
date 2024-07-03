/* 

   test_dihedral.c : Unit tests for the code in pd_dihedral.c.


*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#include<stdio.h>
#include<string.h>

#ifdef HAVE_STDINT_H
   #include<stdint.h>
#endif

#include<stdlib.h>

#ifdef HAVE_GSL_GSL_PERMUTATION_H
   #include<gsl/gsl_permutation.h>
#endif

#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

#include<assert.h>

#include<plcTopology.h>
#include<pd_multidx.h>
#include<pd_dihedral.h>

#ifdef HAVE_ARGTABLE2_H
   #include<argtable2.h> /* We use a local copy of argtable */
#endif



bool test_ndihedral(pd_idx_t n, bool print)

{
  pd_dihedral_t *d;
  pd_idx_t i;

  printf("Testing dihedral group on %d elements.\n",n);
  printf("--------------------------------------\n");
  printf("nvals check ...                       ");

  d=pd_new_dihedral(&n);

  unsigned int nvals = pd_ndihedrals(d);

  if (nvals != 2*n) { 

    fprintf(stderr,"FAIL. %d instead of %d in pd_ndihedrals(%d) \n",nvals,2*n,n);
    return false;

  } else {

    printf("pass (%d == %d).\n",nvals,2*n);

  }

  printf("iteration test ...                    ");

  if (print) { printf("\n"); }

  pd_dihedral_t **d_buf;
  d_buf = calloc(nvals,sizeof(pd_dihedral_t *));
  assert(d_buf != NULL);

  for(i=0;i<nvals;i++,pd_increment_dihedral(d)) {

    if (!pd_dihedral_ok(d)) {

      pd_printf("FAIL. %DIHEDRAL does not pass pd_dihedral_ok.\n",NULL,d);
      return false;

    }

    /* Store this element of the group */

    d_buf[i] = pd_copy_dihedral(d);

    /* print, if desired. */

    if (print) {
      
      pd_printf("\t %d %DIHEDRAL \n",NULL,i,d);

    }

  }
  
  if (print) {
      printf("iteration test ...                    ");
  }
  printf("pass (all elts ok).\n");


  printf("checking return to identity ...       ");

  pd_dihedral_t *e;
  e = pd_new_dihedral(&n);
  
  if (pd_dihedral_cmp(&e,&d) != 0) {

    pd_printf("FAIL. %DIHEDRAL != %DIHEDRAL (e).\n",NULL,d,e);
    return false;

  } 

  pd_printf("pass. (%DIHEDRAL == %DIHEDRAL (e)) \n",NULL,d,e);

  /* Check that all dihedrals unique */

  printf("checking all %3d dihedrals unique ... ",nvals);

  if (!pd_dihedrals_unique(nvals,d_buf)) { 
    
    return false; 

  }

  printf("pass (all unique).\n");

  pd_free_dihedral((void **)(&d));
  pd_free_dihedral((void **)(&e));
  pd_free_dihedral((void **)(&d));

  for(i=0;i<nvals;i++) { pd_free_dihedral((void **)(&(d_buf[i]))); }
  free(d_buf);

  printf("--------------------------------------\n"
	 "%d-element dihedral group PASS. \n\n",n);

  return true;

}

bool test_dihedral() {

  if (!test_ndihedral(1,true)) { return false; }
  if (!test_ndihedral(2,true)) { return false; }
  if (!test_ndihedral(3,true)) { return false; }
  if (!test_ndihedral(4,true)) { return false; }
  if (!test_ndihedral(5,false)) { return false; }
  if (!test_ndihedral(6,false)) { return false; }
  
  return true;

}

int main() {

  printf("test_dihedral (%s)\n",PACKAGE_STRING);
  printf("---------------------------------------\n"
	 "Unit tests for pd_dihedral.c\n"
	 "=======================================\n");

  if (!test_dihedral()) {

    printf("=====================================\n");
    printf("test_dihedral:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================\n");
    printf("test_dihedral:  PASS.\n");
    exit(0);

  }

  return 0;

}
  
