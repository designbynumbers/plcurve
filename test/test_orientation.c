/* 

   test_orientation.c : Unit tests for the code in pd_orientation.c.


*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#include<stdio.h>
#include<string.h>

#ifdef HAVE_STDINT_H
#include<stdint.h>
#endif

#include<math.h>
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
#include<pd_orientation.h>

#ifdef HAVE_ARGTABLE2_H
   #include<argtable2.h> /* We use a local copy of argtable */
#endif



bool test_norientation(pd_idx_t n, bool print)

{
  pd_orientation_t *d;
  pd_idx_t i;

  printf("Testing orientation group on %d elements.\n",n);
  printf("--------------------------------------\n");
  printf("nvals check ...                       ");

  d=pd_new_orientation(&n);

  unsigned int nvals = pd_norientations(d);

  if (nvals != (int)(pow(2,n))) { 

    fprintf(stderr,"FAIL. %d instead of %d in pd_norientations(%d) \n",nvals,2*n,n);
    return false;

  } else {

    printf("pass (%d == %d).\n",nvals,(int)(pow(2,n)));

  }

  printf("iteration test ...                    ");

  if (print) { printf("\n"); }

  pd_orientation_t **d_buf;
  d_buf = calloc(nvals,sizeof(pd_orientation_t *));
  assert(d_buf != NULL);

  for(i=0;i<nvals;i++,pd_increment_orientation(d)) {

    if (!pd_orientation_ok(d)) {

      pd_printf("FAIL. %ORIENTATION does not pass pd_orientation_ok.\n",NULL,d);
      return false;

    }

    /* Store this element of the group */

    d_buf[i] = pd_copy_orientation(d);

    /* print, if desired. */

    if (print) {
      
      pd_printf("\t %d %ORIENTATION \n",NULL,i,d);

    }

  }
  
  if (print) {
      printf("iteration test ...                    ");
  }
  printf("pass (all elts ok).\n");


  printf("checking return to identity ...       ");

  pd_orientation_t *e;
  e = pd_new_orientation(&n);
  
  if (pd_orientation_cmp(&e,&d) != 0) {

    pd_printf("FAIL. %ORIENTATION != %ORIENTATION (e).\n",NULL,d,e);
    return false;

  } 

  pd_printf("pass. (%ORIENTATION == %ORIENTATION (e)) \n",NULL,d,e);

  /* Check that all orientations unique */

  printf("checking all %3d orientations unique ... ",nvals);

  if (!pd_orientations_unique(nvals,d_buf)) { 
    
    return false; 

  }

  printf("pass (all unique).\n");

  pd_free_orientation((void **)(&d));
  pd_free_orientation((void **)(&e));
  pd_free_orientation((void **)(&d));

  for(i=0;i<nvals;i++) { pd_free_orientation((void **)(&(d_buf[i]))); }
  free(d_buf);

  printf("--------------------------------------\n"
	 "%d-element orientation group PASS. \n\n",n);

  return true;

}

bool test_orientation() {

  if (!test_norientation(1,true)) { return false; }
  
  return true;

}

int main() {

  printf("test_orientation (%s)\n",PACKAGE_STRING);
  printf("---------------------------------------\n"
	 "Unit tests for pd_orientation.c\n"
	 "=======================================\n");

  if (!test_orientation()) {

    printf("=====================================\n");
    printf("test_orientation:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================\n");
    printf("test_orientation:  PASS.\n");
    exit(0);

  }

  return 0;

}
  
