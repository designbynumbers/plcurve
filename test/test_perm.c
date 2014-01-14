/* 

   test_perm.c : Unit tests for the code in pd_perm.c.


*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#ifdef HAVE_STDIO_H
   #include<stdio.h>
#endif 

#ifdef HAVE_STRING_H
   #include<string.h>
#endif

#ifdef HAVE_STDINT_H
   #include<stdint.h>
#endif

#ifdef HAVE_STDLIB_H
   #include<stdlib.h>
#endif

#ifdef HAVE_GSL_GSL_PERMUTATION_H
   #include<gsl/gsl_permutation.h>
#endif

#ifdef HAVE_STDINT_H
  #include<stdint.h>
#endif

#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

#ifdef HAVE_ASSERT_H
  #include<assert.h>
#endif

#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>
#include<pd_isomorphisms.h>

#include"argtable2.h" /* We use a local copy of argtable */

int PD_VERBOSE=50;

bool test_nperm(pd_idx_t n, bool print)

{
  pd_perm_t *perm = pd_new_perm(&n);
  unsigned int i;

  printf("Testing perm group on %d elements.\n",n);
  printf("--------------------------------------\n");
  printf("nvals check ...                       ");

  unsigned int nvals = pd_nperms(perm);
  unsigned int expect;

  if (n == 0) { expect = 0; } 
  else {

    for(expect=1,i=1;i<=n;i++) { expect *= i; }
    
  }

  if (nvals != expect) { 

    printf("FAIL. %d instead of %d in pd_perm_nvals(%d) \n",nvals,expect,n);
    return false;

  } else {

    printf("pass (%d == %d).\n",nvals,expect);

  }

  printf("new perm is e test ...                ");
  if (!pd_perm_is_e(perm)) { 
    printf("FAIL. (new perm does not pass pd_perm_is_e).\n");
    return false;
  }
  printf("pass.\n");

  pd_perm_t *e;
  
  printf("copy of e is e test ...               ");
  e = pd_copy_perm(perm);
  if (!pd_perm_is_e(e)) { 
    printf("FAIL. (copy of e does not pass pd_perm_is_e).\n");
    return false;
  }
  printf("pass.\n");
  
  pd_perm_t *ee;
  printf("e * e == e test ...                   ");
  ee = pd_compose_perms(e,e);
   if (!pd_perm_is_e(ee)) { 
    printf("FAIL. (e*e does not pass pd_perm_is_e).\n");
    return false;
  }
  printf("pass.\n");
  pd_free_perm((void **)(&ee));

  pd_perm_t *epow;
  printf("e *= e == e test ...                  ");
  epow = pd_copy_perm(e);
  pd_stareq_perm(epow,e);
  if (!pd_perm_is_e(epow)) { 
    printf("FAIL. (e *= e does not pass pd_perm_is_e).\n");
    return false;
   }
  pd_free_perm((void **)(&epow));
  printf("pass.\n");

  printf("period of e is 1 test ...             ");
  unsigned int periode;
  periode = pd_perm_period(e);
  if (periode != 1) { 
    printf("FAIL. pd_perm_period(e) == %d != 1.\n",periode);
    return false;
  }
  printf("pass.\n");

  if (n >= 2) { /* Doesn't make sense to have cycles for a 1-element group */

    pd_perm_t *cycle2;
    cycle2 = pd_new_perm(&n);
    cycle2->map[0] = 1; cycle2->map[1] = 0;
    pd_regenerate_pcidx(cycle2);
    
    printf("(12) * e == (12) test ...             ");
    pd_perm_t *times_e;
    times_e = pd_compose_perms(cycle2,e);
    if (!pd_perms_eq(cycle2,times_e)) { 
      printf("FAIL. (12) * e != (12).\n");
      return false;
    }
    pd_free_perm((void **)(&times_e));
    printf("pass.\n");
    
    printf("(12) * (12) == e test ...             ");
    pd_perm_t *times_12;
    times_12 = pd_compose_perms(cycle2,cycle2);
    if (!pd_perm_is_e(times_12)) { 
      printf("FAIL. (12) * (12) != e.\n");
      return false;
    }
    pd_free_perm((void **)(&times_12));
    printf("pass.\n");
    
    printf("period of (12) is 2 test ...          ");
    unsigned int period12;
    period12 = pd_perm_period(cycle2);
    if (period12 != 2) { 
      printf("FAIL. pd_perm_period(12) == %d != 2.\n",period12);
      return false;
    }
    printf("pass.\n");
    pd_free_perm((void **)(&cycle2));
    
    printf("generating cycle (1 2 ... %2d) ...     ",n);
    pd_perm_t *cyclen;
    cyclen = pd_new_perm(&n);
    for(i=0;i<n;i++) {
      cyclen->map[i] = (i+1)%n; 
    }
    pd_regenerate_pcidx(cyclen);
    
    if(!pd_perm_ok(cyclen)) { 
      printf("FAIL. (1 2 ... %d) not ok.\n",n);
      return false;
    }
    printf("pass.\n");
    
    printf("period of (1 ... %d) is %d test ...     ",n,n);
    unsigned int periodn;
    periodn = pd_perm_period(cyclen);
    if (periodn != n) { 
      printf("FAIL. pd_perm_period(1 2 .. %d) == %d != %d.\n",n,periodn,n);
      return false;
    }
    printf("pass.\n");
    pd_free_perm((void **)(&cyclen));

  }
    
  printf("iteration test ...                    ");

  if (print) { printf("\n"); } 
  pd_perm_t **perm_buf;
  perm_buf = calloc(nvals,sizeof(pd_perm_t *));
  assert(perm_buf != NULL);

  for(i=0;i<nvals;i++,pd_increment_perm(perm)) {

    if (!pd_perm_ok(perm)) {

      pd_printf("FAIL. %PERM does not pass pd_perm_ok.\n",NULL,perm);
      return false;

    }

    perm_buf[i] = pd_copy_perm(perm);   /* Store this element of the group */

    /* Test regenerate_pcidx */

    pd_idx_t pc_idx = perm->pc_idx;
    perm->pc_idx = (pd_idx_t)(-1); /* Trash the pc index */
    pd_regenerate_pcidx(perm); 
    
    if (perm->pc_idx != pc_idx) {

      pd_printf("FAIL. Could not regenerate pc_idx for %PERM.\n",NULL,perm);
      return false;

    }

    if (print) {   /* print, if desired. */
      
      pd_printf("\t %d %PERM \n",NULL,i,perm);

    }

  }
  
  if (print) {
      printf("iteration test ...                    ");
  }
  printf("pass (all elts ok).\n");


  printf("checking return to identity ...       ");
  e = pd_new_perm(&n);
  
  if (!pd_perms_eq(e,perm)) {

    pd_printf("FAIL. %PERM != %PERM (e).\n",NULL,perm,e);
    return false;

  } 

  printf("pass.\n");

  /* Check that all perms unique */

  printf("checking all %6d perms unique ... ",nvals);

  if (!pd_perms_unique(nvals,perm_buf)) {

    printf("FAIL. Repeated permutation in list.\n");
    return false;

  }

  printf(" pass (all unique).\n");

  pd_free_perm((void **)(&perm));
  pd_free_perm((void **)(&e));
  pd_free_perm((void **)(&perm));

  for(i=0;i<nvals;i++) { pd_free_perm((void **)(&(perm_buf[i]))); }
  free(perm_buf);

  printf("--------------------------------------\n"
	 "%d-element perm group PASS. \n\n",n);

  return true;

}

bool test_perm() {

  if (!test_nperm(1,true)) { return false; }
  if (!test_nperm(2,true)) { return false; }
  if (!test_nperm(3,true)) { return false; }
  if (!test_nperm(4,true)) { return false; }
  if (!test_nperm(5,false)) { return false; }
  if (!test_nperm(6,false)) { return false; }
  if (!test_nperm(7,false)) { return false; }  
  /* if (!test_nperm(8,false)) { return false; } */
  /* if (!test_nperm(9,false)) { return false; } */

  return true;

}

int main() {

  printf("test_perm (%s)\n",PACKAGE_STRING);
  printf("---------------------------------------\n"
	 "Unit tests for pd_perm.c\n"
	 "=======================================\n");

  if (!pd_perm_pcdata_ok(true) || !test_perm()) {

    printf("=====================================\n");
    printf("test_perm:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================\n");
    printf("test_perm:  PASS.\n");
    exit(0);

  }

  return 0;

}
  
