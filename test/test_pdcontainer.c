/* 

   test_pdcontainer.c : Unit tests for the code in pd_container.c.


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

#ifdef HAVE_STDBOOL_H
   #include<stdbool.h>
#endif

#include<ordie.h>
#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>

#include<pd_isomorphisms.h>

#include<pd_container.h>

#include"argtable2.h" /* We use a local copy of argtable */

int PD_VERBOSE=50;

bool init_and_read_tests(pd_contidx_t init,pd_contidx_t toadd,pd_container_t **rcont)

{
  pd_container_t *cont;
  *rcont = NULL;

  printf("new container (%6d elts)...",init);
  cont = pd_new_container(init);
  
  if (!pd_container_ok(cont)) { 

    printf("FAIL (not ok).\n");
    return false;

  } else {

    printf("           pass (container ok).\n");

  }

  printf("Adding %6d pd_code_t to container...",toadd);

  pd_contidx_t i;

  for(i=0;i<toadd;i++) {

    pd_code_t *pd;
    pd = pd_code_new(14); /* Create a 14-crossing pd code */
    pd->uid = i; /* Tag it with a unique id for testing */
    pd_addto_container(cont,pd);

  }

  if (!pd_container_ok(cont)) {

    printf("FAIL (container not ok).\n");
    return false;

  } else {

    printf("  pass (container ok).\n");
    
  }

  printf("Testing number of elements in container...");
  
  if (pd_container_nelts(cont) != toadd) { 

    printf("FAIL. (pd_container_nelts reports %d != %d elts expected).\n",
	   pd_container_nelts(cont),toadd);
    return false;

  } 

  printf("pass (%d elts == %d expected)\n",pd_container_nelts(cont),toadd);

  printf("Reading %6d elements back from cont...",toadd);

  for(i=0;i<toadd;i++) { 

    pd_code_t *pd;
    pd = (pd_code_t *)(pd_container_elt(cont,i));

    if (pd == NULL) {

      printf("FAIL. element %d NULL.\n",i);
      return false;

    }

    if(pd->uid != i) {

      printf("FAIL. element %d corrupted.\n",i);
      return false;

    }

  }

  printf("pass (all elements returned).\n");
  
  *rcont = cont;
  return true;

}

bool free_with_null_test(pd_container_t **cont)
{

  printf("Freeing container (with NULL) ... ");

  pd_free_container(cont,NULL);
  pd_free_container(cont,NULL);

  if (*cont != NULL) {

    printf("FAIL (didn't set cont to NULL).\n");
    return false;

  }

  printf("       pass\n\n");

  return true;

}

bool free_with_pd_code_eltfree_test(pd_container_t **cont)
{

  printf("Freeing container (with pd_code_free) ... ");

  pd_free_container(cont,pd_code_eltfree);
  pd_free_container(cont,pd_code_eltfree);

  if (*cont != NULL) {

    printf("FAIL (didn't set cont to NULL).\n");
    return false;

  }

  printf("       pass\n\n");

  return true;

}

typedef struct composite_struct {

    unsigned int id;
    int *buf;

} comp_t;

void comp_free(void **input) {

  comp_t **comp = (comp_t **)(input);

  if ((*comp)->buf != NULL) {

    free((*comp)->buf);
    (*comp)->buf = NULL;

  }
  
  free(*comp);
  *comp = NULL;

}

comp_t *new_comp_t(pd_contidx_t id)
{
  comp_t *comp;
  comp = calloc(1,sizeof(comp_t));
  comp->id = id;
  comp->buf = calloc(10,sizeof(int));

  return comp;
}

bool init_and_read_with_comp_tests(pd_contidx_t init,pd_contidx_t toadd,pd_container_t **rcont)

{
  pd_container_t *cont;
  *rcont = NULL;

  printf("new container (%6d elts)...",init);
  cont = pd_new_container(init);
  
  if (!pd_container_ok(cont)) { 

    printf("FAIL (not ok).\n");
    return false;

  } else {

    printf("            pass (container ok).\n");

  }

  printf("Adding %6d comp_t to container...",toadd);

  pd_contidx_t i;

  for(i=0;i<toadd;i++) {

    comp_t *comp;
    comp = new_comp_t(i);
    pd_addto_container(cont,comp);

  }

  if (!pd_container_ok(cont)) {

    printf("FAIL (container not ok).\n");
    return false;

  } else {

    printf("      pass (container ok).\n");
    
  }

  printf("Reading %6d elements back from cont...",toadd);

  for(i=0;i<toadd;i++) { 

    comp_t *comp;
    comp = pd_container_elt(cont,i);

    if (comp->id != i || comp->buf == NULL) {

      printf("FAIL. element %d corrupted.\n",i);
      return false;

    }

  }

  printf(" pass (all elements returned).\n");
  
  *rcont = cont;
  return true;

}

bool free_with_compfree_test(pd_container_t **cont)
{

  printf("Freeing container (with comp_free) ... ");

  pd_free_container(cont,comp_free);
  pd_free_container(cont,comp_free);

  if (*cont != NULL) {

    printf("FAIL (didn't set cont to NULL).\n");
    return false;

  }

  printf("   pass\n\n");
  return true;

}

bool push_and_pop_with_comp_tests(pd_contidx_t init,pd_contidx_t toadd)

{
  pd_container_t *cont;

  printf("new container (%6d elts)...",init);
  cont = pd_new_container(init);
  
  if (!pd_container_ok(cont)) { 

    printf("FAIL (not ok).\n");
    return false;

  } else {

    printf("            pass (container ok).\n");

  }

  printf("Adding %6d comp_t to container...",toadd);

  pd_contidx_t i;

  for(i=0;i<toadd;i++) {

    comp_t *comp;
    comp = new_comp_t(i);
    pd_addto_container(cont,comp);

  }

  if (!pd_container_ok(cont)) {

    printf("FAIL (container not ok).\n");
    return false;

  } else {

    printf("      pass (container ok).\n");
    
  }

  printf("Popping %6d elements back from cont...",toadd);

  for(i=0;i<toadd;i++) { 

    pd_contidx_t nelts = pd_container_nelts(cont);

    comp_t *comp;
    comp = pd_pop_container(cont);

    if (pd_container_nelts(cont) != nelts-1) { 

      printf("FAIL. (# elts did not go down after pop).\n");
      return false;

    }

    if (comp->id != (toadd-1)-i || comp->buf == NULL) {

      printf("FAIL. element %d corrupted.\n",i);
      return false;

    }

    comp_free((void **)(&comp)); /* Must get rid of this to avoid memory leak */

  }

  printf(" pass (all elements popped).\n");

  printf("testing 0 elements remain...");
  if (pd_container_nelts(cont) != 0) { 

    printf("FAIL. (%d elements remain)",pd_container_nelts(cont));

  }
  printf("              pass. (0 elts left)\n");

  printf("testing free with pd_free_fake()... ");

  pd_free_container(&cont,pd_free_fake);
  pd_free_container(&cont,pd_free_fake);

  printf("      pass\n");

  return true;

}
 
bool test_container() {

  pd_container_t *cont;

  printf("pd_container test suite\n"
	 "--------------------------------\n");

  if (!init_and_read_tests(100,756,&cont)) { return false; }
  if (!free_with_pd_code_eltfree_test(&cont)) { return false; }

  if (!init_and_read_tests(1000,756,&cont)) { return false; }
  if (!free_with_pd_code_eltfree_test(&cont)) { return false; }

  //if (!init_and_read_tests(100000,234567,&cont)) { return false; }
  //if (!free_with_pd_code_eltfree_test(&cont)) { return false; }

  printf("\n"
	 "Advanced freeing tests\n" 
	 "Should be run under valgrind's memcheck to be sure.\n\n");

  if (!init_and_read_with_comp_tests(100,756,&cont)) { return false; }
  if (!free_with_compfree_test(&cont)) { return false; }

  if (!init_and_read_with_comp_tests(1000,756,&cont)) { return false; }
  if (!free_with_compfree_test(&cont)) { return false; }

  //if (!init_and_read_with_comp_tests(10000,2367,&cont)) { return false; }
  //if (!free_with_compfree_test(&cont)) { return false; }

  printf("\n"
	 "Push and Pop tests.\n\n");

  if(!push_and_pop_with_comp_tests(100,756)) { return false; }
  if(!push_and_pop_with_comp_tests(1000,756)) { return false; }
  //  if(!push_and_pop_with_comp_tests(10000,75601)) { return false; }

  printf("-----------------------------------\n"
	 "pd_container test suite        PASS\n");

  return true;

}
 
int main() {

  printf("test_pdcontainer (%s)\n",PACKAGE_STRING);
  printf("---------------------------------------\n"
	 "Unit tests for pd_containers.c\n"
	 "=======================================\n");

  if (!test_container()) {

    printf("=====================================\n");
    printf("test_pdcontainer:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================\n");
    printf("test_pdcontainer:  PASS.\n");
    exit(0);

  }

  return 0;

}
	 
  
  
  
