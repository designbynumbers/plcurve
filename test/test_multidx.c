/* 

   test_multidx.c : Unit tests for the code in pd_multidx.c.


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

#include"argtable2.h" /* We use a local copy of argtable */

int PD_VERBOSE=50;


/* We are responsible for implementing the following interface for the iterator in 
   a multi-index: */

/* typedef struct pd_iterops_struct { */

/*   void *(*new)(void *);   /\* Given a pointer to initial data, create a new object *\/ */
/*   void  (*free)(void **); /\* Free everything allocated by new_obj, set pointer to NULL. *\/ */

/*   char *(*print)(void *); /\* Allocates a string representing obj. *\/ */
/*   void *(*copy)(void *);  /\* Allocates a new-memory copy of obj. *\/ */

/*   void  (*incr)(void *);  /\* Increments (in place) an object *\/ */
/*   unsigned int (*nvals)(void *); /\* After this many calls to incr, obj should return to initial state. *\/ */

/*   bool  (*ok)(void *);    /\* Performs internal consistency checks on obj. *\/ */
/*   int   (*cmp)(const void *,const void *); /\* Compares POINTERS to obj for searching and *\/ */
/*                                            /\* sorting arrays of POINTERS to obj types. *\/ */
/* } pd_iterops_t; */

void *pd_new_mtt(void *initdata);
void  pd_free_mtt(void **mttp);
char *pd_print_mtt(void *mttp);
void *pd_copy_mtt(void *mttp);
void  pd_incr_mtt(void *mttp);
unsigned int pd_nmtts(void *mttp);
bool  pd_mtt_ok(void *mttp);
int   pd_mtt_cmp(const void *mttppA,const void *mttppB);

pd_iterops_t mtt_ops = {pd_new_mtt,pd_free_mtt,pd_print_mtt,
			pd_copy_mtt,pd_incr_mtt,pd_nmtts,
			pd_mtt_ok,pd_mtt_cmp};

typedef struct multidx_test_struct {

  unsigned int state;
  unsigned int nvals;
  double       *test_buf;

} multidx_test_t;

/* We now give a very minimal iterator implementation. */

void *pd_new_mtt(void *np)
{
  multidx_test_t *mtt;
  pd_idx_t        n = *(pd_idx_t *)(np);
 
  mtt = calloc(1,sizeof(multidx_test_t));
  assert(mtt != NULL);
  
  mtt->nvals = n;
  mtt->state = 0;

  mtt->test_buf = calloc(n,sizeof(double)); assert(mtt->test_buf != NULL);
  return mtt;
}

void pd_free_mtt(void **mttp)
{
  multidx_test_t **mtt = (multidx_test_t **)(mttp);

  if (*mtt == NULL) { return; }
  if ((*mtt)->test_buf != NULL) { free((*mtt)->test_buf); }
  free(*mtt);
  *mtt = NULL;

}

char *pd_print_mtt(void *mttp)
{
  multidx_test_t *mtt = (multidx_test_t *)(mttp);
  char *ps;
  
  ps = calloc(32,sizeof(char)); assert(ps != NULL);
  sprintf(ps,"%.4d/%.4d",mtt->state,mtt->nvals);

  return ps;
}
  
void *pd_copy_mtt(void *mttp) {

  multidx_test_t *mtt = (multidx_test_t *)(mttp);
  multidx_test_t *new_mtt;

  new_mtt = pd_new_mtt(&(mtt->nvals));
  
  new_mtt->state = mtt->state;
  new_mtt->nvals = mtt->nvals;
  
  memcpy(new_mtt->test_buf,mtt->test_buf,mtt->nvals*sizeof(double));
  
  return (void *)(new_mtt);

}

void pd_incr_mtt(void *mttp)
{
  multidx_test_t *mtt = (multidx_test_t *)(mttp);
  mtt->state = (mtt->state + 1) % mtt->nvals;
}

unsigned int pd_nmtts(void *mttp)
{
  multidx_test_t *mtt = (multidx_test_t *)(mttp);
  return (unsigned int)(mtt->nvals);
}

bool pd_mtt_ok(void *mttp) 
{
  multidx_test_t *mtt = (multidx_test_t *)(mttp);
  if (mtt->state < mtt->nvals) { return true; }
  else {return false; }
}

int  pd_mtt_cmp(const void *mttpA,const void *mttpB) 
{
  multidx_test_t *mttA = *(multidx_test_t **)(mttpA);
  multidx_test_t *mttB = *(multidx_test_t **)(mttpB);

  return (int)(mttA->state) - (int)(mttB->state);
}

/* We now write a pretty general test based on the idea
   that we'll take an array of pointers to pd_idx_t as 
   a vector of initial data. */

bool test_nmultidx(pd_idx_t nobj,void **nvals,pd_iterops_t ops,bool print, char *typestring)

{
  pd_multidx_t *idx;
  pd_idx_t i;

  printf("(");
  for(i=0;i<nobj;i++) {
    printf("%d",*(pd_idx_t *)(nvals[i]));
    if (i != nobj-1) { printf("-"); }
  }
  printf(")");

  printf(" %s multidx test suite",typestring);

  printf("\n-----------------------------------\n");

  printf("nvals check ... ");
  
  idx = pd_new_multidx(nobj,(void **)(nvals),ops);
  unsigned int returned_nvals = pd_multidx_nvals(idx);

  printf("expect %d vals from this multidx.\n",returned_nvals);


  printf("iteration test.... ");

  if (print) {

    printf("\n");

  }

  pd_multidx_t **idx_buf;
  idx_buf = calloc(returned_nvals,sizeof(pd_multidx_t *));
  assert(idx_buf != NULL);
  
  for(i=0;i<returned_nvals;i++,pd_increment_multidx(idx)) {
    
    if (!pd_multidx_ok(idx)) {
      
      pd_printf("FAIL. %MULTIDX fails pd_multidx_ok.\n",NULL,idx);
      return false;

    }

    if (print) {

      pd_printf("\t %MULTIDX \n",NULL,idx);

    }

    idx_buf[i] = pd_copy_multidx(idx);

  }
  
  for(i=0;i<idx->nobj;i++) {

    if (idx->i[i] != 0) { 

      pd_printf("FAIL. %MULTIDX index vals did not return to all 0's\n"
		"after pd_multidx_nvals calls to pd_increment_multidx.\n",
		NULL,idx);
      return false;

    }

  }

  pd_printf("pass.\n\t(%MULTIDX == zeros) \n",NULL,idx);

  printf("checking all %4d elements unique ... ",returned_nvals);

  if (!pd_multidxs_unique(returned_nvals,idx_buf)) {

    printf("FAIL. (buffer of multidxs contains repeats).\n");
    return false;

  }

  printf(" pass.\n");

  pd_free_multidx(&idx);
  pd_free_multidx(&idx);

  for(i=0;i<returned_nvals;i++) {

    pd_free_multidx(&(idx_buf[i]));

  }

  free(idx_buf);

  printf("-----------------------------\n");

  printf("(");
  for(i=0;i<nobj;i++) {
    printf("%d",*(pd_idx_t *)(nvals[i]));
    if (i != nobj-1) { printf("-"); }
  }
  printf(")");

  printf(" %s multidx test suite",typestring);
  printf(" PASS\n\n");

  return true;

}

bool test_multidx() {

  void*     idata[10];
  pd_idx_t  nvals[10],i;

  for(i=0;i<10;i++) { idata[i] = &(nvals[i]); }

  nvals[0] = 3; 
  if (!test_nmultidx(1,idata,mtt_ops,true,"multidx_test_t")) { return false; }

  nvals[0] = 3; nvals[1] = 5;
  if (!test_nmultidx(2,idata,mtt_ops,true,"multidx_test_t")) { return false; }

  nvals[0] = 3; nvals[1] = 5; nvals[2] = 2; nvals[3] = 1; nvals[4] = 7;
  if (!test_nmultidx(5,idata,mtt_ops,false,"multidx_test_t")) { return false; }

  nvals[0] = 3; nvals[1] = 5; nvals[2] = 2; nvals[3] = 3; nvals[4] = 7;
  if (!test_nmultidx(5,idata,mtt_ops,false,"multidx_test_t")) { return false; }

  nvals[0] = 3;
  if (!test_nmultidx(1,idata,perm_ops,true,"perm_t")) { return false; }

  nvals[0] = 3; nvals[1] = 2;
  if (!test_nmultidx(2,idata,perm_ops,true,"perm_t")) { return false; }

  nvals[0] = 3; nvals[1] = 3; nvals[2] = 2; nvals[3] = 1; nvals[4] = 2;
  if (!test_nmultidx(5,idata,perm_ops,false,"perm_t")) { return false; }

  nvals[0] = 7; nvals[1] = 2;
  if (!test_nmultidx(2,idata,perm_ops,false,"perm_t")) { return false; }

  nvals[0] = 3;
  if (!test_nmultidx(1,idata,dihedral_ops,true,"dihedral_t")) { return false; }

  nvals[0] = 3; nvals[1] = 2;
  if (!test_nmultidx(2,idata,dihedral_ops,true,"dihedral_t")) { return false; }

  nvals[0] = 3; nvals[1] = 3; nvals[2] = 2; nvals[3] = 1; nvals[4] = 2;
  if (!test_nmultidx(5,idata,dihedral_ops,false,"dihedral_t")) { return false; }

  nvals[0] = 7; nvals[1] = 2;
  if (!test_nmultidx(2,idata,dihedral_ops,false,"dihedral_t")) { return false; }

  return true;
    
}

int main() {

  printf("test_multidx (%s)\n",PACKAGE_STRING);
  printf("---------------------------------------\n"
	 "Unit tests for pd_multidx.c\n"
	 "=======================================\n");

  if (!test_multidx()) {

    printf("=====================================\n");
    printf("test_multidx:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================\n");
    printf("test_multidx:  PASS.\n");
    exit(0);

  }

  return 0;

}
  
