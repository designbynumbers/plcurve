/* 

   pd_perm.c : The basic data type is pd_perm_t, 
   defined in pd_perm.h. This is a basic implementation 
   of the perm group.

*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#ifdef HAVE_SYS_STAT_H
  #include<sys/stat.h>
#endif

#ifdef HAVE_SYS_TYPES_H
  #include<sys/types.h>
#endif

#ifdef HAVE_ASSERT_H
  #include<assert.h>
#endif

#ifdef HAVE_STDINT_H
  #include<stdint.h>
#endif

#ifdef HAVE_STDIO_H
  #include<stdio.h>
#endif

#ifdef HAVE_STDLIB_H
  #include<stdlib.h>
#endif

#ifdef HAVE_STDINT_H
  #include<stdint.h>
#endif

#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

#ifdef HAVE_STRING_H
  #include<string.h>
#endif

#include"plcTopology.h"
#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_permdata.h> /* Precomputed permutation data */

pd_iterops_t perm_ops = {pd_new_perm,pd_free_perm,
			 pd_print_perm,pd_copy_perm,
			 pd_increment_perm,pd_nperms,
			 pd_perm_ok,pd_perm_cmp};

void *pd_new_perm(void *np)
{
  pd_perm_t *perm;
  pd_idx_t n = *(pd_idx_t *)(np);

  perm = calloc(1,sizeof(pd_perm_t));
  assert(perm != NULL);

  perm->n = n;
  perm->map = calloc(n,sizeof(pd_idx_t));
  assert(perm->map != NULL);

  pd_idx_t i;
  for(i=0;i<n;i++) { perm->map[i] = i; }

  perm->pc_idx = 0; /* Whether or not group precomputed, we clean this. */
    
  return perm;
}

void pd_free_perm(void **permp)
{
  pd_perm_t *perm = (pd_perm_t *)(*permp);

  if (perm != NULL) {

    if (perm->map != NULL) { free(perm->map); perm->map = NULL; }  
    perm->n = 0;
  
    free(perm);
    *permp = NULL;

  }
    
}

void pd_increment_perm(void *permP) 

/* Generates (in place) the next value of perm (in lex order). 
   To speed this up, we use precomputed data when we can. */

{
  pd_perm_t *perm = (pd_perm_t *)(permP);

  if (perm->n <= PD_MAX_PC_PERM) { 

    /* This is precomputed, so use pd_pc_perm. */

    perm->pc_idx = (perm->pc_idx + 1) % pd_n_pc_perms[perm->n];
    
    pd_idx_t i;

    for(i=0;i<perm->n;i++) {

      perm->map[i] = pd_pc_perm[perm->n][perm->pc_idx].p[i];

    }

  } else { 
    
    /* NOT precomputed, so update perm to the next (lex order) perm using Knuth
       algorithm:

       1. Find the largest index k such that a[k] < a[k + 1]. If
       no such index exists, the permutation is the last
       permutation.
       
       2. Find the largest index l such that a[k] < a[l]. Since
       k + 1 is such an index, l is well defined and
       satisfies k < l.  
       
       3. Swap a[k] with a[l].  
       
       4. Reverse the sequence from a[k + 1] up to and including
       the final element a[n-1]. */    

    perm->pc_idx = 0; /* Make sure this DOESN'T increment. */
    
    int  i,j,k,l,swap;
    bool lastperm = true;
    
    /*1*/ for(k=(perm->n)-2;k>=0;k--) { if(perm->map[k] < perm->map[k+1]) { lastperm=false; break;} } 
    
    if (lastperm) { 
      
      for(i=0;i<(perm->n);i++) { perm->map[i] = i; } /* Return to identity perm */
      return;

    }
    
    /*2*/ for(l=perm->n-1;l>k+1 && perm->map[k] >= perm->map[l];l--); 
    
    /*3*/ swap = perm->map[k]; perm->map[k] = perm->map[l]; perm->map[l] = swap;
    
    /*4*/ for(i=k+1,j=(perm->n)-1;i<j;i++,j--) { 
      
      swap = perm->map[i]; perm->map[i] = perm->map[j]; perm->map[j] = swap;
      
    }
    
  }
 
}

void   *pd_copy_perm(void *permP) 
/* Make a new-memory copy of perm */
{
  pd_perm_t *perm = (pd_perm_t *)(permP);
  pd_perm_t *perm_cpy;

  perm_cpy = pd_new_perm(&(perm->n));
  memcpy(perm_cpy->map,perm->map,perm->n*sizeof(pd_idx_t));
  perm_cpy->pc_idx = perm->pc_idx;

  return (void *)(perm_cpy);
}
  
unsigned int pd_nperms(void *permP) 

/* Count the number of unique values that perm can take (n!) */

{
  pd_perm_t *perm = (pd_perm_t *)(permP);

  if (perm->n <= PD_MAX_PC_PERM) {

    return pd_n_pc_perms[perm->n];

  }

  if (perm->n == 0) return 0;

  unsigned int np = 1;
  pd_idx_t i;

  for(i=1;i<=perm->n;i++) { np *= i; }

  return np;
}


bool pd_perm_ok(void *permP) 

/* Verifies that this is a meaningful element of the perm group and matches the */
/* appropriate precomputed data if applicable. */

{
  pd_perm_t *perm = (pd_perm_t *)(permP);
  pd_idx_t i,j;

  if (perm->n <= PD_MAX_PC_PERM) { /* Check perm->pc_idx legal. */

    if (perm->pc_idx > pd_nperms(perm)) {

      return pd_error(SRCLOC,"%PERM precomputed idx # %d > # %d-element precomp perms == %d",NULL,
		      perm,perm->pc_idx,pd_nperms(perm));

    }

    /* Now check that the map matches precomputed data with this index */

    for(i=0;i<perm->n;i++) {

      if (perm->map[i] != pd_pc_perm[perm->n][perm->pc_idx].p[i]) {

	return pd_error(SRCLOC,"%PERM != %PERM with index %d\n",NULL,
			perm,&(pd_pc_perm[perm->n][perm->pc_idx]),perm->pc_idx);

      }

    }

    /* Now return to standard tests. */

  }

  for(i=0;i<perm->n;i++) {

    if (perm->map[i] >= perm->n) { 

      return pd_error(SRCLOC,"%PERM contains illegal index %d at position %d",NULL,
		      perm->map[i],i);

    }

  }

  for(i=0;i<perm->n;i++) {

    for(j=i+1;j<perm->n;j++) {

      if (perm->map[i] == perm->map[j]) {

	return pd_error(SRCLOC,"%PERM contains repeated index (%d) at positions %d and %d\n",NULL,
			perm,perm->map[i],i,j);

      }

    }

  }

  /* Now we know that perm contains a valid permutation of 0..n-1 */

  return true;

}
  
bool pd_perms_eq(pd_perm_t *perm1,pd_perm_t *perm2)

{
  pd_idx_t i,n;

  if (perm1->n != perm2->n) { return false; }
  n = perm1->n;

  for(i=0;i<n;i++) {

    if (perm1->map[i] != perm2->map[i]) { return false; }

  }

  return true;

}

int pdint_perm_lex_cmp(pd_idx_t n,pd_idx_t *A,pd_idx_t *B)
/* Compare two buffers of n pd_idx_t variables in lex order. */
{
  pd_idx_t i;

  for(i=0;i<n;i++) {

    if (A[i] != B[i]) { return A[i] - B[i]; }

  }

  return 0;
}

int pd_pcperm_cmp_glob_nelts;

int pd_pcperm_cmp(const void *A,const void *B)
/* Compare two pd_pc_perm_t variables */
{
  pd_pc_perm_t *pcpA = (pd_pc_perm_t *)(A);
  pd_pc_perm_t *pcpB = (pd_pc_perm_t *)(B);

  return pdint_perm_lex_cmp(pd_pcperm_cmp_glob_nelts,pcpA->p,pcpB->p);
}

int pd_perm_cmp(const void *A,const void *B) 
/* Compare two permutations (in lex order) */
{
  pd_perm_t *pA = *(pd_perm_t **)(A);
  pd_perm_t *pB = *(pd_perm_t **)(B);

  pd_idx_t n;
  
  assert(pA->n == pB->n);
  n = pA->n;

  return pdint_perm_lex_cmp(n,pA->map,pB->map);

}

bool pd_perms_unique(unsigned int nperms,pd_perm_t **perm_buf)

{
  assert(perm_buf != NULL);

  /* Check to see if the permutations are unique */

  qsort(perm_buf,nperms,sizeof(pd_perm_t *),pd_perm_cmp);

  unsigned int i;

  for(i=0;i<nperms-1;i++) {

    if (pd_perms_eq(perm_buf[i],perm_buf[i+1])) {

      return pd_error(SRCLOC,"perm_buf contains %PERM == %PERM at positions %d, %d \n",NULL,
		      perm_buf[i],perm_buf[i+1],i,i+1);

    }

  }

  return true;

}

char *pd_print_perm(void *permP)

/* Generates a nice pretty-printed form for perms or returns NULL. */

{
  pd_perm_t *perm = (pd_perm_t *)(permP);
  
  char *buf;
  
  if (perm->n <= PD_MAX_PC_PERM && perm->pc_idx < pd_n_pc_perms[perm->n]) { 
    
    /* pp form exists in precomputed data and pc_idx is set correctly. */
    /* the pre-printed strings are all short, so a small buffer suffices */

    buf = malloc(32*sizeof(char)); assert(buf != NULL);
    strncpy (buf, pd_pc_perm_print[perm->n][perm->pc_idx], 32);
    
  } else {  

    /* We don't have a cycle decomposition printed, so we
       just print in permutation form. A larger buffer may
       be needed. */
    
    size_t bufsize = 256;
    buf = calloc(bufsize,sizeof(char)); assert(buf != NULL); 

    pd_idx_t i;
    size_t chars_ptd = 0;

    chars_ptd = snprintf(buf,bufsize-chars_ptd,"(");

    for(i=0;i<perm->n && chars_ptd < bufsize;i++) {

      chars_ptd += snprintf(buf + chars_ptd,bufsize-chars_ptd,"%d",perm->map[i]); 
      if (i != perm->n-1) { chars_ptd += snprintf(buf + chars_ptd,bufsize-chars_ptd," "); }

    }

    if (chars_ptd > bufsize) { 

      snprintf(buf,bufsize,"[PERM > %d CHARS]",(int)(bufsize));

    }

  }

  return buf;
  
}

void pd_regenerate_pcidx(pd_perm_t *perm) 
/* Recomputes the index of this perm's map in precomputed data.*/
{
  pd_pc_perm_t *pcpos,key;

  if (perm->n > PD_MAX_PC_PERM) { perm->pc_idx = 0; return; } /* This is correct. */

  pd_pcperm_cmp_glob_nelts = perm->n;

  memcpy(key.p,perm->map,perm->n*sizeof(pd_idx_t)); /* Make up a key */
  pcpos = bsearch(&key,pd_pc_perm[perm->n],pd_n_pc_perms[perm->n],sizeof(pd_pc_perm_t),pd_pcperm_cmp);

  if (pcpos == NULL) {

    pd_error(SRCLOC,"Could not find %PERM in corresponding precomputed data.\n",NULL,perm);
    exit(1);

  }

  perm->pc_idx = pcpos - pd_pc_perm[perm->n];

}

bool pd_perm_pcdata_ok(bool print) /* Self-tests on precomputed data. */

{
  if (print) {

    printf("pd_permdata.h test suite.\n"
	   "-----------------------------------------------------------\n");
    printf("Have precomputed data for perms of up to %d elts.\n",PD_MAX_PC_PERM);
    printf("Testing pd_perm_cmp order of precomputed data.... \n");
    
  }

  pd_idx_t nelts,perm;

  for(nelts=1;nelts<=PD_MAX_PC_PERM;nelts++) {

    if (print) { printf("\t %d element perm group ... ",nelts); }

    for(perm=0;perm<pd_n_pc_perms[nelts]-1;perm++) {

      if (pdint_perm_lex_cmp(nelts,(pd_pc_perm[nelts][perm].p),(pd_pc_perm[nelts][perm+1].p)) >= 0) {

	if (print) { printf("FAIL. %d >= %d (pdint_perm_lex_cmp) in pd_pc_perm[%d].\n",perm,perm+1,nelts); }
	return false;

      }

    }
    
    if (print) { printf("pass (all elts in order).\n"); }

  }

  if (print) {

    printf("-----------------------------------------------------------\n");
    printf("pd_permdata.h test suite PASS.\n\n");

  }
    
  return true;
}
  
bool  pd_perm_is_e(pd_perm_t *perm) /* Check for identity perm (for testing) */
{
  pd_idx_t i;

  for(i=0;i<perm->n;i++) { 

    if (perm->map[i] != i) { return false; }

  }

  return true;
}

pd_perm_t *pd_compose_perms(pd_perm_t *A,pd_perm_t *B) 
/* product automorphism (A * B)(x) = A(B(x)). */

{
  assert(A->n == B->n);
  pd_idx_t n = A->n;

  pd_perm_t *AstarB;
  AstarB = pd_new_perm(&n); 
  
  pd_idx_t i;
  for(i=0;i<n;i++) { 

    AstarB->map[i] = A->map[B->map[i]]; 
    /* i -> B->map[i] under B, then B->map[i] -> A->map[B->map[i]] under A */

  }

  pd_regenerate_pcidx(AstarB);

  return AstarB;
}

void pd_stareq_perm(pd_perm_t *A,pd_perm_t *B)  
/* A *= B (updates A in-place) */
{
  assert(A->n == B->n);
  pd_idx_t n = A->n;

  pd_idx_t *newmap;
  newmap = calloc(n,sizeof(pd_idx_t)); assert(newmap != NULL);

  pd_idx_t i;
  for(i=0;i<n;i++) { newmap[i] = A->map[B->map[i]]; }
  
  memcpy(A->map,newmap,n*sizeof(pd_idx_t));
  free(newmap);
}
  
unsigned int pd_perm_period(pd_perm_t *A) 
/* computes period of A in permutation group. */
{
  pd_perm_t *Apow;
  pd_idx_t  period,failsafe = 10000;

  for(period=1,Apow = pd_copy_perm(A);
      !pd_perm_is_e(Apow) && period<failsafe;
      period++,pd_stareq_perm(Apow,A));
 
  assert(period < failsafe);
  pd_free_perm((void **)(&Apow));
  
  return period;
}

pd_perm_t    *pd_inverse_perm(pd_perm_t *perm)
             /*inverts the permutation perm */
{
  pd_perm_t *iperm = pd_new_perm(&(perm->n));
  pd_idx_t i;

  for(i=0;i<perm->n;i++) {

    iperm->map[perm->map[i]] = i;

  }

  pd_regenerate_pcidx(iperm); /* Just in case we need it. */

  return iperm;

}
