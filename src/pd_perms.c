/* 

   pd_perms.c : This file currently unused. 

*/

struct pd_orperm_struct {

  pd_idx n;
  pd_idx_t *perm;        /* Both arrays of n elements. */
  pd_or_t  *or;

}


void pd_next_perm(pd_perm_t ) 

 

}



unsigned int pdint_perms_generate(pd_idx_t n,pd_idx_t *stuff,pd_idx_t **permstuff)

/* Generates all permutations of the n indices in stuff, 
   returning an array of pointers to n-element arrays of 
   pd_idx_t variables in permstuff representing the various
   permutations. 

   Returns the number of permutations generated. */

{
  unsigned int numperms = 1;
  pd_idx_t i,j,k;
  pd_idx_t *perm;

  for(i=2;i<=n;i++) { numperms *= i; }
  *permstuff = calloc(numperms,sizeof(pdx_idx_t *));
  assert(permstuff != NULL);

  perm = calloc(n,sizeof(pd_idx_t));
  assert(perm != NULL);
  for(i=0;i<n;i++) { perm[i] = i; }

  /* Now loop over the permutations. */

  for(i=0;i<numperms;i++,pdint_nextperm(n,perm)) {

    /* Make room for a new perm of stuff. */
    (*permstuff)[i] = calloc(n,sizeof(pd_idx_t));
    assert((*permstuff)[i] != NULL);

    /* Apply the current perm to stuff. */
    for(j=0;j<n;j++) { ((*permstuff)[i])[j] = stuff[perm[j]]; }

  } 

  return numperms;

}

unsigned int pdint_orientations_generate(pd_idx_t n,pd_or_t **ors)

/* Generates a list of all possible orientations for n objects */
/* and returns the list in ors. */  

{
  
}

void pdint_perms_free(pd_idx_t nperms,pd_idx_t ***permstuff)

{
  pd_idx_t i;

  if (*permstuff == NULL) { return; }

  for(i=0;i<nperms;i++) {

    if ((*permstuff)[i] != NULL) { free((*permstuff)[i]); }

  }

  free(*permstuff);
  *permstuff = NULL;
}
