/* 

   pd_isomorphisms.c : This code finds all the
   isomorphisms between any pair of pdcodes where
   no more than MAXISOMCOMPONENTS components have
   the same number of edges.

   (This number can be increased by regenerating
   pd_orperms.h using the Mathematica notebook
   /data/generate_orperms.nb.)

   The basic data type here is the 

   pd_isom_t 

   defined in pd_isomorphisms.h. 

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

#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

#ifdef HAVE_STRING_H
  #include<string.h>
#endif

#include<plcTopology.h>
#include<pd_container.h>

#include<pd_multidx.h>
#include<pd_dihedral.h>
#include<pd_cyclic.h>
#include<pd_perm.h>
  
#include<pd_isomorphisms.h>

bool pd_compgrp_ok(pd_code_t *pd,pd_compgrp_t *cgrp)

/* Checks that indices are meaningful with respect to the pdcode. As usual, */
/* prints an error message and terminates the program if PD_VERBOSE > 10.   */

{
  if (cgrp->ncomps > pd->ncomps) {

    return pd_error(SRCLOC,"cgrp ncomps %d > pd->ncomps %d in pd %PD",pd,cgrp->ncomps,pd->ncomps);

  } 

  if (cgrp->comp == NULL) {

    return pd_error(SRCLOC,"cgrp comp list not allocated\n",pd);
    
  }

  pd_idx_t comp;

  for(comp=0;comp<cgrp->ncomps;comp++) {

    if (cgrp->comp[comp] >= pd->ncomps) {
      
      return pd_error(SRCLOC,"cgrp->comp[%d] (%d) >= pd->ncomps %d",pd,comp,cgrp->ncomps,pd->ncomps);

    } 

  }

  pd_idx_t compB;

  for(comp=0;comp<cgrp->ncomps;comp++) {

    for(compB=comp+1;compB<cgrp->ncomps;compB++) {

      if (cgrp->comp[comp] == cgrp->comp[compB]) {

	return pd_error(SRCLOC,"cgrp entries %d and %d are equal in %COMPGRP",pd,comp,compB,cgrp);

      }

    }

  }

  return true;

}

pd_compgrp_t *pd_build_compgrps(pd_code_t *pdA,pd_code_t *pdB,pd_idx_t *ngrps) 

/* Build groups of components with the same number of edges, 
   using the fact that the components of pdA and pdB are sorted
   by number of edges. 

   Returns NULL if the number of components with a given number
   of edges is different in pdA and pdB, and terminates the program
   if the number of components in any compgrp is > MAXISOMCOMPONENTS. 

   Allocates child memory in the compgrp types as needed.
   Returns the number of component groups in *ngrps. */

{
  pd_compgrp_t *gps;
  pd_idx_t      comp,ncomps,gp;  

  /* First, make sure that this can work by checking that the vector of 
     numbers of edges per component matches between pdA and pdB */
  
  if (pdA->ncomps != pdB->ncomps) { *ngrps = 0; return NULL; } else { ncomps = pdA->ncomps; }
  
  for(comp=0;comp<ncomps;comp++) {

    if (pdA->comp[comp].nedges != pdB->comp[comp].nedges) { 
       
      *ngrps = 0;
      return NULL;
    
    } 

  }

  /* Now we need to count the number of groups, using the fact
     that the number of edges in each component are sorted. */

  *ngrps = 1;
  for(comp=1;comp<pdA->ncomps;comp++) { 
    if (pdA->comp[comp].nedges != pdA->comp[comp-1].nedges) { (*ngrps)++; }
  }

  gps = calloc(*ngrps,sizeof(pd_compgrp_t));
  assert(gps != NULL); 

  /* We now count the number of components in each group on pass 2. */

  gps[0].ncomps=1; gps[0].nedges = pdA->comp[0].nedges; /* Start by adding component 0 to group 0 */
  for(gp=0,comp=1;comp<pdA->ncomps;comp++) { 
    if (pdA->comp[comp].nedges != pdA->comp[comp-1].nedges) { gp++; gps[gp].nedges = pdA->comp[comp].nedges; }
    gps[gp].ncomps++;
  }

  /* Knowing the number of components in each group, we can allocate space for them. */

  for(gp=0;gp<(*ngrps);gp++) { 
    gps[gp].comp = calloc(gps[gp].ncomps,sizeof(pd_idx_t));
    assert(gps[gp].comp != NULL);
    gps[gp].ncomps = 0; /* This is now the number of components actually written in */
  }

  /* On pass 3 (through the components), we can actually assign them to groups. */

  gps[0].comp[0] = 0; gps[0].ncomps = 1;
  for(gp=0,comp=1;comp<pdA->ncomps;comp++) { 
    if (pdA->comp[comp].nedges != pdA->comp[comp-1].nedges) { gp++; }
   
    assert(gps[gp].nedges == pdA->comp[comp].nedges);
    gps[gp].comp[gps[gp].ncomps] = comp;
    gps[gp].ncomps++;
  }
  
  /* We now check that the total number of components is right. */
  
  pd_idx_t totalcomps = 0;
  for(gp=0;gp<(*ngrps);gp++) { totalcomps += gps[gp].ncomps; }
  assert(totalcomps == pdA->ncomps);

  /* We now check that every component has been assigned uniquely to a group. */

#define PD_UNSET_GP -1

  int *gp_assigned;
  gp_assigned = calloc(pdA->ncomps,sizeof(int));
  assert(gp_assigned != NULL);
  
  for(comp=0;comp<pdA->ncomps;comp++) { gp_assigned[comp] = PD_UNSET_GP; }

  for(gp=0;gp<(*ngrps);gp++) { 
    for(comp=0;comp<gps[gp].ncomps;comp++) {
      assert(gp_assigned[gps[gp].comp[comp]] == PD_UNSET_GP);
      assert(pdA->comp[gps[gp].comp[comp]].nedges == gps[gp].nedges);
      gp_assigned[gps[gp].comp[comp]] = gp;
    }
  }

  for(comp=0;comp<pdA->ncomps;comp++) { 
    assert(gp_assigned[comp] != PD_UNSET_GP);
  }

  free(gp_assigned);

  /* We have now grouped the components correctly. Return the list. */

  return gps;

}

void pd_free_compgrps(pd_compgrp_t *grps,pd_idx_t ngrps) 
{

  if (grps == NULL || ngrps == 0) { return; }
  pd_idx_t gp;

  for(gp=0;gp<ngrps;gp++) { 

    if (grps[gp].comp != NULL) { free(grps[gp].comp); grps[gp].comp = NULL; }
    grps[gp].ncomps = 0;
    grps[gp].nedges = 0;
    
  }

  free(grps);

}

/******************* comp_perms **************************/

pd_perm_t **pd_build_compperms(pd_idx_t ngrps, pd_compgrp_t *compgrps,
			       unsigned int *ncompperms)

{
  assert(ngrps > 0 && compgrps != NULL);
  assert(ncompperms != NULL);

  /* 0. Set up a multi-idx with permutations to iterate over. */
  /*    Compute total number of components as a side effect.  */

  pd_multidx_t *idx;
  pd_idx_t*     *grpsizes;
  pd_idx_t      grp;
  pd_idx_t      ncomps = 0;

  grpsizes = calloc(ngrps,sizeof(pd_idx_t *));
  assert(grpsizes != NULL);

  for(grp=0;grp<ngrps;grp++) {

    grpsizes[grp] = &(compgrps[grp].ncomps);
    ncomps += *(grpsizes[grp]);
  
  }

  idx = pd_new_multidx(ngrps,(void **)(grpsizes),perm_ops);

  /* 1. Allocate the output permutations. */
  
  pd_perm_t **comp_perms;
  unsigned int i;

  *ncompperms = pd_multidx_nvals(idx);
  comp_perms = calloc(*ncompperms,sizeof(pd_perm_t *)); assert(comp_perms != NULL);
  
  for(i=0;i<*ncompperms;i++) {

    comp_perms[i] = pd_new_perm(&ncomps);

  }

  /* 2. Loop over the multi-idx and actually generate the perms. */
  
  for(i=0;i<*ncompperms;i++,pd_increment_multidx(idx)) {

    for(grp=0;grp<idx->nobj;grp++) { 

      pd_perm_t    *this_perm = (pd_perm_t *)(idx->obj[grp]);
      pd_compgrp_t *this_grp = &(compgrps[grp]);
      pd_idx_t      comp;

      /* Apply this_perm to the group of comp #s in this_grp. */

      for(comp=0;comp<this_grp->ncomps;comp++) {

	comp_perms[i]->map[this_grp->comp[comp]] = this_grp->comp[this_perm->map[comp]];

      }

      /* Recompute the pc_idx (since we didn't build comp_perms[i] by iteration, this is unset */

      pd_regenerate_pcidx(comp_perms[i]);

    } /* We should have built an entire permutation of all ncomps components. */

  }

  pd_free_multidx(&idx);
  free(grpsizes);

  return comp_perms;

}

void pd_free_compperms(unsigned int ncomp_perms,pd_perm_t ***comp_permsP)

/* Free all memory associated with the buffer of pd_perm_t pointers *comp_permsP */

{
  pd_perm_t **comp_perms = *(comp_permsP);

  if (comp_perms == NULL || ncomp_perms == 0) { return; }

  unsigned int i;

  for(i=0;i<ncomp_perms;i++) {

    pd_free_perm((void **)(&(comp_perms[i]))); 

  }

  free(comp_perms);
  *comp_permsP = NULL;

}

bool pd_compperms_ok(unsigned int ncomp_perms,pd_perm_t **comp_perms)

/* Check the buffer of component-wise permutations. */

{
  unsigned int i;

  for(i=0;i<ncomp_perms;i++) {

    if (!pd_perm_ok(comp_perms[i])) {

      return pd_error(SRCLOC,"Component permutation %PERM fails pd_perm_ok.",NULL,comp_perms[i]);

    }

  }

  if (!pd_perms_unique(ncomp_perms,comp_perms)) {

    return pd_error(SRCLOC,"List of %d component permutations contains duplicates.\n",NULL,ncomp_perms);

  }

  return true;

}

/********************** edgemaps ***********************/

pd_edgemap_t  *pd_new_edgemap(pd_idx_t *nedges)
/* Allocate (cleared) memory for new edgemap */
{
  pd_edgemap_t *edgemap;

  edgemap = calloc(1,sizeof(pd_edgemap_t)); assert(edgemap != NULL);
  edgemap->perm = (pd_perm_t *) pd_new_perm(nedges);
  
  edgemap->or = calloc(*nedges,sizeof(pd_or_t));
  assert(edgemap->or != NULL);
  pd_idx_t edge;
  for(edge=0;edge<*nedges;edge++) { edgemap->or[edge] = PD_POS_ORIENTATION; }
    
  return edgemap;
}
  
void pd_free_edgemap(pd_edgemap_t **edgemapP) 
/* Free all memory associated with edgemap */
{
  pd_edgemap_t *edgemap = *(edgemapP);
  if (edgemap == NULL) { return; }

  pd_free_perm((void **)(&(edgemap->perm)));
  if (edgemap->or != NULL) { free(edgemap->or); edgemap->or = NULL; }

  free(edgemap);
  *edgemapP = NULL;
}

char  *pd_print_edgemap(pd_edgemap_t *edgemap) 
/* Returns a new-memory character string containing printed rep of edgemap */
{
  char *str,*strptr;
  size_t bufsize, printed;

  bufsize = 5*edgemap->perm->n + 25; 
  str = calloc(bufsize,sizeof(char)); assert(str != NULL);
  strptr = str;

  printed = snprintf(strptr,bufsize,"edgemap (");
  assert(printed < bufsize);
  bufsize -= printed;
  strptr += printed;

  pd_idx_t edge;
  
  for(edge=0;edge<edgemap->perm->n;edge++) {
    
    printed = snprintf(strptr,bufsize,"%c%d",pd_print_or(edgemap->or[edge]),edgemap->perm->map[edge]);
    assert(printed < bufsize);
    bufsize -= printed;
    strptr += printed;

    if (edge != edgemap->perm->n - 1) { 

      printed = snprintf(strptr,bufsize," "); 
      assert(printed < bufsize);
      bufsize -= printed;
      strptr += printed;

    }
    
  }

  printed = snprintf(strptr,bufsize,")");
  assert(printed < bufsize);
  bufsize -= printed;
  strptr += printed;

  return str;
}

void *pd_copy_edgemap(pd_edgemap_t *edgemap)
/* Make a new-memory copy of edgemap */
{
  pd_edgemap_t *newmap;

  newmap = calloc(1,sizeof(pd_edgemap_t));  assert(newmap != NULL);
  newmap->perm = pd_copy_perm(edgemap->perm);
  newmap->or = calloc(edgemap->perm->n,sizeof(pd_idx_t));
  assert(newmap->or != NULL);

  memcpy(newmap->or,edgemap->or,edgemap->perm->n*sizeof(pd_or_t));

  return newmap;
}


bool pd_edgemap_ok(pd_edgemap_t *edgemap)
/* Check whether the permutation is ok */
{
  if (!pd_perm_ok(edgemap->perm)) { return false; }
  pd_idx_t edge;

  for(edge=0;edge<edgemap->perm->n;edge++) {

    if (!(edgemap->or[edge] == PD_POS_ORIENTATION || edgemap->or[edge] == PD_NEG_ORIENTATION)) { 

      return pd_error(SRCLOC,"%EDGEMAP has illegal orientation for edge %d",NULL,edgemap,edge);

    }

  }

  return true;

}


void pd_free_edgemaps(unsigned int nedgemaps,pd_edgemap_t ***edgemapsP) 
/* Free a buffer of (pointers to) edgemaps */
{
  if (nedgemaps == 0 || *edgemapsP == NULL) { return; }

  pd_edgemap_t **edgemaps = *edgemapsP;
  unsigned int i;

  for(i=0;i<nedgemaps;i++) { pd_free_edgemap(&(edgemaps[i])); }
  free(edgemaps);

  *edgemapsP = NULL;
}

int  pd_edgemap_cmp(const void *A,const void *B) 

/* Compare permutations (lex order) and orientations */

{
  pd_edgemap_t *edgemapA = *(pd_edgemap_t **)(A);
  pd_edgemap_t *edgemapB = *(pd_edgemap_t **)(B);

  int permcmp = pd_perm_cmp(&(edgemapA->perm),&(edgemapB->perm));
  if (permcmp != 0) { return permcmp; }

  assert(edgemapA->perm->n == edgemapB->perm->n);
  pd_idx_t n = edgemapA->perm->n;

  pd_idx_t edge;
  int orcmp;

  for(edge=0;edge<n;edge++) { 

    orcmp = pd_or_cmp(&(edgemapA->or[edge]),&(edgemapB->or[edge]));
    if (orcmp != 0) { return orcmp; }

  }
    
  return 0;
}
  
bool pd_edgemaps_ok(unsigned int nedgemaps,pd_edgemap_t **edgemap_buf)
/* Check that the buffer of edgemaps contains unique, ok edgemaps. */
{
  assert(edgemap_buf != NULL);

  /* 0th. Check that the buffer contains valid edgemaps. */

  unsigned int i;
  
  for(i=0;i<nedgemaps;i++) {
    
    if (!pd_edgemap_ok(edgemap_buf[i])) {
      
      return pd_error(SRCLOC,"%EDGEMAP at position %d of edgemap_buf not ok",NULL,edgemap_buf[i],i);
      
    }

  }

  /* First, we copy the old buffer of pointers. */

  pd_edgemap_t **emb_copy;
  emb_copy = calloc(nedgemaps,sizeof(pd_edgemap_t *)); assert(emb_copy != NULL);
  memcpy(emb_copy,edgemap_buf,nedgemaps*sizeof(pd_edgemap_t *));

  /* Now we sort the copy */

  qsort(emb_copy,(size_t)(nedgemaps),sizeof(pd_edgemap_t *),pd_edgemap_cmp);

  /* Now make sure all the entries are unique */

  for(i=0;i<nedgemaps-1;i++) {

    if (pd_edgemap_cmp(&(emb_copy[i]),&(emb_copy[i+1])) == 0) {

      /* Now search for these pointers in the original edgemap_buffer. */

      unsigned int j,k;
      bool ifound = false, iplusfound = false;

      for(j=0;j<nedgemaps;j++) { 

	if (emb_copy[i] == edgemap_buf[j]) { ifound = true; break; }

      }

      for (k=0;k<nedgemaps;k++) { 

	if (emb_copy[i+1] == edgemap_buf[k]) { iplusfound = true; break; }

      }

      assert(ifound && iplusfound);

      return pd_error(SRCLOC,
		      "edgemap_buf contains\n"
		      "%EDGEMAP == \n"
		      "%EDGEMAP at positions %d, %d \n"
		      "(sorted) emb_copy buf positions are %d, %d\n"
		      ,NULL,
		      emb_copy[i],emb_copy[i+1],j,k,i,i+1);

    }

  }

  /* Housekeeping */

  free(emb_copy);

  return true;
  
}  

pd_edgemap_t **pd_build_edgemaps(pd_code_t *pdA,pd_code_t *pdB,pd_perm_t *comp_perm,unsigned int *nedgemaps)

/* Given a component permutation, we can build the set of
   edgemaps from components of pdA to components of pdB by
   iterating over a multi-index composed of dihedral
   groups. */

{
  assert(comp_perm != NULL);
  assert(nedgemaps != NULL);
  assert(pdA != NULL);
  assert(pdB != NULL);

  assert(pdA->ncomps == pdB->ncomps && pdB->ncomps == comp_perm->n);
  pd_idx_t ncomps = pdA->ncomps;

  assert(pdA->nedges == pdB->nedges);
  pd_idx_t nedges = pdA->nedges;

  /* 0. Set up a multi-idx with dihedral groups to iterate over. */

  pd_multidx_t *idx;
  pd_idx_t*     *compsizes;
  pd_idx_t      comp;
  
  compsizes = calloc(ncomps,sizeof(pd_idx_t *));
  assert(compsizes != NULL);

  for(comp=0;comp<ncomps;comp++) {

    compsizes[comp] = &(pdA->comp[comp].nedges);
  
  }

  idx = pd_new_multidx(ncomps,(void **)(compsizes),dihedral_ops);

  /* 1. Allocate the output edgemaps. */
  
  pd_edgemap_t **edgemaps;
  unsigned int i;

  *nedgemaps = pd_multidx_nvals(idx);
  edgemaps = calloc(*nedgemaps,sizeof(pd_edgemap_t *)); assert(edgemaps != NULL);
  
  for(i=0;i<*nedgemaps;i++) {

    edgemaps[i] = pd_new_edgemap(&nedges);

  }

  /* 2. Loop over the multi-idx and actually generate the edgemaps. */
  
  for(i=0;i<*nedgemaps;i++,pd_increment_multidx(idx)) {

    for(comp=0;comp<idx->nobj;comp++) { 

      pd_dihedral_t *this_dihedral = (pd_dihedral_t *)(idx->obj[comp]);
      pd_idx_t       edge;

      /* Map the edge #s in component comp to the dihedral-group images
	 of these edge numbers in the target component given by comp_perm. */

      assert(pdA->comp[comp].nedges == pdB->comp[comp_perm->map[comp]].nedges);
      pd_idx_t this_comp_edges = pdA->comp[comp].nedges;

      for(edge=0;edge<this_comp_edges;edge++) {
	
	pd_idx_t edgefrom = pdA->comp[comp].edge[edge];
	pd_idx_t edgeto   = pdB->comp[comp_perm->map[comp]].edge[this_dihedral->map[edge]];

	edgemaps[i]->perm->map[edgefrom] = edgeto;
	edgemaps[i]->or[edgefrom] = this_dihedral->or;

      }

    } /* We should have built an entire permutation of all nedges edges of the pd_code. */

    /* Recompute the pc_idx, now that we've got the whole permutation in place. */

    pd_regenerate_pcidx(edgemaps[i]->perm);

  }

  pd_free_multidx(&idx);
  free(compsizes);

  return edgemaps;

}


bool pd_edgemap_consistent(pd_code_t *pdA,pd_code_t *pdB,pd_edgemap_t *edgemap)

/* Checks to see that the targets of the edges are legal,
   unique, and take components to components in an orientation
   consistent way. */

{
  assert(pdA != NULL && pdB != NULL);
  assert(edgemap != NULL);

  if (edgemap->perm->n != pdA->nedges || edgemap->perm->n != pdB->nedges) {

    return pd_error(SRCLOC,"%EDGEMAP is not consistent with nedges (%d) in pdA or nedges (%d) in pdB.\n",
		    NULL,edgemap,pdA->nedges,pdB->nedges);

  }

  if (pdA->ncomps != pdB->ncomps) { 

    return pd_error(SRCLOC,"pdA->ncomps (%d) != pdB->ncomps (%d)",NULL,pdA->ncomps,pdB->ncomps);

  }

  pd_idx_t ncomps = pdA->ncomps;
  pd_idx_t comp;

  for(comp=0;comp<ncomps;comp++) {
    
    pd_idx_t Aedge = pdA->comp[comp].edge[0];
    pd_idx_t target_comp;

    pd_component_and_pos(pdB,edgemap->perm->map[Aedge],&target_comp,NULL);
    pd_or_t  target_or=edgemap->or[Aedge];

    assert(pdA->comp[comp].nedges == pdB->comp[comp].nedges);
    pd_idx_t compedges = pdA->comp[comp].nedges;

    pd_idx_t edge;
    for(edge=1;edge<compedges;edge++) {

      Aedge = pdA->comp[comp].edge[edge];

      pd_idx_t this_comp;
      pd_component_and_pos(pdB,edgemap->perm->map[Aedge],&this_comp,NULL);
      
      if (this_comp != target_comp) {

  	return pd_error(SRCLOC,"%EDGEMAP takes edges in comp %d of pdA to comps %d and %d of pdB",
  			pdA,edgemap,comp,target_comp,this_comp);

      }

      pd_or_t this_or = edgemap->or[Aedge];

      if (this_or != target_or) {

  	return pd_error(SRCLOC,"%EDGEMAP takes edges in comp %d of pdA to orientations %OR, %OR in pdB",
  			pdA,edgemap,comp,
			target_or,
			this_or);

      }
      
    }

  }
  
  return true;

}

pd_edgemap_t  *pd_compose_edgemaps(pd_edgemap_t *A,pd_edgemap_t *B)
/* When the edgemap maps a pd to itself, we can iterate. 
   Creates a new memory (A * B)(pd) = A(B(pd)) */
{
  assert(A->perm->n == B->perm->n);
  pd_idx_t n = A->perm->n;
  
  pd_edgemap_t *AstarB;
  pd_idx_t edge;

  AstarB = pd_new_edgemap(&n);
  
  for(edge=0;edge<n;edge++) { 

    /* First step:

       edge -> B->perm->map[edge] 
       while its orientation is changed by B->or[edge].
       
       Next step:

       B->perm->map[edge] -> A->map[B->perm->map[edge]] 
       while its orientation is changed by A->or[B->perm->map[edge]] */

    AstarB->perm->map[edge] = A->perm->map[B->perm->map[edge]];
    AstarB->or[edge] = pd_compose_or(B->or[edge],A->or[B->perm->map[edge]]);

  }

  return AstarB;

}

void           pd_stareq_edgemap(pd_edgemap_t *A,pd_edgemap_t *B)
/* Compose A with B in-place. */
{
  pd_edgemap_t *AstarB;
  AstarB = pd_compose_edgemaps(A,B);
  memcpy(A->perm->map,AstarB->perm->map,A->perm->n*sizeof(pd_idx_t));
  memcpy(A->or,AstarB->or,A->perm->n*sizeof(pd_or_t));
  pd_free_edgemap(&AstarB);
}


void pd_apply_edgemap(pd_code_t *pd, pd_edgemap_t *edgemap)
/* Apply the transformation in edgemap to the pd, changing references to edges in 
   component, face, and crossing data, reorienting edges as needed. */

{
  /* Zeroth, we need to check that this edgemap is really fully ok. 
     If not, we're going to do Very Bad Things below. */

  assert(pd_edgemap_consistent(pd,pd,edgemap));

  pd_edge_t *new_edgebuf = calloc(pd->nedges,sizeof(pd_edge_t));
  pd_idx_t   e;

  /* First, we do the actual (possibly reversed) copy of edges
     to the corresponding new places in the edge buffer. */

  for(e=0;e<pd->nedges;e++) { 
   
    pd_reorient_edge(pd,edgemap->perm->map[e],edgemap->or[e]);
    new_edgebuf[e] = pd->edge[edgemap->perm->map[e]];

  }

  for(e=0;e<pd->nedges;e++) {
    pd->edge[e] = new_edgebuf[e];
  }

  free(new_edgebuf);

  /* Next, we replace REFERENCES to the old edgenumbers in the
     crossing and face data with corresponding NEW edgenumbers where
     they occur. 

     Remember that we've already checked that components go to
     components in an orientation-consistent way, so we don't have to
     worry about orientation (as the edges occur in components) or 
     about swapping edge REFERENCES in components: we have swapped 
     the actual edges underneath.

     We DO have to swap orientations in face references if needed.*/

  pd_idx_t c,f;

  /* We need to be careful here. Suppose that our edgeperm is (1 2 0), 
     meaning that the new edge 0 is the old edge 1, the new edge 1 is the 
     old edge 2, and the new edge 2 is the old edge 0. 

     When we see a REFERENCE to edge 1, we should replace that number
     with 0 (the new index of edge 1), not with 2 (the old edge which 
     is currently in position 1). That is, we need to compute the INVERSE
     permutation of the one we've just applied:

  */
  pd_perm_t *iperm = pd_inverse_perm(edgemap->perm);

  for(c=0;c<pd->ncross;c++) {

    pd_idx_t p;
    for(p=0;p<4;p++) { 

      pd->cross[c].edge[p] = iperm->map[pd->cross[c].edge[p]];
      
    }

  }

  for(f=0;f<pd->nfaces;f++) {

    for(e=0;e<pd->face[f].nedges;e++) {

      pd->face[f].edge[e] = iperm->map[pd->face[f].edge[e]];
      pd->face[f].or[e] = pd_compose_or(edgemap->or[pd->face[f].edge[e]],pd->face[f].or[e]);

    }

  }

  pd_free_perm((void **)(&iperm));

  /* All this has put the crossings and faces out of order, 
     and may have destroyed the hash, too so we regenerate them. */

  pd_regenerate_crossings(pd);
  pd_regenerate_faces(pd);
  pd_regenerate_hash(pd);

}


  

/************************ crossmaps *****************************/

pd_crossmap_t *pd_new_crossmap(pd_idx_t *ncross) 
/* Allocate new crossmap */
{
  assert(ncross != NULL);
  pd_crossmap_t *crossmap;

  crossmap = calloc(1,sizeof(pd_crossmap_t)); assert(crossmap != NULL);
  crossmap->perm = pd_new_perm(ncross);
  crossmap->or   = PD_POS_ORIENTATION;

  return crossmap;
}

void pd_free_crossmap(pd_crossmap_t **crossmapP)
/* Free all memory associated with crossmap */
{
  pd_crossmap_t *crossmap = *crossmapP;
  
  if (crossmap == NULL) { return; }

  pd_free_perm((void **)(&crossmap->perm));
  free(crossmap);
  *crossmapP = NULL;
}

char           *pd_print_crossmap(pd_crossmap_t *crossmap)
/* Produce a printed representation of crossmap in new memory. */
{
  char *buf, *pstring;
  size_t bufsize, printed;

  pstring = pd_print_perm(crossmap->perm);
  
  bufsize = strlen(pstring) + 32;
  buf = calloc(bufsize,sizeof(char));
  printed = snprintf(buf,bufsize,"crossmap %c %s",pd_print_or(crossmap->or),pstring);
  assert(printed < bufsize);

  free(pstring);
  return buf;
}


void           *pd_copy_crossmap(pd_crossmap_t *crossmap)
/* A new-memory copy of crossmap, or pass through a NULL pointer. */
{
  if (crossmap == NULL) { return NULL; }

  pd_crossmap_t *newcrossmap;

  newcrossmap = calloc(1,sizeof(pd_crossmap_t)); assert(newcrossmap != NULL);
  newcrossmap->perm = pd_copy_perm(crossmap->perm);
  newcrossmap->or = crossmap->or;

  return newcrossmap;
}

bool pd_crossmap_ok(pd_crossmap_t *crossmap) 
{
  if (crossmap == NULL) { return true; }
  return pd_perm_ok(crossmap->perm);
}

int  pd_crossmap_cmp(const void *A,const void *B) 
/* Compares **pd_crossmap_t */
{
  pd_crossmap_t *crossmapA = *(pd_crossmap_t **)(A);
  pd_crossmap_t *crossmapB = *(pd_crossmap_t **)(B);

  int cmp = pd_or_cmp(&(crossmapA->or),&(crossmapB->or));
  if (cmp != 0) { return cmp; }

  cmp = pd_perm_cmp(&(crossmapA->perm),&(crossmapB->perm));
  return cmp;
}
  

pd_crossmap_t *pdint_build_oriented_crossmap(pd_code_t *pdA,pd_code_t *pdB,
					     pd_edgemap_t *edgemap,pd_or_t or)

/* Try to generate a crossing map _with the given orientation_. Return crossmap or NULL. */


{
  assert(pdA->ncross == pdB->ncross);
  pd_idx_t ncross = pdA->ncross;

  pd_crossmap_t *crossmap;
  crossmap = pd_new_crossmap(&ncross);

  pd_idx_t cross;
  pd_crossing_t mapcross;

  /* Now we can try to generate the crmap with this orientation. */

  crossmap->or = or;

  for(cross=0;cross<ncross;cross++) {

    /* Step 0. Build the crossing which this would map to. */

    pd_pos_t pos;
    
    for(pos=0;pos<4;pos++) {

      mapcross.edge[pos] = edgemap->perm->map[pdA->cross[cross].edge[pos]];

    }

    pd_canonorder_cross(&mapcross,or);

    /* Step 1. Search for that crossing in pdB's crossing list */

    pd_crossing_t *bptr;
    bptr = bsearch(&mapcross,pdB->cross,pdB->ncross,sizeof(pd_crossing_t),pd_cross_cmp);

    /* Step 2. If found, translate to index, continue, if not found, go on to negative orientation. */

    if (bptr == NULL) {

      pd_free_crossmap(&crossmap);
      return NULL;

    }

    crossmap->perm->map[cross] = (pd_idx_t)(bptr - pdB->cross);
    /* Pointer arithmetic reveals the position of the
       crossing in the pdB->cross array. */
  
  }

  /* We have now generated the map for the permutation, but not 
     identified it with a precomputed perm. We try: */

  pd_regenerate_pcidx(crossmap->perm);

  return crossmap;

}

pd_crossmap_t *pdint_build_oriented_signed_crossmap(pd_code_t *pdA,pd_code_t *pdB,
						    pd_edgemap_t *edgemap,pd_or_t or)

/* Try to generate a crossing map _with the given orientation_. Return crossmap or NULL. */


{
  assert(pdA->ncross == pdB->ncross);
  pd_idx_t ncross = pdA->ncross;

  pd_crossmap_t *crossmap;
  crossmap = pd_new_crossmap(&ncross);

  pd_idx_t cross;
  pd_crossing_t mapcross;

  /* Now we can try to generate the crmap with this orientation. */

  crossmap->or = or;

  for(cross=0;cross<ncross;cross++) {

    /* Step 0. Build the crossing which this would map to. */

    pd_pos_t pos;
    
    for(pos=0;pos<4;pos++) {

      mapcross.edge[pos] = edgemap->perm->map[pdA->cross[cross].edge[pos]];

    }

    pd_canonorder_cross(&mapcross,or);

    /* Step 1. Search for that crossing in pdB's crossing list */

    pd_crossing_t *bptr;
    bptr = bsearch(&mapcross,pdB->cross,pdB->ncross,sizeof(pd_crossing_t),pd_cross_cmp);

    /* Step 2. If found, translate to index, continue, if not found, return (and we'll try the other orientation). */

    if (bptr == NULL) {

      pd_free_crossmap(&crossmap);
      return NULL;

    }

    crossmap->perm->map[cross] = (pd_idx_t)(bptr - pdB->cross);
    /* Pointer arithmetic reveals the position of the
       crossing in the pdB->cross array. */

    if (pdA->cross[cross].sign != pdB->cross[crossmap->perm->map[cross]].sign) { 

      pd_free_crossmap(&crossmap);
      return NULL;

    }
  
  }

  /* We have now generated the map for the permutation, but not 
     identified it with a precomputed perm. We try: */

  pd_regenerate_pcidx(crossmap->perm);

  return crossmap;

}


pd_crossmap_t **pd_build_crossmaps(pd_code_t *pdA,pd_code_t *pdB,
				   pd_edgemap_t *emap,unsigned int *ncrmaps)

/* 
   This is where the rubber meets the road. IF POSSIBLE, use the
   edgemap to build a map from the crossings of pdA to the crossings
   of pdB. The crossings in pdA and pdB should be sorted, so we can
   use the gcc searching functions to look for matching crossings.

   If the map reverses orientation in the plane (globally), it will be
   reflected in the cyclic ordering of the edges at each crossing,
   which will switch from clockwise to counterclockwise.

   A buffer to crossing maps generated (if 1 or 2) and NULL
   otherwise. Returns number of crmaps generated in ncrmaps. If there
   are no crossings, we return a buffer consisting of a single NULL
   pointer.

*/

{
  pd_crossmap_t **cross_buf;
  cross_buf = calloc(2,sizeof(pd_crossmap_t *)); assert(cross_buf != NULL);
  *ncrmaps = 0;

  if (pdA->ncross == 0 && pdB->ncross == 0) { 

    *ncrmaps = 1;
    return cross_buf;

  }

  cross_buf[0] = pdint_build_oriented_crossmap(pdA,pdB,emap,PD_POS_ORIENTATION);
  if (cross_buf[0] != NULL) { (*ncrmaps)++; }
  cross_buf[*ncrmaps] = pdint_build_oriented_crossmap(pdA,pdB,emap,PD_NEG_ORIENTATION);
  if (cross_buf[*ncrmaps] != NULL) { (*ncrmaps)++; }

  return cross_buf;
}

pd_crossmap_t **pd_build_signed_crossmaps(pd_code_t *pdA,pd_code_t *pdB,
					  pd_edgemap_t *emap,unsigned int *ncrmaps)

/* 
   This is where the rubber meets the road. IF POSSIBLE, use the
   edgemap to build a map from the crossings of pdA to the crossings
   of pdB. The crossings in pdA and pdB should be sorted, so we can
   use the gcc searching functions to look for matching crossings.

   If the map reverses orientation in the plane (globally), it will be
   reflected in the cyclic ordering of the edges at each crossing,
   which will switch from clockwise to counterclockwise.

   A buffer to crossing maps generated (if 1 or 2) and NULL
   otherwise. Returns number of crmaps generated in ncrmaps. If there
   are no crossings, we return a buffer consisting of a single NULL
   pointer.

*/

{
  pd_crossmap_t **cross_buf;
  cross_buf = calloc(2,sizeof(pd_crossmap_t *)); assert(cross_buf != NULL);
  *ncrmaps = 0;

  if (pdA->ncross == 0 && pdB->ncross == 0) { 

    *ncrmaps = 1;
    return cross_buf;

  }

  cross_buf[0] = pdint_build_oriented_signed_crossmap(pdA,pdB,emap,PD_POS_ORIENTATION);
  if (cross_buf[0] != NULL) { (*ncrmaps)++; }
  cross_buf[*ncrmaps] = pdint_build_oriented_signed_crossmap(pdA,pdB,emap,PD_NEG_ORIENTATION);
  if (cross_buf[*ncrmaps] != NULL) { (*ncrmaps)++; }

  return cross_buf;
}

void pd_free_crossmaps(unsigned int ncrmaps,pd_crossmap_t ***crossmap_bufP)

/* Free all (well, none, one or both) crossmaps. */

{
  pd_crossmap_t **crossmap_buf = *crossmap_bufP;

  if (crossmap_buf == NULL) { return; }

  pd_idx_t i;

  for(i=0;i<ncrmaps;i++) {

    pd_free_crossmap(&(crossmap_buf[i]));

  }

  free(crossmap_buf);
  *crossmap_bufP = NULL;

}

pd_crossmap_t  *pd_compose_crossmaps(pd_crossmap_t *crossmapA,pd_crossmap_t *crossmapB)
/* When a crossmap maps pd->pd, we can iterate. Makes a new-memory (A*B)(pd) = A(B(pd)). */
{
  pd_crossmap_t *AstarB;
  AstarB = calloc(1,sizeof(pd_crossmap_t)); assert(AstarB != NULL); 

  AstarB->perm = pd_compose_perms(crossmapA->perm,crossmapB->perm);
  AstarB->or   = pd_compose_or(crossmapA->or,crossmapB->or);

  return AstarB;
}

void pd_stareq_crossmap(pd_crossmap_t *crossmapA,pd_crossmap_t *crossmapB)
/* Compose A with B in-place. */
{
  pd_stareq_perm(crossmapA->perm,crossmapB->perm);
  crossmapA->or = pd_compose_or(crossmapA->or,crossmapB->or);
}

/**************************** face maps *****************************/

pd_facemap_t *pd_new_facemap(pd_idx_t *nfaces) 
/* Allocate new facemap */
{
  assert(nfaces != NULL);
  pd_facemap_t *facemap;

  facemap = calloc(1,sizeof(pd_facemap_t)); assert(facemap != NULL);
  facemap->perm = pd_new_perm(nfaces);
  facemap->or   = PD_POS_ORIENTATION;

  return facemap;
}

void pd_free_facemap(pd_facemap_t **facemapP)
/* Free all memory associated with facemap */
{
  pd_facemap_t *facemap = *facemapP;
  
  if (facemap == NULL) { return; }

  pd_free_perm((void **)(&facemap->perm));
  free(facemap);
  *facemapP = NULL;
}

char           *pd_print_facemap(pd_facemap_t *facemap)
/* Produce a printed representation of facemap in new memory. */
{
  char *buf, *pstring;
  size_t bufsize, printed;

  pstring = pd_print_perm(facemap->perm);
  
  bufsize = strlen(pstring) + 32;
  buf = calloc(bufsize,sizeof(char));
  printed = snprintf(buf,bufsize,"facemap %c %s",pd_print_or(facemap->or),pstring);
  assert(printed < bufsize);

  free(pstring);
  return buf;
}


void           *pd_copy_facemap(pd_facemap_t *facemap)
/* A new-memory copy of facemap. */
{
  pd_facemap_t *newfacemap;

  newfacemap = calloc(1,sizeof(pd_facemap_t)); assert(newfacemap != NULL);
  newfacemap->perm = pd_copy_perm(facemap->perm);
  newfacemap->or = facemap->or;

  return newfacemap;
}


bool pd_facemap_ok(pd_facemap_t *facemap) 
{
  return pd_perm_ok(facemap->perm);
}

int  pd_facemap_cmp(const void *A,const void *B) 
/* Compares **pd_facemap_t */
{
  pd_facemap_t *facemapA = *(pd_facemap_t **)(A);
  pd_facemap_t *facemapB = *(pd_facemap_t **)(B);

  int cmp = pd_or_cmp(&(facemapA->or),&(facemapB->or));
  if (cmp != 0) { return cmp; }

  cmp = pd_perm_cmp(&(facemapA->perm),&(facemapB->perm));
  return cmp;
}
  
pd_facemap_t *pdint_build_oriented_facemap(pd_code_t *pdA,pd_code_t *pdB,
					   pd_edgemap_t *edgemap,pd_or_t or)

/* Try to generate a facemap _with the given orientation_. Return facemap or NULL. */

{
  assert(pdA->nfaces == pdB->nfaces);
  pd_idx_t nfaces = pdA->nfaces;

  pd_facemap_t *facemap;
  facemap = pd_new_facemap(&nfaces);

  pd_idx_t face;
  pd_face_t mapface;

  /* Now we can try to generate the facemap with this orientation. */

  facemap->or = or;

  for(face=0;face<nfaces;face++) {

    /* Step 0. Build the face which this would map to. */

    pd_idx_t edge;

    mapface.nedges = pdA->face[face].nedges;
    mapface.edge = calloc(mapface.nedges,sizeof(pd_idx_t));
    mapface.or = calloc(mapface.nedges,sizeof(pd_or_t));
    assert(mapface.edge != NULL && mapface.or != NULL);
    
    for(edge=0;edge<pdA->face[face].nedges;edge++) {

      mapface.edge[edge] = edgemap->perm->map[pdA->face[face].edge[edge]];

    }

    pd_canonorder_face(&mapface,or);

    /* Step 1. Search for that face in pdB's face list. This is known to be sorted */
    /* if this pd code passed pd_faces_ok, so it's ok to bsearch. */

    pd_face_t *bptr;
    bptr = bsearch(&mapface,pdB->face,pdB->nfaces,sizeof(pd_face_t),pd_face_cmp);

    /* Now dispose of the mapface memory. */

    free(mapface.edge); free(mapface.or);
    mapface.edge = NULL; mapface.or = NULL;
    mapface.nedges = 0;

    /* Step 2. If found, translate to index, continue, if not found, go on to negative orientation. */

    if (bptr == NULL) {

      pd_free_facemap(&facemap);
      return NULL;

    }

    facemap->perm->map[face] = (pd_idx_t)(bptr - pdB->face);
    /* Pointer arithmetic reveals the position of the
       face in the pdB->face array. */
  
  }

  /* We have now generated the map for the permutation, but not 
     identified it with a precomputed perm. We try: */

  pd_regenerate_pcidx(facemap->perm);

  return facemap;

}

pd_facemap_t **pd_build_facemaps(pd_code_t *pdA,pd_code_t *pdB,
				 pd_edgemap_t *emap,unsigned int *nfacemaps)

/* 
   This is where the rubber meets the road. IF POSSIBLE, use the edgemap to build 
   a map from the faces of pdA to the faces of pdB. The faces in pdA and pdB
   should be sorted, so we can use the gcc searching functions to look for matching 
   faces. 

   If the map reverses orientation in the plane (globally), it will be reflected in 
   the cyclic ordering of the edges in each face, which will reverse.

   A buffer to face maps generated (if 1 or 2) and NULL otherwise. Returns
   number of facemaps generated in nfacemaps.

*/

{
  pd_facemap_t **face_buf;
  face_buf = calloc(2,sizeof(pd_facemap_t *)); assert(face_buf != NULL);
  *nfacemaps = 0;

  if (pdA->ncross == 0 && pdB->ncross == 0) { 

    *nfacemaps = 1;
    face_buf[0] = pd_new_facemap(&(pdA->nfaces)); /* That is, 2 faces */
    face_buf[0]->or = emap->or[0];     /* Reverses <=> we reverse the edge */
    face_buf[0]->perm->map[0] = 0;     /* Identity permutation. */
    face_buf[0]->perm->map[1] = 1;
    return face_buf;

  }

  face_buf[0] = pdint_build_oriented_facemap(pdA,pdB,emap,PD_POS_ORIENTATION);
  if (face_buf[0] != NULL) { (*nfacemaps)++; }
  face_buf[*nfacemaps] = pdint_build_oriented_facemap(pdA,pdB,emap,PD_NEG_ORIENTATION);
  if (face_buf[*nfacemaps] != NULL) { (*nfacemaps)++; }

  return face_buf;
}

void pd_free_facemaps(unsigned int nfacemaps,pd_facemap_t ***facemap_bufP)

/* Free all (well, none, one or both) facemaps. */

{
  pd_facemap_t **facemap_buf = *facemap_bufP;

  if (facemap_buf == NULL) { return; }

  pd_idx_t i;

  for(i=0;i<nfacemaps;i++) {

    pd_free_facemap(&(facemap_buf[i]));

  }

  free(facemap_buf);
  *facemap_bufP = NULL;

}
    
pd_facemap_t  *pd_compose_facemaps(pd_facemap_t *facemapA,pd_facemap_t *facemapB)
/* When a facemap maps pd->pd, we can iterate. Makes a new-memory (A*B)(pd) = A(B(pd)). */
{
  pd_facemap_t *AstarB;
  AstarB = calloc(1,sizeof(pd_facemap_t)); assert(AstarB != NULL); 

  AstarB->perm = pd_compose_perms(facemapA->perm,facemapB->perm);
  AstarB->or   = pd_compose_or(facemapA->or,facemapB->or);

  return AstarB;
}

void pd_stareq_facemap(pd_facemap_t *facemapA,pd_facemap_t *facemapB)
/* Compose A with B in-place. */
{
  pd_stareq_perm(facemapA->perm,facemapB->perm);
  facemapA->or = pd_compose_or(facemapA->or,facemapB->or);
}

/************************* isomorphisms ***************************/
pd_iso_t **pdint_build_isos(pd_code_t *pdA,pd_code_t *pdB,
			    unsigned int *nisos,bool FAILEARLY);

pd_iso_t **pd_build_isos(pd_code_t *pdA,pd_code_t *pdB,unsigned int *nisos)
{
  return pdint_build_isos(pdA,pdB,nisos,false);
}

pd_iso_t **pdint_build_isos(pd_code_t *pdA,pd_code_t *pdB,
			    unsigned int *nisos,bool FAILEARLY)
{
  /* For the pdcodes to even have a chance of being
     isomorphic, the hash codes, # crossings, # edges, #
     faces, edge # vector (over components), and edge #
     vect (over faces) must all match. We check these
     first, hoping to fail out early (and fast). */

  bool abortflag = false;
  
  *nisos = 0;

  pd_idx_t i;  
  for(i=0;i<PD_HASHSIZE;i++) { if (pdA->hash[i] != pdB->hash[i]) { return NULL; } }

  pd_idx_t ncross,nedges,ncomps,nfaces;

  if (pdA->ncross != pdB->ncross) { return NULL; } else {ncross = pdA->ncross;}
  if (pdA->nedges != pdB->nedges) { return NULL; } else {nedges = pdA->nedges;}
  if (pdA->ncomps != pdB->ncomps) { return NULL; } else {ncomps = pdA->ncomps;}
  if (pdA->nfaces != pdB->nfaces) { return NULL; } else {nfaces = pdA->nfaces;}
  
  pd_idx_t comp,face;

  /* Remember that # of edges is a primary sort criterion
     for both comps and faces, and that these buffers are
     maintained in sorted order if pdA and pdB pass
     pd_ok. So it's safe to just compare the buffers
     element-for-element in their existing order. */

  for(comp=0;comp<ncomps;comp++) {
    if (pdA->comp[comp].nedges != pdB->comp[comp].nedges) return NULL; }
  for(face=0;face<nfaces;face++) {
    if (pdA->face[face].nedges != pdB->face[face].nedges) return NULL; }

  /* We could put in weirder topological comparisons here
     like "number of edges in faces adjacent to face i" or
     "number of edges in faces adjacent to component j" if
     we liked. There's a danger of losing isomorphisms
     here if those have bugs, so there's some risk to
     implementing this.  */

  pd_compgrp_t   *compgrp;
  pd_idx_t        ngrps;

  compgrp = pd_build_compgrps(pdA,pdB,&ngrps);

  if (compgrp == NULL) { /* No compgrps, no isomorphism. */

    return NULL;

  }

  assert(ngrps != 0);  /* If the buffer isn't null, but
			  there are no groups, something
			  is very wrong. */

  /* At this point, there may BE isomorphisms and we'll
     have to check all the possibilities before ruling
     anything out. So we start by allocating a container
     to hold the isomorphisms */

  pd_container_t *isos; 
  isos = pd_new_container((pd_contidx_t)(1000));
  /* A colossal overestimate, but only 8K of memory */

  /* The next step is to generate component permutations. */

  pd_perm_t **comp_perms;
  unsigned int ncomp_perms,comp_perm;

  comp_perms = pd_build_compperms(ngrps,compgrp,&ncomp_perms);
  assert(pd_compperms_ok(ncomp_perms,comp_perms));

  for(comp_perm=0;comp_perm < ncomp_perms && !abortflag;comp_perm++) { 

    pd_edgemap_t **edgemaps;
    unsigned int  nedgemaps,edgemap;

    edgemaps = pd_build_edgemaps(pdA,pdB,comp_perms[comp_perm],&nedgemaps);
    assert(pd_edgemaps_ok(nedgemaps,edgemaps));

    for(edgemap=0;edgemap < nedgemaps && !abortflag;edgemap++) { 

      pd_crossmap_t **crossmaps;
      pd_facemap_t  **facemaps;
      unsigned int   ncrossmaps,nfacemaps,map,nmaps;

      crossmaps = pd_build_crossmaps(pdA,pdB,edgemaps[edgemap],&ncrossmaps);
      facemaps  = pd_build_facemaps(pdA,pdB,edgemaps[edgemap],&nfacemaps);

      assert(ncrossmaps == nfacemaps);
      nmaps = ncrossmaps;

      for(map=0;map<nmaps && !abortflag;map++) { /* We have actually generated an isomorphism! */

	pd_iso_t *new_iso;
	new_iso = calloc(1,sizeof(pd_iso_t)); assert(new_iso != NULL); 

	new_iso->compperm = pd_copy_perm(comp_perms[comp_perm]);
	new_iso->edgemap  = pd_copy_edgemap(edgemaps[edgemap]);
	new_iso->crossmap = pd_copy_crossmap(crossmaps[map]);
	new_iso->facemap  = pd_copy_facemap(facemaps[map]);

	assert(pd_iso_consistent(pdA,pdB,new_iso)); /* Check ok-ness as we generate */
	pd_addto_container(isos,new_iso);

	if (FAILEARLY) { abortflag = true; } /* We only need to know if there is ONE iso */

      } /* End of crossmap/facemap loop */

      pd_free_crossmaps(ncrossmaps,&crossmaps);
      pd_free_facemaps(nfacemaps,&facemaps);

    } /* End of edgemap loop */

    pd_free_edgemaps(nedgemaps,&edgemaps);

  } /* End of main (compperm) loop */

  pd_free_compperms(ncomp_perms,&comp_perms);
  pd_free_compgrps(compgrp,ngrps);
  
  /* Now we need to move the data from the container into a simple array */

  pd_iso_t **isobuf;
  
  *nisos = pd_container_nelts(isos);
  
  if (*nisos != 0) { /* This is not guaranteed, since pdA and pdB may not BE isomorphic. */

    isobuf = calloc(*nisos,sizeof(pd_iso_t *)); assert(isobuf != NULL);

    unsigned int i;
    for(i=0;i<*nisos;i++) { isobuf[i] = (pd_iso_t *)(pd_pop_container(isos)); }
    /* Transfer the isos (and their associated memory) to isobuf. */ 

    assert(pd_isos_unique(*nisos,isobuf));

  } else {

    isobuf = NULL;

  }
  
  /* When we're done this loop, there should be nothing in the container. */
  assert(pd_container_nelts(isos) == 0);
  pd_free_container(&isos,pd_free_fake); /* pd_free_fake should never be called. */
  
  return isobuf;    
}


bool pd_isomorphic(pd_code_t *pdA,pd_code_t *pdB) 

/* This is basically just a wrapper for pd_build_isos. */

{
  pd_iso_t **iso_buf;
  unsigned int nisos;

  assert(pdA != NULL); assert(pdB != NULL);
  /* If one of these has been pd_free'd by mistake, die here. */

  iso_buf = pdint_build_isos(pdA,pdB,&nisos,true); /* Run in FAILEARLY configuration */

  assert((nisos == 0 && iso_buf == NULL) || (nisos != 0 && iso_buf != NULL));

  if (iso_buf != NULL) { 

    pd_free_isos(&nisos,&iso_buf);
    return true;

  } 

  return false;
  
}

/****************** pd_iso standard primitives ********************/

void pd_free_iso(pd_iso_t **isoP)
/* Free all memory associated with iso. */
{
  pd_iso_t *iso = *isoP;

  if (iso == NULL) { return; } /* Already freed */

  pd_free_perm((void **)(&(iso->compperm)));
  pd_free_edgemap(&(iso->edgemap));
  pd_free_crossmap(&(iso->crossmap));
  pd_free_facemap(&(iso->facemap));

  free(iso);
  *isoP = NULL;
}

pd_iso_t *pd_copy_iso(pd_iso_t *iso) 
/* Make a new-memory copy of iso. */
{
  pd_iso_t *iso_copy;

  iso_copy = calloc(1,sizeof(pd_iso_t)); assert(iso_copy != NULL);

  iso_copy->compperm  = pd_copy_perm(iso->compperm);
  iso_copy->edgemap   = pd_copy_edgemap(iso->edgemap);
  iso_copy->crossmap  = pd_copy_crossmap(iso->crossmap);
  iso_copy->facemap   = pd_copy_facemap(iso->facemap);

  return iso_copy;
}

char *pd_print_iso(pd_iso_t *iso)
/* Returns a new-memory string describing iso */
{
  char *compperm_print,*edgemap_print,*crossmap_print,*facemap_print;
  char *buf;
  size_t bufsize,printed;
  
  compperm_print = pd_print_perm(iso->compperm);
  edgemap_print = pd_print_edgemap(iso->edgemap);
  crossmap_print = pd_print_crossmap(iso->crossmap);
  facemap_print = pd_print_facemap(iso->facemap);

  bufsize = strlen(compperm_print) + strlen(edgemap_print) + strlen(crossmap_print) + strlen(facemap_print) + 256;
  buf = calloc(bufsize,sizeof(char)); assert(buf != NULL);
  
  pd_idx_t ncross,nedges,nfaces,ncomps;
  ncross = iso->crossmap->perm->n;
  nedges = iso->edgemap->perm->n;
  nfaces = iso->facemap->perm->n;
  ncomps = iso->compperm->n;

  printed = snprintf(buf,bufsize,
		     "iso cr %d e %d f %d cmps %d\n"
		     "\t compperm %s \n"
		     "\t %s \n"
		     "\t %s \n"
		     "\t %s \n",
		     ncross,nedges,nfaces,ncomps,
		     compperm_print,
		     edgemap_print,
		     crossmap_print,
		     facemap_print);
  
  assert(printed < bufsize);
  free(compperm_print); free(edgemap_print); free(crossmap_print); free(facemap_print);

  return buf;
}

bool pd_iso_is_e(pd_iso_t *iso) /* Check whether iso is identity (for testing) */
{
  pd_idx_t edge;

  if (!pd_perm_is_e(iso->compperm))       { return false; }
  if (!pd_perm_is_e(iso->edgemap->perm))  { return false; }
  if (!pd_perm_is_e(iso->crossmap->perm)) { return false; }
  if (!pd_perm_is_e(iso->facemap->perm))  { return false; }

  if (iso->facemap->or != PD_POS_ORIENTATION) { return false; }
  if (iso->crossmap->or != PD_POS_ORIENTATION) { return false; }
  
  for(edge=0;edge<iso->edgemap->perm->n;edge++) {

    if (iso->edgemap->or[edge] != PD_POS_ORIENTATION) { return false; }

  }

  return true;

}

pd_iso_t *pd_compose_isos(pd_iso_t *A,pd_iso_t *B) 
/* product automorphism (A * B)(x) = A(B(x)). */

{
  pd_iso_t *AstarB;
  AstarB = calloc(1,sizeof(pd_iso_t)); assert(AstarB != NULL); 

  AstarB->compperm = pd_compose_perms(A->compperm,B->compperm);
  AstarB->edgemap  = pd_compose_edgemaps(A->edgemap,B->edgemap);
  AstarB->crossmap = pd_compose_crossmaps(A->crossmap,B->crossmap);
  AstarB->facemap  = pd_compose_facemaps(A->facemap,B->facemap);

  return AstarB;
}

void pd_stareq_iso(pd_iso_t *A,pd_iso_t *B)  
/* A *= B (updates A in-place) */
{
  pd_stareq_perm(A->compperm,B->compperm);
  pd_stareq_edgemap(A->edgemap,B->edgemap);
  pd_stareq_crossmap(A->crossmap,B->crossmap);
  pd_stareq_facemap(A->facemap,B->facemap);
}
  
unsigned int pd_iso_period(pd_iso_t *A) 
/* computes period of A in isomorphism group. */
{
  pd_iso_t *Apow;
  pd_idx_t  period,failsafe = 10000;

  for(period=1,Apow = pd_copy_iso(A);
      !pd_iso_is_e(Apow) && period<failsafe;
      period++,pd_stareq_iso(Apow,A));
 
  assert(period < failsafe);
  pd_free_iso(&Apow);
  
  return period;
}

void  pd_free_isos(unsigned int *nisos,pd_iso_t ***isobufP)
/* Free a buffer of pointers to pd_iso_t. */
{
  pd_iso_t **isobuf = *isobufP;

  assert((*nisos == 0 && isobuf == NULL) || (*nisos != 0 && isobuf != NULL));
  if (isobuf == NULL && *nisos == 0) { return; }

  unsigned int i;

  for(i=0;i<*nisos;i++) {

    pd_free_iso(&(isobuf[i]));
    
  }

  free(isobuf);
  *isobufP = NULL;
  *nisos = 0;
}
  
int pd_iso_cmp(const void *A,const void *B)
/* Compare two isomorphisms; we compare the individual elements */
/* using their "cmp" primitives. */
{
  pd_iso_t *isoA = *(pd_iso_t **)(A);
  pd_iso_t *isoB = *(pd_iso_t **)(B);

  int cmp;

  cmp = pd_perm_cmp(&(isoA->compperm),&(isoB->compperm));
  if (cmp != 0) { return cmp; }

  cmp = pd_edgemap_cmp(&(isoA->edgemap),&(isoB->edgemap));
  if (cmp != 0) { return cmp; }

  cmp = pd_crossmap_cmp(&(isoA->crossmap),&(isoB->crossmap));
  if (cmp != 0) { return cmp; }
  
  cmp = pd_facemap_cmp(&(isoA->facemap),&(isoB->facemap));
  return cmp;
}

bool pd_iso_ok(pd_iso_t *iso) 
/* 
   This is just a wrapper for the "element_ok" primitives. 
*/

{
  if (!pd_perm_ok(iso->compperm))     { return false; }
  if (!pd_edgemap_ok(iso->edgemap))   { return false; }
  if (!pd_crossmap_ok(iso->crossmap)) { return false; }
  if (!pd_facemap_ok(iso->facemap))   { return false; }

  return true;
}

bool pdint_in_cyclic_order(pd_idx_t a,pd_idx_t b,pd_idx_t n,pd_or_t or)
/* Returns true if a and b are consective in the oriented
   (by or) cyclic order on n elements. */
{
  assert(pd_or_ok(or));

  if (or == PD_POS_ORIENTATION) { /* Should be climbing */

    if (b == a + 1) { return true; }
    if (b == 0 && a == n-1) { return true; }

    return false;

  } else { /* Should be descending */

    if (b == a - 1) { return true; }
    if (b == n - 1 && a == 0) { return true; }
    
    return false;

  }

}

bool pdint_buffers_cycliceq(pd_idx_t *bufA,pd_idx_t *bufB,pd_idx_t n,pd_or_t or)
/* Decide whether buffers A and B (of length n) are
   (cyclic) rotations of one another if or is PD_POS_ORIENTATION 
   or orientation-reversed rotations if or is PD_NEG_ORIENTATION. */
{
  pd_idx_t i,ofs;
  pd_idx_t *bcopy;

  bcopy = calloc(n,sizeof(pd_idx_t));

  if (or == PD_POS_ORIENTATION) { 

    for(i=0;i<n;i++) { bcopy[i] = bufB[i]; }

  } else {

    for(i=0;i<n;i++) { bcopy[i] = bufB[(n-1)-i]; }

  }

  /* Now we try all offsets. */

  for(ofs=0;ofs<n;ofs++) { /* Just try all the offsets. */

    bool ofs_ok = true;
      
    for(i=0;i<n;i++) {
      
      if (bufA[i] != bcopy[(i+ofs)%n]) { /* This offset won't work! */
	
	ofs_ok = false;
	break; 

      }

    }

    if (ofs_ok) { free(bcopy); return true; }

  } 

  /* We tried all the offsets with no luck. */

  free(bcopy);
  return false;

}

bool pd_iso_consistent(pd_code_t *pdA,pd_code_t *pdB,pd_iso_t *iso)
/*
  This is a heavyweight check which checks each primitive against 
  the others to make sure that all the data agrees (and makes sense
  with respect to the two pd codes in question).
  
  A fair amount of searching and sorting is involved, but an iso
  which passes this qualifies for "certificate" status.
*/
{
  if (!pd_ok(pdA)) {

    pd_error(SRCLOC,"pd_iso_consistent called with !ok pd code (pdA) %PD",pdA);
    exit(1);

  }
  
  if (!pd_ok(pdB)) {

    pd_error(SRCLOC,"pd_iso_consistent called with !ok pd code (pdB) %PD",pdB);
    exit(1);

  }

  if (!pd_iso_ok(iso)) { return false; } /* An iso which is not ok cannot be consistent */
  
  /* We start by assembling some information about the pd codes */

  pd_idx_t ncross,nedges,nfaces,ncomps;

  assert(pdA->ncross == pdB->ncross);
  ncross = pdA->ncross;

  if (ncross != 0) { 

    if (ncross != iso->crossmap->perm->n) { 

      return pd_error(SRCLOC,"# elements in iso->crossmap->perm->n (%d) != ncross (%d) for pd codes",NULL,
		      iso->crossmap->perm->n,ncross);

    }

  }

  assert(pdA->nedges == pdB->nedges);
  nedges = pdA->nedges;

  if (nedges != iso->edgemap->perm->n) {

    return pd_error(SRCLOC,"# elements in iso->edgemap->perm->n (%d) != nedges (%d) for pd codes",NULL,
		    iso->edgemap->perm->n,nedges);

  } 
 
  assert(pdA->nfaces == pdB->nfaces);
  nfaces = pdA->nfaces;

  if (nfaces != iso->facemap->perm->n) {

    return pd_error(SRCLOC,"# elements in iso->facemap->perm->n (%d) != nfaces (%d) for pd codes",NULL,
		    iso->facemap->perm->n,nfaces);

  } 

  assert(pdA->ncomps == pdB->ncomps);
  ncomps = pdA->ncomps;

  if (ncomps != iso->compperm->n) {

    return pd_error(SRCLOC,"# elements in iso->compperm->n (%d) != ncomps (%d) for pd codes",NULL,
		    iso->compperm->n,ncomps);

  }
 
  pd_idx_t cross,edge,face,comp;

  /* The 0th check is whether the component permutation
     actually takes components with the same number of
     edges to each other. */

  for(comp=0;comp<ncomps;comp++) {

    pd_idx_t before_nedges = pdA->comp[comp].nedges;
    pd_idx_t after_comp    = iso->compperm->map[comp]; 
    pd_idx_t after_nedges  = pdB->comp[after_comp].nedges;

    if (before_nedges != after_nedges) { 

      pd_printf("%PERM takes %d edge %COMP",pdA,
		iso->compperm,before_nedges,comp);

      pd_printf("-> %d edge %COMP\n",
		pdB,after_nedges,after_comp);

      return pd_error(SRCLOC,"error, since %d != %d\n",NULL,before_nedges,after_nedges);

    }

  }

  /* The first check is edgemap consistent with compperm. */

  for(edge=0;edge<nedges;edge++) {

    pd_idx_t before_edge,after_edge;
    pd_idx_t before_comp,after_comp;

    before_edge = edge; 
    pd_component_and_pos(pdA,before_edge,&before_comp,NULL);
    
    after_edge = iso->edgemap->perm->map[before_edge];
    pd_component_and_pos(pdB,after_edge,&after_comp,NULL);

    if (after_comp != iso->compperm->map[before_comp]) { 

      return pd_error(SRCLOC,"%EDGEMAP takes %d (comp %d) -> %d (comp %d), but compperm takes comp %d -> comp %d",
		      NULL,iso->edgemap,before_edge,before_comp,after_edge,after_comp,
		      before_comp,iso->compperm->map[before_comp]);

    }

  }

  /* We have now checked that components are being mapped
     to one another as sets.  We check that all the edges
     on a given component are mapped in cyclic order to
     the other component with orientation consistent with
     the orientation records in the edgemap. */

  pd_or_t *comp_or;
  comp_or = calloc(ncomps,sizeof(pd_or_t));
  
  for(comp=0;comp<ncomps;comp++) { 

    pd_idx_t start_edge = pdA->comp[comp].edge[0];
    pd_idx_t compedge;

    comp_or[comp] = iso->edgemap->or[start_edge]; 

    for(compedge=1;compedge<pdA->comp[comp].nedges;compedge++) {
      
      pd_idx_t this_edge = pdA->comp[comp].edge[compedge];

      if (iso->edgemap->or[this_edge] != comp_or[comp]) { 

	return pd_error(SRCLOC,
			"edge %d (edge %d of comp %d) edgemap or %OR does \n"
			"not match edge %d (edge 0 of comp %d) or %OR",NULL,
			this_edge,compedge,comp,iso->edgemap->or[this_edge],
			start_edge,comp,comp_or[comp]);

      }

    }

  }

  /* We now know that the edges in each component are all
     mapped positively (or negatively), and we have used that
     to define the array comp_or of component orientations. 

     Do these component-wise orientations actually match
     with the mapping in iso->edgemap->perm? */

  for(comp=0;comp<pdA->ncomps;comp++) {

    pd_idx_t compedge;

    for(compedge=0;compedge<pdA->comp[comp].nedges;compedge++) { 

      pd_idx_t edge1 = pdA->comp[comp].edge[compedge];
      pd_idx_t edge2 = pdA->comp[comp].edge[(compedge+1) % pdA->comp[comp].nedges];

      pd_idx_t after1 = iso->edgemap->perm->map[edge1];
      pd_idx_t after2 = iso->edgemap->perm->map[edge2];

      pd_idx_t comp1,pos1,comp2,pos2;

      pd_component_and_pos(pdB,after1,&comp1,&pos1);
      pd_component_and_pos(pdB,after2,&comp2,&pos2);

      /* We have already checked comp1 == comp2, so
	 there's something weird if this is no longer
	 true. */

      assert(comp1 == comp2);

      if (!pdint_in_cyclic_order(pos1,pos2,pdB->comp[comp1].nedges,comp_or[comp])) {

	return pd_error(SRCLOC,
			"along comp %d of pdA, edges %d and %d (global edges %d and %d) map to\n"
			"edges %d and %d of comp %d of pdB (global edges %d and %d) which are not\n"
			"in cyclic order %OR",NULL,
			comp,compedge,(compedge+1)%pdA->comp[comp].nedges,
			edge1,edge2,pos1,pos2,comp1,after1,after2,&(comp_or[comp]));

      }

    }

  }

  free(comp_or); /* We can now dispose of comp_or, and we do in order to avoid
		    forgetting to do it later on. */
  comp_or = NULL; 

  /* We have now checked that the edgemap is internally
     consistent and also that it's consistent with the
     component permutation. This is as much as we can
     determine about the edgemap.

     We now turn to the crossmap. Here we can check the
     construction: the (transformed) collection of edges
     (according to the edgemap) really is (cyclically)
     equivalent to the collection of edgemaps in the
     target crossing.

     We deliberately check cyclic equivalence exhaustively
     rather than rotating into canonical order because we
     want the check to depend on two different bodies of
     code. */
  
  
  for(cross=0;cross<ncross;cross++) { 

    pd_crossing_t after_edgemap_cross; /* Crossing after edgemap */
    pd_pos_t      pos;

    for(pos=0;pos<4;pos++) { 

      after_edgemap_cross.edge[pos] = iso->edgemap->perm->map[pdA->cross[cross].edge[pos]]; 

    }

    pd_crossing_t *after_crossmap_cross = &(pdB->cross[iso->crossmap->perm->map[cross]]);

    if (!pdint_buffers_cycliceq(&(after_edgemap_cross.edge[0]),
				&(after_crossmap_cross->edge[0]),
				4,iso->crossmap->or)) { 

      return pd_error(SRCLOC,
		      "%CROSSMAP takes %d -> %OR %d, but \n"
		      "%EDGEMAP takes %CROSSPTR -> %CROSSPTR != %CROSSPTR\n"
		      "with orientation %OR.\n",
		      NULL,iso->crossmap,cross,iso->crossmap->or,iso->crossmap->perm->map[cross],
		      iso->edgemap,&(pdA->cross[cross]),&after_edgemap_cross,
		      after_crossmap_cross);
      
    }

  }

  /* We have now checked that the crossing map is
     compatible with the edgemap.  We now turn to checking
     the facemap. The first thing to check, of course is 
     that we're mapping faces to faces with the same number 
     of edges. */

  for(face=0;face<nfaces;face++) {

    pd_idx_t before_nedges = pdA->face[face].nedges;
    pd_idx_t after_face    = iso->facemap->perm->map[face]; 
    pd_idx_t after_nedges  = pdB->face[after_face].nedges;

    if (before_nedges != after_nedges) { 

      pd_printf("%FACEMAP takes %d edge %FACE",pdA,
		iso->facemap,before_nedges,face);

      pd_printf("-> %d edge %FACE\n",
		pdB,after_nedges,after_face);

      return pd_error(SRCLOC,"error, since %d != %d\n",NULL,before_nedges,after_nedges);

    }

  }
		       
  /* The next thing to check is that the facemap is
     compatible with the edgemap as well. This is
     basically the same code, but with longer buffers. */

  for(face=0;face<nfaces;face++) { 

    pd_face_t after_edgemap_face; /* Faceing after edgemap */
    pd_idx_t  facepos;

    after_edgemap_face.nedges = pdA->face[face].nedges;
    after_edgemap_face.edge = calloc(after_edgemap_face.nedges,sizeof(pd_idx_t));
    after_edgemap_face.or = calloc(after_edgemap_face.nedges,sizeof(pd_or_t));
    assert(after_edgemap_face.edge != NULL && after_edgemap_face.or != NULL);

    for(facepos=0;facepos<pdA->face[face].nedges;facepos++) { 

      after_edgemap_face.edge[facepos] = iso->edgemap->perm->map[pdA->face[face].edge[facepos]]; 

    }
    
    pd_face_t *after_facemap_face = &(pdB->face[iso->facemap->perm->map[face]]);

    assert(after_facemap_face->nedges == after_edgemap_face.nedges);
    pd_idx_t nedges = after_facemap_face->nedges;

    if (!pdint_buffers_cycliceq(&(after_edgemap_face.edge[0]),
				&(after_facemap_face->edge[0]),
				nedges,iso->facemap->or)) { 

      return pd_error(SRCLOC,
		      "%FACEMAP takes %d -> %OR %d, but \n"
		      "%EDGEMAP takes %FACEPTR -> %FACEPTR != %FACEPTR\n"
		      "with orientation %OR.\n",
		      NULL,iso->facemap,face,iso->facemap->or,iso->facemap->perm->map[face],
		      iso->edgemap,&(pdA->face[face]),&after_edgemap_face,
		      after_facemap_face);
      
    }

    /* Now free the after_edgemap face */
    free(after_edgemap_face.edge); free(after_edgemap_face.or);
    after_edgemap_face.edge = NULL; after_edgemap_face.or = NULL;
    after_edgemap_face.nedges = 0;

  }

  /* If we've gotten to this point, the facemap and
     crossmap are consistent with the edgemap, which is in
     turn consistent with the compperm. There's nothing
     else to check: this looks like a legitimate
     isomorphism between pdA and pdB. */

  return true;

}

bool      pd_isos_unique(unsigned int nisos,pd_iso_t **iso_buf) 

/* First, we copy the old buffer of pointers. */

{
  pd_iso_t **isobuf_copy;
  isobuf_copy = calloc(nisos,sizeof(pd_iso_t *)); assert(isobuf_copy != NULL);
  memcpy(isobuf_copy,iso_buf,nisos*sizeof(pd_iso_t *));

  /* Now we sort the copy */

  qsort(isobuf_copy,(size_t)(nisos),sizeof(pd_iso_t *),pd_iso_cmp);

  /* Now make sure all the entries are unique */

  unsigned int i;

  for(i=0;i<nisos-1;i++) {

    if (pd_iso_cmp(&(isobuf_copy[i]),&(isobuf_copy[i+1])) == 0) {

      /* Now search for these pointers in the original iso_buffer. */

      unsigned int j,k;
      bool ifound = false, iplusfound = false;

      for(j=0;j<nisos;j++) { 

	if (isobuf_copy[i] == iso_buf[j]) { ifound = true; break; }

      }

      for (k=0;k<nisos;k++) { 

	if (isobuf_copy[i+1] == iso_buf[k]) { iplusfound = true; break; }

      }

      assert(ifound && iplusfound);

      return pd_error(SRCLOC,
		      "iso_buf contains\n"
		      "%ISO == \n"
		      "%ISO at positions %d, %d \n"
		      "(sorted) isobuf_copy buf positions are %d, %d\n"
		      ,NULL,
		      isobuf_copy[i],isobuf_copy[i+1],j,k,i,i+1);

    }

  }

  /* Housekeeping */

  free(isobuf_copy);

  return true;

}

bool pd_is_diagram_isotopy(pd_iso_t *A, pd_code_t *pdA, pd_code_t *pdB)
/* Checks whether a given pd_isomorphism qualifies as an diagram isotopy. */
{
  pd_idx_t i;

  /* Check whether component tags match */

  for(i=0;i<A->compperm->n;i++) { 

    if (pdA->comp[i].tag != pdB->comp[A->compperm->map[i]].tag) { 

      return false; 

    }

  }

  /* Check whether we're flipping plane or not. */

  assert(A->crossmap->or == A->facemap->or); 
  pd_or_t plane_or = A->crossmap->or;

  /* All the edge orientations should agree (and agree with plane_or). */

  for(i=0;i<A->edgemap->perm->n;i++) { 

    if (A->edgemap->or[i] != plane_or) { 

      return false; 

    }

  }

  /* Finally, all the crossing signs should match. */

  for(i=0;i<A->crossmap->perm->n;i++) { 

    if (pdA->cross[i].sign != pdB->cross[A->crossmap->perm->map[i]].sign) {

      return false;

    }

  }

  /* We can now return true. */

  return true;

}
  

bool pd_diagram_isotopy_ok(pd_iso_t *A, pd_code_t *pdA,pd_code_t *pdB)
/* Checks to see whether A passes pd_iso_ok and pd_is_diagram_isotopy. */
{

  if (!pd_iso_ok(A)) { 

    return pd_error(SRCLOC,"pd_iso_t %ISO doesn't pass pd_iso_ok, so it can't pass pd_diagram_isotopy_ok.\n",pdA,A);

  }

  if (!pd_is_diagram_isotopy(A,pdA,pdB)) { 

    pd_error(SRCLOC,"pd_iso_t %ISO isn't an diagram isotopy between %PD (pdA) and ",pdA,A);
    return pd_error(SRCLOC,"pdB (%PD) ",pdB);

  }

  return true;

}

pd_edgemap_t **pd_build_oriented_edgemaps(pd_code_t *pdA,pd_code_t *pdB,pd_perm_t *comp_perm,unsigned int *nedgemaps)

/* Given a component permutation, we can build the set of edgemaps
   which consistently preserve orientation from ALL components of pdA
   to ALL components of pdB (or consistently reversing orientation
   from ALL components of A to all components of pdB) by iterating
   over a multi-index composed of CYCLIC groups. */

{
  assert(comp_perm != NULL);
  assert(nedgemaps != NULL);
  assert(pdA != NULL);
  assert(pdB != NULL);

  assert(pdA->ncomps == pdB->ncomps && pdB->ncomps == comp_perm->n);
  pd_idx_t ncomps = pdA->ncomps;

  assert(pdA->nedges == pdB->nedges);
  pd_idx_t nedges = pdA->nedges;

  /* 0. Set up a multi-idx with cyclic groups to iterate over. */

  pd_multidx_t *idx;
  pd_idx_t*     *compsizes;
  pd_idx_t      comp;
  
  compsizes = calloc(ncomps,sizeof(pd_idx_t *));
  assert(compsizes != NULL);

  for(comp=0;comp<ncomps;comp++) {

    compsizes[comp] = &(pdA->comp[comp].nedges);
  
  }

  idx = pd_new_multidx(ncomps,(void **)(compsizes),cyclic_ops);

  /* 1. Allocate the output edgemaps. */
  
  pd_edgemap_t **edgemaps;
  unsigned int i;

  *nedgemaps = 2*pd_multidx_nvals(idx);

  edgemaps = calloc(*nedgemaps,sizeof(pd_edgemap_t *)); assert(edgemaps != NULL);
  
  for(i=0;i<*nedgemaps;i++) {

    edgemaps[i] = pd_new_edgemap(&nedges);

  }

  /* 2. Loop over the multi-idx and generate orientation PRESERVING edgemaps. */
  
  for(i=0;i<pd_multidx_nvals(idx);i++,pd_increment_multidx(idx)) {

    for(comp=0;comp<idx->nobj;comp++) { 

      pd_cyclic_t *this_cyclic = (pd_cyclic_t *)(idx->obj[comp]);
      pd_idx_t     edge;

      /* Map the edge #s in component comp to the cyclic-group images
	 of these edge numbers in the target component given by comp_perm. */

      assert(pdA->comp[comp].nedges == pdB->comp[comp_perm->map[comp]].nedges);
      pd_idx_t this_comp_edges = pdA->comp[comp].nedges;

      for(edge=0;edge<this_comp_edges;edge++) {
	
	pd_idx_t edgefrom = pdA->comp[comp].edge[edge];
	pd_idx_t edgeto   = pdB->comp[comp_perm->map[comp]].edge[this_cyclic->map[edge]];
	
	edgemaps[i]->perm->map[edgefrom] = edgeto;
	edgemaps[i]->or[edgefrom] = PD_POS_ORIENTATION;	
     
      } 
      
    }

    /* We should have built an entire permutation of all nedges edges of the pd_code. */
    /* Recompute the pc_idx, now that we've got the whole permutation in place. */

    pd_regenerate_pcidx(edgemaps[i]->perm);
    
  }

  /* 3. Loop over multidx again and generate all orientation REVERSING edgemaps */

  int j;
  
  for(j=0;j<pd_multidx_nvals(idx);j++,pd_increment_multidx(idx),i++) {

    for(comp=0;comp<idx->nobj;comp++) { 

      pd_cyclic_t *this_cyclic = (pd_cyclic_t *)(idx->obj[comp]);
      pd_idx_t     edge;

      /* Map the edge #s in component comp to the cyclic-group images
	 of these edge numbers in the target component given by comp_perm. */

      assert(pdA->comp[comp].nedges == pdB->comp[comp_perm->map[comp]].nedges);
      pd_idx_t this_comp_edges = pdA->comp[comp].nedges;

      for(edge=0;edge<this_comp_edges;edge++) {
	
	pd_idx_t edgefrom = pdA->comp[comp].edge[edge];
	pd_idx_t edgeto   = pdB->comp[comp_perm->map[comp]].edge[
	      (pdB->comp[comp_perm->map[comp]].nedges-1) - this_cyclic->map[edge]];
	/* Note that this is REVERSED */
	
	edgemaps[i]->perm->map[edgefrom] = edgeto;
	edgemaps[i]->or[edgefrom] = PD_NEG_ORIENTATION;	
     
      } 
      
    }

    /* We should have built an entire permutation of all nedges edges of the pd_code. */
    /* Recompute the pc_idx, now that we've got the whole permutation in place. */

    pd_regenerate_pcidx(edgemaps[i]->perm);
    
  }

  /* 4. Housekeeping */
  
  pd_free_multidx(&idx);
  free(compsizes);

  return edgemaps;

}

pd_iso_t **pd_build_diagram_isotopies(pd_code_t *pdA,pd_code_t *pdB,unsigned int *nisos)
/* Returns a buffer of 0, 1, or 2 diagram isotopies between pdA and pdB. */
{
  /* For the pdcodes to even have a chance of being
     isomorphic, the hash codes, # crossings, # edges, #
     faces, edge # vector (over components), and edge #
     vect (over faces) must all match. We check these
     first, hoping to fail out early (and fast). */

  *nisos = 0;

  pd_idx_t i;  
  for(i=0;i<PD_HASHSIZE;i++) { if (pdA->hash[i] != pdB->hash[i]) { return NULL; } }

  pd_idx_t ncross,nedges,ncomps,nfaces;

  if (pdA->ncross != pdB->ncross) { return NULL; } else {ncross = pdA->ncross;}
  if (pdA->nedges != pdB->nedges) { return NULL; } else {nedges = pdA->nedges;}
  if (pdA->ncomps != pdB->ncomps) { return NULL; } else {ncomps = pdA->ncomps;}
  if (pdA->nfaces != pdB->nfaces) { return NULL; } else {nfaces = pdA->nfaces;}
  
  pd_idx_t comp,face;

  /* Remember that # of edges is a primary sort criterion
     for both comps and faces, and that these buffers are
     maintained in sorted order if pdA and pdB pass
     pd_ok. So it's safe to just compare the buffers
     element-for-element in their existing order. */

  for(comp=0;comp<ncomps;comp++) { if (pdA->comp[comp].nedges != pdB->comp[comp].nedges) return NULL; }
  for(face=0;face<nfaces;face++) { if (pdA->face[face].nedges != pdB->face[face].nedges) return NULL; }

  /* We could put in weirder topological comparisons here
     like "number of edges in faces adjacent to face i" or
     "number of edges in faces adjacent to component j" if
     we liked. There's a danger of losing isomorphisms
     here if those have bugs, so there's some risk to
     implementing this.  */

  /* At this point, there may BE isomorphisms and we'll
     have to check all the possibilities before ruling
     anything out. So we start by allocating a container
     to hold the isomorphisms */

  pd_container_t *isos; 
  isos = pd_new_container((pd_contidx_t)(1000)); /* A colossal overestimate, but only 8K of memory */

  /* The next step is to generate component permutations. */
  /* Since we are trying to generate diagram isotopies here, 
     there are basically 0 or 2 possibilities. We're going to 
     build them in-place using a slow (n^2) algorithm in the 
     number of components. */

  pd_perm_t *comp_perm;
  pd_idx_t j;  
  comp_perm = pd_new_perm(&pdA->ncomps);

  for(i=0;i<pdA->ncomps;i++) { 

    /* We now search for the component in B with a matching tag. */

    bool tag_found = false;

    for(j=0;j<pdB->ncomps && !tag_found;j++) { 

      if (pdB->comp[j].tag == pdA->comp[i].tag) {  /* We send this comp i of A to comp j of B. */

	comp_perm->map[i] = j;
	tag_found = true;

      }

    }

    if (!tag_found) { /* Something weird has happened. */

      pd_error(SRCLOC,"couldn't find tag from %COMP of pdA = %PD in",pdA,i);
      pd_error(SRCLOC,"pdB = %PD\n",pdB);
      exit(1);

    }

  }

  /* We SHOULD have assembled a permutation of components. */

  if (!pd_perm_ok(comp_perm)) { 

    pd_error(SRCLOC,"matching tags in pd_diagram_isotopic yielded component permutation %PERM, which is not ok\n",pdA,comp_perm);
    exit(1);

  }

  /* Now, it's entirely possible that the two pd codes have the same number of components,
     and even a valid tagset, but that the permutation is impossible for some easy reason,
     like the number of edges between the different components doesn't match. We check for 
     that now, noting that we'll have to fail out cleanly (without losing memory) if things
     don't pan out. */

  for(i=0;i<pdA->ncomps;i++) { 

    if (pdA->comp[i].nedges != pdB->comp[comp_perm->map[i]].nedges) { 

      pd_free_perm((void **)(&comp_perm));
      *nisos = 0;
      return NULL;

    }

  }

  /* If we've passed this point, we can go ahead and try to generate edgemaps. */

  pd_edgemap_t **edgemaps;
  unsigned int  nedgemaps,edgemap;

  edgemaps = pd_build_oriented_edgemaps(pdA,pdB,comp_perm,&nedgemaps);
  assert(pd_edgemaps_ok(nedgemaps,edgemaps));

  for(edgemap=0;edgemap < nedgemaps;edgemap++) { 

    pd_crossmap_t **crossmaps;
    pd_facemap_t  **facemaps;
    unsigned int   ncrossmaps,nfacemaps,map,nmaps;
    
    crossmaps = pd_build_signed_crossmaps(pdA,pdB,edgemaps[edgemap],&ncrossmaps);

    if (ncrossmaps > 0) { /* We may not actually be able to construct any crossing maps
			     which map crossing signs correctly. In this case, there's 
			     no need to try to extend things to face maps. */

      facemaps  = pd_build_facemaps(pdA,pdB,edgemaps[edgemap],&nfacemaps);
      nmaps = ncrossmaps; /* There may be two acceptable facemaps, but only one crossing map. */

      if (ncrossmaps == nfacemaps) {

	for(map=0;map<nmaps;map++) { /* We have actually generated an isomorphism! */
	
	  pd_iso_t *new_iso;
	  new_iso = calloc(1,sizeof(pd_iso_t)); assert(new_iso != NULL); 
	  
	  new_iso->compperm = pd_copy_perm(comp_perm);
	  new_iso->edgemap  = pd_copy_edgemap(edgemaps[edgemap]);
	  new_iso->crossmap = pd_copy_crossmap(crossmaps[map]);
	  new_iso->facemap  = pd_copy_facemap(facemaps[map]);
	  
	  assert(pd_diagram_isotopy_ok(new_iso,pdA,pdB)); /* Check ok-ness as we generate */
	  pd_addto_container(isos,new_iso);

	}
	
      } else {

	assert(ncrossmaps == 1 && nfacemaps == 2); 

	pd_iso_t *new_iso;
	new_iso = calloc(1,sizeof(pd_iso_t)); assert(new_iso != NULL); 
	
	new_iso->compperm = pd_copy_perm(comp_perm);
	new_iso->edgemap  = pd_copy_edgemap(edgemaps[edgemap]);
	new_iso->crossmap = pd_copy_crossmap(crossmaps[0]);

	if (crossmaps[0]->or == PD_POS_ORIENTATION) { 

	  new_iso->facemap  = pd_copy_facemap(facemaps[0]);

	} else { 

	  new_iso->facemap = pd_copy_facemap(facemaps[1]);

	}
	
	assert(pd_diagram_isotopy_ok(new_iso,pdA,pdB)); /* Check ok-ness as we generate */
	pd_addto_container(isos,new_iso);

      }

      pd_free_facemaps(nfacemaps,&facemaps);
	  
    } /* ends "if (ncrossmaps > 0) { " case where we actually tried to generate facemaps */
    
    pd_free_crossmaps(ncrossmaps,&crossmaps);
      
  } /* End of edgemap loop */

  pd_free_edgemaps(nedgemaps,&edgemaps);
    
  /* Now we need to move the data from the container into a simple array */

  pd_iso_t **isobuf;
  
  *nisos = pd_container_nelts(isos);
  
  if (*nisos != 0) { /* This is not guaranteed, since pdA and pdB may not BE isomorphic. */

    isobuf = calloc(*nisos,sizeof(pd_iso_t *)); assert(isobuf != NULL);

    unsigned int i;
    for(i=0;i<*nisos;i++) { isobuf[i] = (pd_iso_t *)(pd_pop_container(isos)); }
    /* Transfer the isos (and their associated memory) to isobuf. */ 

    assert(pd_isos_unique(*nisos,isobuf));

  } else {

    isobuf = NULL;

  }
  
  /* When we're done this loop, there should be nothing in the container. */
  assert(pd_container_nelts(isos) == 0);
  pd_free_container(&isos,pd_free_fake); /* pd_free_fake should never be called. */
  pd_free_perm((void **)(&comp_perm));

  return isobuf;    
}

bool pd_diagram_isotopic(pd_code_t *pdA,pd_code_t *pdB) 

/* This is basically just a wrapper for pd_build_diagram_isotopies. */

{
  pd_iso_t **iso_buf;
  unsigned int nisos;

  iso_buf = pd_build_diagram_isotopies(pdA,pdB,&nisos);

  assert((nisos == 0 && iso_buf == NULL) || (nisos != 0 && iso_buf != NULL));

  if (iso_buf != NULL) { 

    pd_free_isos(&nisos,&iso_buf);
    return true;

  } 

  return false;
  
}
