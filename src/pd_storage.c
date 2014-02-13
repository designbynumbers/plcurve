/* 

   pd_storage.c : Efficient storage for pd codes. 

*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
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

#ifdef HAVE_STRING_H
  #include<string.h>
#endif

#ifdef HAVE_STDLIB_H
  #include<stdlib.h>
#endif

#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

/* In this version of the world, Judy.h is an installed header
   under libpdcode/, or a local header in randomdiagram/judy/src. */

#include<Judy.h>

/* This depends on a carefully constructed stack of
   headers. */

//#include<thrift/Thrift.h>
//#include<libcassie/cassie.h>

#include<plcTopology.h>
#include<pd_container.h>

#include<pd_multidx.h>
#include<pd_dihedral.h>
#include<pd_perm.h>
  
#include<pd_isomorphisms.h>
#include<pd_storage.h>

struct pdstorage_struct {

  Pvoid_t PJSLArray;

  unsigned int nelts;
  
  char         iter_hash[2*PD_HASHSIZE];      /* Extra-large; should only ever need to contain 32 chars */
  pd_uid_t     iter_uid;
  Word_t       *iter_PValue_hash;  /* Pointer to the current JudyL corresponding to iterator uid */
  
};

// typedef struct pdstorage_struct pd_stor_t;

void pdint_init_pdstor_iterstate(pd_stor_t *pdstor)
{
  strcpy(pdstor->iter_hash,"");
  pdstor->iter_uid = 0;
  pdstor->iter_PValue_hash = NULL;
}

pd_stor_t *pd_new_pdstor()

{
  pd_stor_t *pdstor;

  pdstor = (pd_stor_t *)(malloc(sizeof(pd_stor_t)));
  assert(pdstor != NULL);

  pdstor->PJSLArray = (Pvoid_t) NULL;
  pdstor->nelts = 0;

  /* Initialize the iteration state */
  pdint_init_pdstor_iterstate(pdstor);

  return pdstor;

}

void pd_free_pdstor(pd_stor_t **pdstor)
{
  Word_t *PValue_hash; /* Top Level (Hash indexed JudySL) Pointer */
  char    hash[32];
  Word_t  bytes_freed_all_hashes;

  hash[0] = 0;
  JSLF(PValue_hash,(*pdstor)->PJSLArray,hash); /* Find the first entry in the JudySL of hashes */

  if (PD_VERBOSE > 10) { 

    printf("\n"
	   "pd_free_pdstor: Attempting to free pdstor with %d elts.\n",pd_stor_nelts(*pdstor));

  }
  
  for(;PValue_hash != NULL;) {

    /* At this point, PValue_hash should be a pointer to the head pointer of the JudyL of uids with this hash */

    Pvoid_t PJLArray = (Pvoid_t)(*PValue_hash);
    PWord_t PValue_uid;
    Word_t  uid=0,bytes_freed_this_hash;
    pd_uid_t min_uid,max_uid;
    Word_t   uids_this_hash;

    JLF(PValue_uid,PJLArray,uid); /* Get the first entry in this JudyL of uids */
    min_uid = max_uid = (pd_uid_t)(uid);
    JLC(uids_this_hash,PJLArray,0,-1);

    for(;PValue_uid != NULL;) { /* Iterates in JudyL */

      /* At this point, PValue should be a pointer to a pointer an actual pd_code_t. */
      /* Actually free the data (and clear the pointer too, while we're at it) */
 
      assert((pd_code_t *)(*PValue_uid) != NULL); 
      free((pd_code_t *)(*PValue_uid));
      *PValue_uid = (Word_t)(NULL);
      max_uid = (pd_uid_t)(uid) > max_uid ? (pd_uid_t)(uid) : max_uid;

      JLN(PValue_uid,PJLArray,uid); /* Advance to the next uid */
      
    }

    JLFA(bytes_freed_this_hash,PJLArray); /* Free the actual JudyL memory */
    *PValue_hash = (Word_t)(PJLArray);    /* Propagate updated JudyL pointer back to top JudySL */

    if (PD_VERBOSE > 10) { 

      printf("pd_free_pdstor: \t Freed uids %d-%d (of %d uids total)/%d bytes of memory for hash %s.\n",
	     (int)(min_uid),(int)(max_uid),(int)(uids_this_hash),(int)(bytes_freed_this_hash),hash);

    }

    JSLN(PValue_hash,(*pdstor)->PJSLArray,hash); /* Advance to the next hash string. */

  }

  JSLFA(bytes_freed_all_hashes,(*pdstor)->PJSLArray); /* Free the top-level JudySL memory */

  if (PD_VERBOSE > 10) { 

    printf("pd_stor_free: Freed %d bytes of memory in JudySL for all hashes.\n\n",
	   (int)(bytes_freed_all_hashes));
    
  }

  free(*pdstor);
  *pdstor = NULL;
  
}
    
unsigned int pd_stor_nelts(pd_stor_t *pdstor)
{
  if (pdstor == NULL) { return 0; }
  else { return pdstor->nelts; }
}

void pd_copyinto_pdstor(pd_stor_t *pdstor,pd_code_t *pd)

/* Given a pd_code copy it into the pdstor <=> it's not isomorphic to an existing elt. */ 
/* Assign a uid to any new insert, and update the uid of an existing insert. */

{
  PWord_t PValue_hash;
  Pvoid_t PJLArray;    /* The JudyL corresponding to this hash. */

  /* Step 0: Even calling this (whether or not we succeed) clears the iteration state */
  pdint_init_pdstor_iterstate(pdstor);

  /* Step 1: Try to insert the hash in the top level JudySL. If this hash is 
     already present, we won't interfere with the underlying JudyL. 

     This takes a little unpacking: the Value stored under
     PValue_hash in the PJSLArray attached to this
     pd_stor_t is the root POINTER to the PJLArray of uids
     for this hash. */
  
  JSLI(PValue_hash,pdstor->PJSLArray,/*(uint8_t *)*/(pd->hash));
  PJLArray = (Pvoid_t)(*PValue_hash); 

  /* PJLArray now maps to a JudyL (maybe an empty one!).
 
     We're depending here on the behavior that on a successful 
     insert, libJudy should initialize the value pointed to be
     PValue_hash to zero (NULL). 

     PJLArray should now be a fully functional JudyL.  We
     now need to compare pd with all the pd's in this
     JudyL to see if pd is isomorphic to someone we're
     already storing. */

  Word_t  Index;
  PWord_t PValue_uid;

  if (PJLArray != NULL) { /* If there's a list, we must search it for an isomorphic PD. */

    Index = 0;
        
    JLF(PValue_uid,PJLArray,Index);  /* Initialize the loop at uid (Index) 0 */
    assert(PValue_uid != NULL);      /* If the JudyL is here, it's supposed to contain SOMETHING */    

    Word_t last_ok_Index;

    for(;PValue_uid != NULL;) {

      last_ok_Index = Index; /* This index is present, so record it. */
      
      if (pd_isomorphic((pd_code_t *)(*PValue_uid),pd)) { /* We found it: Quit! */
	
	if (PD_VERBOSE > 20) {
	  
	  printf("\t pd_insert_pdstor: Attempting to insert pd code with hash %s and uid %d.\n"
		 "\t                   Found isomorphic pd code with hash %s and uid %d present.\n"
		 "\t                   Resetting original uid (%d) to correct uid (%d).\n",
		 pd->hash,
		 (int)(pd->uid),
		 ((pd_code_t *)(*PValue_uid))->hash,
		 (int)((pd_code_t *)(*PValue_uid))->uid,
		 (int)(pd->uid),
		 (int)((pd_code_t *)(*PValue_uid))->uid);
	  
	}

	assert( (pd_uid_t)(Index) == ((pd_code_t *)(*PValue_uid))->uid ); /* Check that uid is correctly set. */
	pd->uid = (pd_uid_t)(Index);
	
	return;
	
      }
      
      JLN(PValue_uid,PJLArray,Index); /* Iterate to the next uid. */
      /* This will destroy "Index" if we are at the end of the list */

    }

    /* Having just ended our loop (without finding an isomorphic pd) */
    /* We reset "Index" to the last ok index + 1 in order to find a new uid. */

    Index = last_ok_Index + 1;

  } else {

    Index = 1;

  }

  /* We looped all the way through the JudyL of uids (or
     the JudyL was empty) and didn't find a pd code
     isomorphic to this one.

     It's time to copy the pd code and insert a pointer 
     to the copy. We start by obtaining a pointer which 
     is tied to the new uid. */

  pd_uid_t new_uid;
  new_uid = Index;

  JLI(PValue_uid,PJLArray,(Word_t)(new_uid));
  assert(*PValue_uid == 0); /* If this isn't a clean insert, crash out. */

  /* Now we make a new-memory copy of the pd code */
  
  pd_code_t *pdA;
  pdA = pd_copy(pd);
  pdA->uid = new_uid;

  /* Now we do the insert, storing the POINTER to the new
     pd code in the Value spot POINTED TO by
     PValue_uid. */

  *PValue_uid = (Word_t)(pdA);  

  /* Finally, we have to update the value of PValue_hash to reflect the new PJLArray value */

  *PValue_hash = (Word_t)(PJLArray);

  /* And update the nelts count. */

  pdstor->nelts++;

  /* We also need to report the new uid back to the user. */

  pd->uid = new_uid;

  /* Now report to the user if desired. */

  if (PD_VERBOSE > 20) { 

    printf("\t pd_insert_pdstor: Inserted pd with hash %s and (new) uid %d into storage.\n",
	   pdA->hash,(int)(pdA->uid));

  }

}

pd_code_t *pd_stor_firstelt(pd_stor_t *pdstor)

/* Find the first element in the storage structure. */

{
  /* Step 0: Initialize the iteration state */
  pdint_init_pdstor_iterstate(pdstor);

  /* Step 1: Find the first hash */

  JSLF(pdstor->iter_PValue_hash,pdstor->PJSLArray,
       /*(uint8_t *)*/(pdstor->iter_hash));  /* Find the first hash */

  if (pdstor->iter_PValue_hash == NULL) { /* The PJSLArray is actually empty! */ 

    pdstor->iter_uid = 0; 
    return NULL;

  } else { 

    /* Step 2: Find the first uid in the first hash */

    PWord_t PValue_uid;
    Pvoid_t PJLArray = (Pvoid_t)(*(pdstor->iter_PValue_hash)); /* This is now the JudyL of uids. */
    Word_t  Index = 0;
    
    JLF(PValue_uid,PJLArray,Index); /* Find the first uid. */

    assert(PValue_uid != NULL); /* If the hash is present in the JudySL,
				   it should contain at least one uid! */
    pdstor->iter_uid = (pd_uid_t)(Index);

    return (pd_code_t *)(*PValue_uid); /* PValue_uid points to a pointer to this pd_code. */

  }

}

pd_code_t *pd_stor_nextelt(pd_stor_t *pdstor)

/* Finds the next element in the pdstor (using the current iterstate of the pdstor) */
/* or return NULL. */

{
  /* First, we make sure that firstelt has been called successfully. */

  assert(pdstor->iter_PValue_hash != NULL);
  assert(pdstor->iter_hash[0] != 0); 

  /* Now we search for the next uid in the current hash's JudyL. */

  PWord_t PValue_uid;
  Word_t  Index = (Word_t)(pdstor->iter_uid);
  Pvoid_t PJLArray = (Pvoid_t)(*(pdstor->iter_PValue_hash));

  JLN(PValue_uid,PJLArray,Index);

  if (PValue_uid == NULL) { 

    /* We've reached the end of the current hash's JudyL and need to increment hashes */

    Word_t *PValue_hash;
    JSLN(PValue_hash,pdstor->PJSLArray,/*(uint8_t *)*/(pdstor->iter_hash));

    if (PValue_hash == NULL) {

      /* That was actually the last hash! Reset the iteration state and return NULL. */

      pdint_init_pdstor_iterstate(pdstor);
      return NULL;

    } else {

      /* We've found a new hash. Update the iterstate. */
      pdstor->iter_PValue_hash = PValue_hash;
      PJLArray = (Pvoid_t)(*pdstor->iter_PValue_hash);

      /* Now go ahead and find the first uid in this JudyL */

      Index = 0;
      JLF(PValue_uid,PJLArray,Index);

      assert(PValue_uid != NULL);  /* If the hash is present in the JudySL, 
				      it shouldn't be an empty JudyL. */
      
      pdstor->iter_uid = (pd_uid_t)(Index);
      return (pd_code_t *)(*PValue_uid);
  
    }

  } else { /* We found another uid in this hash's JudyL */

    pdstor->iter_uid = (pd_uid_t)(Index);
    return (pd_code_t *)(*PValue_uid);

  }

  assert(1 == 0); /* We shouldn't get here, so make sure we crash. */

}

pd_code_t *pd_search_pdstor_by_isomorphism(pd_stor_t *pdstor,pd_code_t *pd,
					   pd_iso_t ***pisos,unsigned int *nisos) 

/* Return a pointer to the actual stored copy of the
   unique pd in pdstor which is isomorphic to the given
   pd, along with a buffer of (pointers to) all the
   isomorphisms from pd to the stored copy. Returns NULL
   if we can't find such a pd in storage. */
{
  PWord_t PValue_hash;
  Pvoid_t PJLArray;

  JSLG(PValue_hash,pdstor->PJSLArray,/*(uint8_t *)*/(pd->hash)); /* Search for a matching hash */

  if (PValue_hash == NULL) { /* We didn't find one */

    *pisos = NULL; *nisos = 0; 
    return NULL;

  }

  PJLArray = (Pvoid_t)(*PValue_hash); /* This should now point to the JudyL of pds with matching hash */
  assert(PJLArray != NULL);           /* If the hash is stored, it's supposed to store SOMETHING. */

  Word_t  Index = 0;
  PWord_t PValue_uid;

  JLF(PValue_uid,PJLArray,Index);
  assert(PValue_uid != NULL);    /* As before, we are supposed to have some index present. */
  
  for(;PValue_uid != NULL;) {

    pd_code_t    *stored_pd = (pd_code_t *)(*PValue_uid);
    pd_iso_t     **isos;

    isos = pd_build_isos(pd,stored_pd,nisos);
    
    if (*nisos != 0) { /* We found it! */

      *pisos = isos;
      return stored_pd;

    }

    JLN(PValue_uid,PJLArray,Index); /* Iterate to the next uid. */
    /* This will destroy "Index" if we are at the end of the list */

  }

  /* If we got here, it means that we exhausted the hash
     without finding an isomorphic pd. In this case, the
     search has failed and we're returning NULLs. */

  *pisos = NULL; *nisos = 0;
  return NULL;

}


void pd_stor_stats(pd_stor_t *pdstor,unsigned int *nhashes,unsigned int *nelts)

/* Returns some statistics on pdstor gained by examining the data. */
/* nhashes contains the number of distinct hashes present in the top JudySL */
/* nelts contains a count of the number of elements gained from actual inspection */

{
  Word_t  *PValue_hash; 
  char     Index[2*PD_HASHSIZE];
  Index[0] = 0;
  
  JSLF(PValue_hash,pdstor->PJSLArray,Index);
  
  if (PValue_hash == NULL) { /* The entire thing is empty */

    *nhashes = 0;
    *nelts = 0;

    return;

  }

  for(*nhashes=0,*nelts=0;
      PValue_hash != NULL;) {

    (*nhashes)++;

    Word_t  Rc_word;
    Pvoid_t PJLArray = (Pvoid_t)(*PValue_hash); /* Point to the current JudyL */
    JLC(Rc_word,PJLArray,0,-1);
    *nelts += (unsigned int)(Rc_word);

    JSLN(PValue_hash,pdstor->PJSLArray,Index); /* Iterate to the next hash */
    
  }

}

void pd_display_pdstor(FILE *stream,pd_stor_t *pdstor) /* Prints a representation to stream. */

{
  unsigned int nelts,nhashes;

  pd_stor_stats(pdstor,&nhashes,&nelts);

  fprintf(stream,
	  "pdstor \n"
	  "nelts %u/%u (claimed/actual) nhashes %u \n",pd_stor_nelts(pdstor),nelts,nhashes);

  /* Now list each hash and uid separately using iterators. */

  pd_code_t *pd;
  char       cur_hash[2*PD_HASHSIZE]; 
  pd_uid_t   uid_min,uid_max;

  pd = pd_stor_firstelt(pdstor);
  strncpy(cur_hash,pd->hash,PD_HASHSIZE);
  uid_min = pd->uid;
  uid_max = pd->uid;

  fprintf(stream,"\t hash %s uids %u-",cur_hash,(unsigned int)(uid_min));

  for(;pd != NULL;pd = pd_stor_nextelt(pdstor)) { 

    if (!strcmp(cur_hash,pd->hash)) { /* The hash hasn't changed, so update uid_max */

      uid_max = (pd->uid > uid_max) ? pd->uid : uid_max;

    } else { 

      fprintf(stream,"%u \n",(unsigned int)uid_max);
      
      strncpy(cur_hash,pd->hash,32);
      uid_min = pd->uid;
      uid_max = pd->uid;
      
      fprintf(stream,"\t hash %s uids %u-",cur_hash,(unsigned int) uid_min);

    }

  }

  fprintf(stream,"%u \n",(unsigned int) uid_max);
      
  /* Now wrap up the display. */

}
  
void pd_write_pdstor(FILE *stream,pd_stor_t *pdstor) /* Prints a representation to stream. */

{
  assert(stream != NULL);
  assert(pdstor != NULL);

  /* First, generate some metadata at the top of the file */

  /* fprintf(stream, */
  /* 	  "# generated by %s (build %s, %s) ", */
  /* 	  PACKAGE_VERSION, */
  /* 	  __DATE__,__TIME__); */

  /* char svntag[1024]; */

  /* sprintf(svntag,"%s",SVNVERSION); */
  /* if (!strstr("exported",svntag)) {  /\* We were built from svn *\/ */
  /*   fprintf(stream,"svn version %s\n",SVNVERSION); */
  /* } */

  /* Next, describe the actual data to come. */

  unsigned int nelts,nhashes;
  pd_stor_stats(pdstor,&nhashes,&nelts);

  fprintf(stream,
	  "pdstor \n"
	  "nelts %u/%u (claimed/actual) nhashes %u \n\n",pd_stor_nelts(pdstor),nelts,nhashes);

  if (nelts == 0) { return; }
 
  /* Now iterate through the pdstor, writing each pd to disk */

  pd_code_t *pd;
  for(pd = pd_stor_firstelt(pdstor);pd != NULL;pd = pd_stor_nextelt(pdstor)) {

    pd_write(stream,pd);
    fprintf(stream,"\n");

  }

  if (PD_VERBOSE > 20) { 

    printf("pd_write_pdstor: Wrote %d element, %d hash pdstor to stream.\n",
	   nelts,nhashes);

  }

}

pd_stor_t *pd_read_pdstor(FILE *stream) 

{
  assert(stream != NULL);

  /* First, try to clear the metadata line */

  /* Now read the header */

  unsigned int nelts_claimed,nelts_actual;
  unsigned int nhashes; 
  
  if (fscanf(stream,"pdstor nelts %u/%u (claimed/actual) nhashes %u",&nelts_claimed,&nelts_actual,&nhashes) != 3) { 

    if (PD_VERBOSE > 20) {

      printf("pd_read_pdstor: Couldn't parse header.\n");
      
    }

    return NULL;

  }

  unsigned int nelts;

  if (nelts_claimed != nelts_actual) { 

    if (PD_VERBOSE > 20) { printf("pd_read_pdstor: Header of pdstor file appears corrupt.\n"); return NULL; }

  }

  nelts = nelts_claimed;

  /* Now we go into a loop to build the pdstor */

  pd_stor_t *pdstor = pd_new_pdstor();

  unsigned int elt;

  for(elt=0;elt<nelts && !feof(stream);elt++) {

    if (PD_VERBOSE > 50) { /* Actually draw a progress bar */

      printf("pd_read_pdstor: Read %u of %u pdcodes.\r",elt,nelts);

    }

    pd_code_t pd;

    if (!pd_read(stream,&pd)) { 

      if (PD_VERBOSE > 20) { printf("pd_read_pdstor: Couldn't parse pd code %d.\n",elt); }
      pd_free_pdstor(&pdstor);
      return NULL;

    } else {

      pd_copyinto_pdstor(pdstor,&pd);

    }

  }

  if (elt != nelts) { 

    if (PD_VERBOSE > 20) { printf("pd_read_pdstor: Only read %u of %u elements from file.\n",elt,nelts); }
    pd_free_pdstor(&pdstor);
    return NULL;

  }

  if (PD_VERBOSE > 20) { 

    printf("pd_read_pdstor: Read %u elements from stream.\n",nelts);
    
  }

  /* Check stats. */

  unsigned int built_nhashes, built_nelts;
  pd_stor_stats(pdstor,&built_nhashes,&built_nelts);

  assert(built_nhashes == nhashes);
  assert(built_nelts == nelts);

  return pdstor;

}
  
void pd_copyinto_cass(pd_code_t *pd)
{

  //int com_size = 21 + PD_MAXVERTS*4 + 1 + 32;
  char *com_string[110];
  char *current_char[1];
  
  sprintf(com_string, "python shad_stor.py ");
  
  int i;
  int j;
  
  for(i=0; i<pd->ncross; i++){
    for(j=0; j<4; j++){
      sprintf(current_char, "%c",(char)(pd->cross[i].edge[j]+65));
      strcat(com_string,current_char);
    }
    
  }
  
  strcat(com_string," ");
  strcat(com_string,pd->hash);
  
  printf("%s \n", com_string );  
  system(com_string);
}
