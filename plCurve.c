/*
 *  @COPYRIGHT@
 * 
 *  Routines to create, destroy, read and write links (and plines)
 * 
 *  $Id: plCurve.c,v 1.1 2003-12-24 00:38:54 ashted Exp $
 *
 */

#include <c.h>
#include <stdio.h>
#ifdef I_MALLOC
#include <malloc.h>
#endif
#include "config.h"
#include "ropelength.h"

/*
 * Set up a new pline.  Pl should point to an *ALREADY ALLOCATED* pline (but
 * with an unallocated space for vertices).  The number of vertices is given in
 * nv and acyclic is set to TRUE or FALSE depending on whether the pline is
 * open or closed.                                        
 *
 */
static inline void pline_new(pline *Pl,int nv, int acyclic) {
  const char fn[10] = "pline_new";
  int   i;

  /* Sanity checking */

  if (nv < 1) {
    fprintf(stderr,"Can't create a pline with %d vertices.\n",nv);
    exit(-1);
  }

  Pl->acyclic = acyclic;
  Pl->nv = nv;
  if ((Pl->vt = (vector *)malloc(nv*sizeof(vector))) == NULL) {
    fprintf(stderr,"Can't allocate space for %d vertices in pline_new.\n",nv);
    exit(-1);
  }
}

/*
 * Procedure allocates memory for a new link. The number of components is given
 * by components. The number of vertices in each component shows up in the
 * buffer pointed to be nv.  The closed/open nature of each pline is given in
 * the array pointed to by acyclic.                           
 *
 */
link *link_new(int components, int *nv, int *acyclic) {
  const char fn[10] = "link_new";
  link *L;
  int   i;

  /* First, we check to see that the input values are reasonable. */

  if (components < 1) {
    fprintf(stderr,"Can't create a link with %d components.",components);
    exit(-1);
  }

  if (nv == NULL || acyclic == NULL) {
    fprintf(stderr,"The vertex or acyclic status list is empty.");
    exit(-1);
  }

  /* Now we attempt to allocate space for these components. */
  
  if ((L = (link *)malloc(sizeof(link))) == NULL) {
    fprintf(stderr,"Could not allocate space for link in link_new.\n");
    exit(-1);
  }
  L->nc = components;
  if ((L->cp = (pline *)malloc(L->nc*sizeof(pline))) == NULL) {
    fprintf(stderr,"Can't allocate array of pline ptrs in link_new.\n");
    exit(-1);
  }

  for (i = 0; i < L->nc; i++) {
    pline_new(&L->cp[i],nv[i],acyclic[i]);
  }

  return L;
}

/*
 * Free the memory used to hold vertices in a given pline (not the memory of
 * the pline itself).  We then set all the values in the link data structure to
 * reflect the fact that the memory has been freed.  We can call pline_free
 * twice on the same pline without fear. 
 *
 */ 
void pline_free(pline *Pl) {
  const char fn[10] = "pline_free";
  
  if (Pl == NULL) {
    return;
  }

  Pl->nv = 0;
  if (Pl->vt != NULL) {
    free(Pl->vt);
    Pl->vt == NULL;
  }
} /* pline_free */

/*
 * Free the memory associated with a given link.  We then set all the values in
 * the link data structure to reflect the fact that the memory has been freed.
 * We can call link_free twice on the same link without fear. 
 *
 */ 
void link_free(link *L) {
  const char fn[10] = "link_free";
  int i;

  /* First, we check the input. */
  
  if (L == NULL) {
    return; /* Move along, nothing to see here */
  }

  if (L->nc < 0) {
    fprintf(stderr,"Link appears corrupted. L.nc = %d.",L->nc);
    exit(-1);
  }

  /* Now we can get to work. */

  for (i=0; i<L->nc; i++) {
    pline_free(&L->cp[i]);
  }

  free(L->cp);
  L->nc = 0;

  free(L);
  L = NULL;
} /* link_free */

/*
 * Writes the link to a file in Geomview VECT format.  The file format is:
 *
 * VECT                           # mandatory keyword
 * Ncomponents Nvertices Ncolors  # total number of components and vertices
 * Nv[0] ... Nv[NPolylines-1]     # number of vertices in each polyline 
 *                                # closed polylines use negative numbers
 * Nc[0] ... Nc[NPolylines-1]     # number of colors for each polyline 
 *
 * Vert[0] ... Vert[Nvertices-1]  # All the vertices, as triples of doubles
 * Color[0] ... Color[NColors]    # All the colors, in RGBA format
 *
 * Comments begin with #, and proceed to the end of the line. They are allowed
 * wherever a newline is allowed.
 *
 * We assume that file is open for writing.  Colors are arbitrarily assigned.
 *
 * Returns TRUE if successful write, FALSE otherwise.
 *
 */
int link_write(FILE *file,link *L) {
  const char fn[10] = "link_write";
  int i;              /* Counter for the for loops */ 
  int totNumVert = 0; /* Total number of vertices of all components */

  /************************************************/
  fprintf(stderr,"link_write has not yet been rewritten");
  exit(-1);
  /************************************************/

//  for(i=0;i<L->nc;i++)                      /* Loop through each component */   
//    totNumVert += L->nv[i];                 /* Running count of total number of vertices in compts */
//  
//  fprintf(file, "VECT\n"); 
//  fprintf(file, "%d %d %d \n",L->nc, totNumVert, L->nc);
//  /* write the number of compts, */
//  /* total number of verts in all components, and */
//  /* total number of colors, with one color per component */
//
//  for(i=0;i<L->nc;i++)                      /* Loop through the total number of components*/
//    fprintf(file, "%d ",L->nv[i]);          /* write the number of vertices in each component */
//
//  fprintf(file, "\n");
//
//  for(i=0;i<L->nc;i++)                      /* Loop through the total number of components */
//    fprintf(file, "%d ", 1);               /* write the number of colors (one) in each component */
//
//  fprintf(file, "\n");
//
//  for(i=0;i<totNumVert;i++)                /* Loop through each vertex */
//    fprintf(file, "%f %f %f\n",L->vt[i]->c[0],L->vt[i]->c[1],L->vt[i]->c[2]);   /* write the coordinates of each vertex */
//
//  fprintf(file, "\n");
//
//  for(i=0;i<L->nc;i++) /* Loop through each component*/
//    fprintf(file, ",%f,%f,%f,", 1.0 - (1.0/i), 0.0 + (1.0/i), 1.0 - (1.0/i));
//  /* assign the color 1/ith of the way down the diagonal of the spectrum */
//  /* to the ith component */
//
//  fprintf(file, "\n");
}

/*
 * Read a Geomview VECT file and create a link.  Color information is not 
 * preserved.  File is assumed to be open for reading.  Returns TRUE on 
 * success, FALSE on failure.
 *
 */
int link_read(FILE *file, link *L) {
  const char fn[10] = "link_read";
  int i;
  int n;
  int trash_int;
  int return_code1;
  int return_code2;
  int span;
  int nc;
  int total_verts;
  int check[] = {FALSE, FALSE};  
  double energy;
  char trash_string[] = "VECT";
  char ch;
  char line[1024];  
  char comment[] = "#";
  fpos_t pos;

  /*****************************************************/
  fprintf(stderr,"link_read has not yet been rewritten");
  exit(-1);
  /*****************************************************/
  
//  if (file == NULL || L == NULL) {
//    fprintf(stderr, "Of the Link and File passed to link_read, one (or both) is a null pointer.\n");
//      return FALSE;
//  }
//
//  while(TRUE)
//    {
//      fgetpos(file, &pos);
//      
//      while(TRUE)
//	{
//	  i = 0;
//	  ch = fgetc(file);
//	  if( ch == '\n' )
//	    {
//	      line[i] = '\0';
//	      break;
//	    }
//	  if( ch == EOF )
//	    {
//	      fprintf(stderr, "NewBreedLink: Error, file does not contain energy and/or link_type data.\n");
//	      return FALSE;
//	    }
//	  line[i++]=ch;
//	}
//      return_code1 = sscanf(line,"# energy %f", energy);    
//      return_code2 = sscanf(line,"# HOMEFLY span %d", span);
//      
//      if (return_code1 == 1) {
//	check[0] = TRUE;
//      }
//
//      if (return_code2 == 1) {
//	fsetpos(file, &pos);
//	fscanf(file, "# HOMFLY span %d coeffs", span);
//	for (n =0; n < span; n++) {
//	  fscanf(file, " %d", link_type[n]);
//	}
//	check[1] = TRUE;
//      } 
//
//      if (check[0] == TRUE && check[1] == TRUE) {
//	break;
//      }
//      
//    }
// 
//  rewind(file);
//  
//  while(TRUE)
//    { 
//      while(TRUE)
//	{
//	  i = 0;
//	  ch = fgetc(file);
//	  if( ch == '\n' )
//	    {
//	      line[i] = '\0';
//	      break;
//	    }
//	  if( ch == EOF )
//	    {
////	      fprintf(stderr, "NewBreedLink: VECT not found. File not proper .vect format");
//	      return FALSE;
//	    }
//	  line[i++]=ch;
//	}
//      if (strcmp(line, trash_string) == TRUE)
//	{
//	    fscanf(file, " %d", nc);
//	    fscanf(file, " %d", total_verts);
//	    fscanf(file, " %d", nc);
//	    for (n = 0; n < nc; n++) {
//	      fscanf(file, " %d", nv[n]);
//	    }
//	    for (n = 0; n < nc; n++) {
//	      fscanf(file, " %d", trash_int);
//	    }
//	    break;
//	}   
//    }
//
//  *L = link_new(nc, nv);
//  L->energy = energy;
//  /* L->span = span; */
//
//  for (i = 0; i < MAX_HOMFLY; i++) {
//    L->link_type[i] = link_type[i];
//  }
//  
//  //TIME TO READ IN THE VERTICES
//
//  for(i=0;i<total_verts;i++) {               /* Loop through each vertex */
//    fgetpos(file, &pos);
//    if (fscanf(file, "%f %f %f\n",L->vt[i]->c[0],L->vt[i]->c[1],L->vt[i]->c[2]) == FALSE){ /* write the coordinates of each vertex */
//      fsetpos(file, &pos);   
//      while(TRUE)
//	{
//	  n = 0;
//	  ch = fgetc(file);
//	  if( ch == '\n' )
//	    {
//	      line[i] = '\0';
//	      break;
//	    }
//	  if( ch == EOF )
//	    {
//	      fprintf(stderr, "NewBreedLink: Error, specified number of vertices is greater than actual number in .vect file.");
//	      return FALSE;
//	    }
//	  line[n++]=ch;
//	}
//      if(strcspn(comment,line) < strlen(line)){
//	i--;
//      } else {
//	fprintf(stderr, "NewBreedLink: Error, invalid information in vertex set at vertex %d.", i);
//	return FALSE;
//      }
//    } 
//  } 
//  return TRUE;
  
}
