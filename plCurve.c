/*
 *  @COPYRIGHT@
 * 
 *  Routines to create, destroy, read and write links (and plines)
 * 
 *  $Id: plCurve.c,v 1.9 2004-01-28 18:56:40 cantarel Exp $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#ifdef I_MALLOC
#include <malloc.h>
#endif
#include "config.h"
#include "stdarg.h"
#include "ropelength.h"

/*
 * Set up a new pline.  Pl should point to an *ALREADY ALLOCATED* pline (but
 * with an unallocated space for vertices).  The number of vertices is given in
 * nv and acyclic is set to TRUE or FALSE depending on whether the pline is
 * open or closed.                                        
 *
 */
static void pline_new(pline *Pl,int nv, int acyclic) {
 
  /*const char fn[10] = "pline_new";
    int   i;*/

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
  /* const char fn[10] = "link_new"; */
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
  /* const char fn[10] = "pline_free";*/
  
  if (Pl == NULL) {
    return;
  }

  Pl->nv = 0;
  if (Pl->vt != NULL) {
    free(Pl->vt);
    Pl->vt = NULL;
  }
} /* pline_free */

/*
 * Free the memory associated with a given link.  We then set all the values in
 * the link data structure to reflect the fact that the memory has been freed.
 * We can call link_free twice on the same link without fear. 
 *
 */ 
void link_free(link *L) {
  /* const char fn[10] = "link_free"; */
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

  /* const char fn[10] = "link_write";*/ 
  int i,j;              /* Counter for the for loops */ 
  int nverts =0;        /* Total number of vertices of all components */

  /* First, do a little sanity checking. */

  if (L == NULL) {

    fprintf(stderr,"link_write: Passed NULL pointer as link. \n");
    exit(2);

  }

  if (file == NULL) {

    fprintf(stderr,"link_write: Passed NULL pointer as file.\n");
    exit(2);

  }

  /* Now we begin work. */
    

  for(i=0;i<L->nc;i++) {

    nverts += L->cp[i].nv;

  }

  /* We are ready to write the link. */

  fprintf(file,"VECT \n");
  fprintf(file,"%d %d %d \n",L->nc,nverts,0);
  
  for(i=0;i<L->nc;i++) {

    if (L->cp[i].acyclic) {

      fprintf(file,"%d ",L->cp[i].nv); 
    
    } else {

      fprintf(file,"%d ",-L->cp[i].nv);

    }

  }
  
  fprintf(file,"\n");

  for(i=0;i<L->nc;i++) {

    fprintf(file,"0 ");

  }

  fprintf(file,"\n");

  /* Now we write the vertex data. */

  for(i=0;i<L->nc;i++) {

    for(j=0;j<L->cp[i].nv;j++) {

      fprintf(file,"%g %g %g \n",
	      L->cp[i].vt[j].c[0],
	      L->cp[i].vt[j].c[1],
	      L->cp[i].vt[j].c[2]);

    }

  }

  /* And we're done. */

  return 0;

}


/* The next section of the library file includes some (private) procedures for reading link
   data reliably from Geomview VECT files. We also add a "color" structure to store the color 
   information that may be present in the files, though we don't do anything with it as yet. */

int skip_whitespace_and_comments(FILE *infile)

     /* Procedure positions the file pointer on */
     /* next non-whitespace character, returning */
     /* FALSE if EOF happens first. We skip anything */
     /* between a # and a newline. */

{
  int thischar,commentflag = {FALSE};

  /* First, we check to make sure that infile looks legit. */

  if (infile == NULL) {

    fprintf(stderr,"skip_whitespace_and_comments: infile is a null pointer.\n");
    exit(2);
    
  }
  
  /* Now we start to work. */

  for(;;) {

    thischar = fgetc(infile);

    if (thischar == EOF) {  /* Reached end of file before a non-space, non-comment */
      
      return 0;
    
    } else if (thischar == '#') { /* Started a comment. */
    
      commentflag = TRUE;

    } else if (thischar == '\n' && commentflag) { /* End a comment. */

      commentflag = FALSE;

    } else if (!isspace(thischar) && !commentflag) { /* Found a hit! */

      ungetc(thischar,infile);
      return 1;

    } /* It must have been a space or a non-space in a comment. */

  }

}

int scandoubles(FILE *infile,int ndoubles, ...)

     /* Procedure scans for nfloats floating point (or double) numbers, ignoring 
	whitespace and comments between them. We expect the variable length arguments 
	to contain a collection of pointers to doubles. If not, there's trouble. */

{

  int nconverted = 0,i;
  va_list ap;
  double *thisdouble;

  /* First, we check for overall sanity. */

  if (infile == NULL) {

    fprintf(stderr,"scandoubles: infile is a null pointer.\n");
    exit(2);

  }

  if (ndoubles < 1) {

    fprintf(stderr,"scandoubles: ndoubles (%d) is less than one.\n",ndoubles);
    exit(2);

  }

  va_start(ap,ndoubles);

  /* Now we're ready to work. */

  for (i=0;i<ndoubles;i++) {	/* We expect to exit from the loop by */
                                /* returning, but this is a safety.   */
    
    if (skip_whitespace_and_comments(infile) == 0) { /* Failed */

      return nconverted;

    }

    thisdouble = va_arg(ap,double *);

    if (fscanf(infile,"%lf",thisdouble) != 1) {	/* We couldn't scan. */

      return nconverted;	/* So give up here */

    } else {			/* Else record our victory */

      nconverted++;

    }

  }

  va_end(ap);

  return nconverted;

}

int scanints(FILE *infile,int nints, ...)

     /* Procedure scans for nints integers, ignoring whitespace and
	comments between them. We expect the variable length arguments
	to contain a collection of pointers to doubles. If not,
	there's trouble. */

{

  int nconverted = 0,i;
  va_list ap;
  int *thisint;

  /* First, we check for overall sanity. */

  if (infile == NULL) {

    fprintf(stderr,"scanints: infile is a null pointer.\n");
    exit(2);

  }

  if (nints < 1) {

    fprintf(stderr,"scanints: nints (%d) is less than one.\n",nints);
    exit(2);

  }

  va_start(ap,nints);

  /* Now we're ready to work. */

  for (i=0;i<nints;i++) {	/* We expect to exit from the loop by */
                                /* returning, but this is a safety.   */
    
    if (skip_whitespace_and_comments(infile) == 0) { /* Failed */

      return nconverted;

    }

    thisint = va_arg(ap,int *);

    if (fscanf(infile,"%d",thisint) != 1) {	/* We couldn't scan. */

      return nconverted;	/* So give up here */

    } else {			/* Else record our victory */

      nconverted++;

    }

  }

  va_end(ap);

  return nconverted;

}

int skipwhitespace(FILE *infile)

     /* Procedure skips whitespace (but NOT newlines). */
     /* Returns 1 if we find something, 0 if EOF. */
{
  int thischar;

  for(;;) {

    thischar = fgetc(infile);
    
    if (thischar == '\n' || !isspace(thischar)) {

      ungetc(thischar,infile);
      return 1;

    } else if (thischar == EOF) {

      return 0;

    }

  }

}

/*
 * Read a Geomview VECT file and create a link.  Color information is not 
 * preserved.  File is assumed to be open for reading. Returns either a 
 * pointer to a newly allocated link structure (don't forget to FREE it!)
 * or NULL on failure. */

link *link_read(FILE *file) 
{

  link *L;
  int nverts, ncomp, ncolors;
  int *nvarray, *acyclic;
  int i, scratch, j;
  
  /* First, we check for the 'VECT' keyword. */

  if (fscanf(file," VECT ") == EOF) {

    return NULL;

  }

  /* Now we read the three integers giving vertices, components, and colors. */

  if (scanints(file,3,&ncomp,&nverts,&ncolors) != 3) {

    return NULL;

  }

  /* We now try to read the array of numbers of vertices. */

  nvarray = (int *)calloc(ncomp,sizeof(int));
  acyclic = (int *)calloc(ncomp,sizeof(int));

  for(i=0;i<ncomp;i++) {

    if (scanints(file,1,&(nvarray[i])) != 1) {

      return NULL;

    }

    if (nvarray[i] < 0) {

      acyclic[i] = FALSE;             /* A negative number of vertices indicates a CLOSED component. */
      nvarray[i] *= -1;

    } else {

      acyclic[i] = TRUE;
    
    }

  }

  /* We have now set nvarray and acyclic, and are ready to read (and discard) the color data. */

  for(i=0;i<ncomp;i++) {

    if (scanints(file,1,&scratch) != 1) {

      return NULL;

    }

  }

  /* We now allocate the link data structure. */

  L = link_new(ncomp,nvarray,acyclic);

  if (L == NULL) {   /* If we don't have this much memory, then return NULL. */

    return NULL;

  }

  /* And get ready to read the actual data. */

  for(i=0;i<ncomp;i++) {

    for(j=0;j<L->cp[i].nv;j++) {

      if (scandoubles(file,3,&L->cp[i].vt[j].c[0],&L->cp[i].vt[j].c[1],&L->cp[i].vt[j].c[2]) != 3) {

	return NULL;

      }

    }

  }

  /* We could now read the colors, but we don't bother. */

  return L;
  
}

int link_verts(link *L)

     /* Procedure returns the total number of verts in L. */

{
  int i,verts = 0;

  for(i=0;i<L->nc;i++) {

    verts += L->cp[i].nv;

  }

  return verts;
}

int comp_edges(link *L, int comp) 

     /* Procedure returns the number of edges in component <comp> of <L>. */

{
  if (comp < 0 || comp > L->nc-1) {

    fprintf(stderr,"comp_edges: Can't access component %d of a %d component link.\n",
	    comp,L->nc);
    exit(1);

  }

  return L->cp[comp].nv - ((L->cp[comp].acyclic == TRUE) ? 1:0);

}


int link_edges(link *L) 

     /* Procedure returns the total number of edges in link. */

{
  int sum = 0, i;

  for(i=0;i<L->nc;i++) {

    sum += comp_edges(L,i);

  }

  return sum;
}


link *hopf_link(int verts_per_comp)

     /* Generates a carefully-built Hopf link with <verts_per_comp> vertices in each ring. */

{
  int i, nv[2], acyclic[2] = {FALSE,FALSE};
  double theta, t_step;
  double pi = 3.14159265358979;
  link   *L;

  nv[0] = nv[1] = verts_per_comp;
  
  L = link_new(2,nv,acyclic);

  for(i=0,t_step = 2.0*pi/(double)(verts_per_comp),theta = t_step/2.0;
      i<verts_per_comp;
      i++,theta += t_step) {

    link_v(L,0,i)->c[1] = (1/cos(t_step/2.0))*cos(theta);
    link_v(L,0,i)->c[0] = (1/cos(t_step/2.0))*sin(theta);
    link_v(L,0,i)->c[2] = 0;

    link_v(L,1,i)->c[0] = 0;
    link_v(L,1,i)->c[1] = -(1/cos(t_step/2.0))*cos(theta) + 1.0;
    link_v(L,1,i)->c[2] = (1/cos(t_step/2.0))*sin(theta);

  }

  return L;

}
