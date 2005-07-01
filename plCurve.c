/*
 *  Routines to create, destroy, read and write links (and plines)
 * 
 *  $Id: plCurve.c,v 1.23 2005-07-01 01:56:33 cantarel Exp $
 *
 */

/* Copyright 2004 The University of Georgia. */

/* This file is part of vecttools.
   
vecttools is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

vecttools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vecttools; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#ifdef HAVE_MATH_H
#include <math.h>
#endif

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_STDARG_H
#include <stdarg.h>
#endif
#ifdef HAVE_CTYPE_H
#include <ctype.h>
#endif

#include <link.h>

extern int  linklib_error_num;
extern char linklib_error_str[80];

/*
 * Set up a new pline.  Pl should point to an *ALREADY ALLOCATED* pline (but
 * with an unallocated space for vertices).  The number of vertices is given in
 * nv and acyclic is set to TRUE or FALSE depending on whether the pline is
 * open or closed.                                        
 *
 * We allocate two extra vertices, at -1 and nv to make "wrap-around" much 
 * simpler.
 */
static void pline_new(linklib_pline *Pl,int nv, int acyclic, int cc) {
 
  if (nv < 1) {
    linklib_error_num = 21;
    sprintf(linklib_error_str,"pline_new: Can't create a pline with %d vertices.\n",nv);
    return;
  }

  Pl->acyclic = acyclic;
  Pl->nv = nv;
  if ((Pl->vt = 
       (linklib_vector *)calloc((nv+2),sizeof(linklib_vector))) == NULL) {
    linklib_error_num = 22;
    sprintf(linklib_error_str,"pline_new: Can't allocate space for %d vertices in pline_new.\n",nv);
    return;
  }
  Pl->vt++; /* so that Pl->vt[-1] is a valid space */

  Pl->cc = cc;
  if ((Pl->clr = (linklib_color *)calloc(cc,sizeof(linklib_color))) == NULL) {
    linklib_error_num = 23;
    sprintf(linklib_error_str,"pline_new: Can't allocate space for %d colors in pline_new.\n",cc);
    return;
  }
}

/*
 * Procedure allocates memory for a new link. The number of components is given
 * by components. The number of vertices in each component shows up in the
 * buffer pointed to be nv.  The closed/open nature of each pline is given in
 * the array pointed to by acyclic.                           
 *
 */
linklib_link *linklib_link_new(int components, const int *nv, 
                               const int *acyclic, const int *cc) 
{
  linklib_link *L;
  int i;

  /* First, we check to see that the input values are reasonable. */

  if (components < 1) {
    linklib_error_num = 31;
    sprintf(linklib_error_str,"linklib_link_new: Can't create a link with %d components.",components);
    return NULL;
  }

  if (nv == NULL || acyclic == NULL || cc == NULL) {
    linklib_error_num = 32;
    sprintf(linklib_error_str,"linklib_link_new: nv, acyclic or cc is NULL.");
    return NULL;
  }

  /* Now we attempt to allocate space for these components. */
  
  if ((L = (linklib_link *)malloc(sizeof(linklib_link))) == NULL) {
    linklib_error_num = 33;
    sprintf(linklib_error_str,"linklib_link_new: Could not allocate space for link in link_new.\n");
    return NULL;
  }
  L->nc = components;
  if ((L->cp = (linklib_pline *)malloc(L->nc*sizeof(linklib_pline))) == NULL) {
    linklib_error_num = 34;
    sprintf(linklib_error_str,"Can't allocate array of pline ptrs in link_new.\n");
    return NULL;
  }

  for (i = 0; i < L->nc; i++) {
    pline_new(&L->cp[i],nv[i],acyclic[i],cc[i]);
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
void pline_free(linklib_pline *Pl) {
  
  if (Pl == NULL) {
    return;
  }

  Pl->nv = 0;
  if (Pl->vt != NULL) {
    Pl->vt--; /* undo our original vt++ (for wraparound) */
    free(Pl->vt);
    Pl->vt = NULL;
  }
 
  Pl->cc = 0;
  if (Pl->clr != NULL) {
    free(Pl->clr);
  }
} /* pline_free */

/*
 * Free the memory associated with a given link.  We then set all the values in
 * the link data structure to reflect the fact that the memory has been freed.
 * We can call link_free twice on the same link without fear. 
 *
 */ 
void linklib_link_free(linklib_link *L) {
  int i;

  /* First, we check the input. */
  if (L == NULL) {
    return; /* Move along, nothing to see here */
  }

  if (L->nc < 0) {
    linklib_error_num = 41;
    sprintf(linklib_error_str,"linklib_link_free: Link appears corrupted. L.nc = %d.",L->nc);
    return;
  }

  /* Now we can get to work. */
  for (i=0; i<L->nc; i++) {
    pline_free(&L->cp[i]);
  }

  free(L->cp);
  L->nc = 0;

  free(L);
  L = NULL;
} /* linklib_link_free */

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

int linklib_link_write(FILE *file, const linklib_link *L) {
  int i,j;              /* Counter for the for loops */ 
  int nverts = 0;       /* Total number of vertices of all components */
  int colors = 0;       /* Total number of colors of all components */

  /* First, do a little sanity checking. */
  if (L == NULL) {
    linklib_error_num = 51;
    sprintf(linklib_error_str,"linklib_link_write: Passed NULL pointer as link. \n");
    return -1;
  }

  if (file == NULL) {
    linklib_error_num = 52;
    sprintf(linklib_error_str,"linklib_link_write: Passed NULL pointer as file.\n");
    return -1;
  }

  /* Now we begin work. */
  for(i=0;i<L->nc;i++) {
    nverts += L->cp[i].nv;
    colors += L->cp[i].cc;
  }

  /* We are ready to write the link. */
  fprintf(file,"VECT \n");
  fprintf(file,"%d %d %d \n",L->nc,nverts,colors);
  
  for(i=0;i<L->nc;i++) {
    if (L->cp[i].acyclic) {
      fprintf(file,"%d ",L->cp[i].nv); 
    } else {
      fprintf(file,"%d ",-L->cp[i].nv);
    }
  }
  fprintf(file,"\n");

  for(i=0;i<L->nc;i++) {
    fprintf(file,"%d ",L->cp[i].cc);
  }
  fprintf(file,"\n");

  /* Now we write the vertex data . . . */
  for(i=0;i<L->nc;i++) {
    for(j=0;j<L->cp[i].nv;j++) {
      fprintf(file,"%.16g %.16g %.16g \n", 
	      L->cp[i].vt[j].c[0], L->cp[i].vt[j].c[1],
	      L->cp[i].vt[j].c[2]);
    }
  }

  /* . . . and the color data. */
  for (i=0; i < L->nc; i++) {
    for (j=0; j < L->cp[i].cc; j++) {
      fprintf(file,"%g %g %g %g\n", L->cp[i].clr[j].r, L->cp[i].clr[j].g,
        L->cp[i].clr[j].b, L->cp[i].clr[j].alpha);
    }
  }
      
  /* And we're done. */
  return 0;
}


/* The next section of the library file includes some (private) procedures *
 * for reading link data reliably from Geomview VECT files. We also add a  *
 * "color" structure to store the color information that may be present in *
 * the files, though we don't do anything with it as yet.                  */

/* Procedure positions the file pointer on next non-whitespace character,   *
 * returning FALSE if EOF happens first. We skip anything between a # and a *
 * newline.                                                                 */
int skip_whitespace_and_comments(FILE *infile)
{
  int thischar,commentflag = {FALSE};

  /* First, we check to make sure that infile looks legit. */
  if (infile == NULL) {
    linklib_error_num = 61;
    sprintf(linklib_error_str,"skip_whitespace_and_comments: infile is a null pointer.\n");
    return -1;
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

/* Procedure scans for nfloats floating point (or double) numbers, ignoring  *
 * whitespace and comments between them. We expect the variable length       *
 * arguments to contain a collection of pointers to doubles. If not, there's *
 * trouble.                                                                  */
int scandoubles(FILE *infile,int ndoubles, ...)
{
  int nconverted = 0,i;
  va_list ap;
  double *thisdouble;

  /* First, we check for overall sanity. */

  if (infile == NULL) {
    linklib_error_num = 71;
    sprintf(linklib_error_str,"scandoubles: infile is a null pointer.\n");
    return -1;
  }

  if (ndoubles < 1) {
    linklib_error_num = 72;
    sprintf(linklib_error_str,"scandoubles: ndoubles (%d) is less than one.\n",ndoubles);
    return -1;
  }

  va_start(ap,ndoubles);

  /* Now we're ready to work. */

  for (i=0;i<ndoubles;i++) {    /* We expect to exit from the loop by */
                                /* returning, but this is a safety.   */
    if (skip_whitespace_and_comments(infile) == 0) { /* Failed */
      return nconverted;
    }

    thisdouble = va_arg(ap,double *);
    if (fscanf(infile,"%lf",thisdouble) != 1) { /* We couldn't scan. */
      return nconverted;        /* So give up here */
    } else {                    /* Else record our victory */
      nconverted++;
    }
  }
  va_end(ap);

  return nconverted;
}

/* Procedure scans for nints integers, ignoring whitespace and     *
 * comments between them. We expect the variable length arguments  *
 * to contain a collection of pointers to ints. If not,            *
 * there's trouble.                                                */
int scanints(FILE *infile,int nints, ...)
{
  int nconverted = 0,i;
  va_list ap;
  int *thisint;

  /* First, we check for overall sanity. */

  if (infile == NULL) {
    linklib_error_num = 73;
    sprintf(linklib_error_str,"scanints: infile is a null pointer.\n");
    return -1;
  }

  if (nints < 1) {
    linklib_error_num = 74;
    sprintf(linklib_error_str,"scanints: nints (%d) is less than one.\n",nints);
    return -1;
  }

  va_start(ap,nints);

  /* Now we're ready to work. */
  for (i=0;i<nints;i++) {       /* We expect to exit from the loop by */
                                /* returning, but this is a safety.   */
    if (skip_whitespace_and_comments(infile) == 0) { /* Failed */
      return nconverted;
    }
    thisint = va_arg(ap,int *);

    if (fscanf(infile,"%d",thisint) != 1) {     /* We couldn't scan. */
      return nconverted;        /* So give up here */
    } else {                    /* Else record our victory */
      nconverted++;
    }
  }
  va_end(ap);
  return nconverted;
}

/* Procedure skips whitespace (but NOT newlines). */
/* Returns 1 if we find something, 0 if EOF.      */
int skipwhitespace(FILE *infile)
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
 * Touchup the "extra" vertices at each end of the component plines which are
 * used to implement "wraparound".
 *
 */
void linklib_link_fix_wrap(const linklib_link *L) {
  int i,nv;

  for (i = 0; i < L->nc; i++) {
    nv = L->cp[i].nv;
    if (L->cp[i].acyclic) {
      /* fold it back on itself: v_{-1} = v_1 and v_{nv} = v_{nv-2} */
      L->cp[i].vt[-1] = L->cp[i].vt[1];
      L->cp[i].vt[nv] = L->cp[i].vt[nv-2];
    } else {
      /* wrap it around: v_{-1} = v_{nv-1} and v_{nv} = v_0 */
      L->cp[i].vt[-1] = L->cp[i].vt[nv-1];
      L->cp[i].vt[nv] = L->cp[i].vt[0];
    }
  }
}

/*
 * Read a Geomview VECT file and create a link.  Color information is not 
 * preserved.  File is assumed to be open for reading. Returns either a 
 * pointer to a newly allocated link structure (don't forget to FREE it!)
 * or NULL on failure. 
 *
 */
linklib_link *linklib_link_read(FILE *file) 
{
  linklib_link *L;
  int nverts, ncomp, ncolors;
  int *nvarray, *acyclic, *ccarray;
  int i, j;
  int nv;
  
#ifdef DEBUG
  if (linklib_debug_level() > 8) {
    printf("linklib_link_read:Reading link from file.\n");
  }
#endif

  /* First, we check for the 'VECT' keyword. */

  if (fscanf(file," VECT ") == EOF) {
  
    linklib_error_num = 81;
    sprintf(linklib_error_str,"linklib_link_read: Couldn't find VECT keyword.");    
    return NULL;
  }

  /* Now we read the three integers giving vertices, components, and colors. */

  if (scanints(file,3,&ncomp,&nverts,&ncolors) != 3) {
  
    linklib_error_num = 82;
    sprintf(linklib_error_str,"linklib_link_read: Couldn't parse <ncomp> <nverts> <ncolors> line");
    return NULL;
  }

  /* We now try to read the array of numbers of vertices. */

  nvarray = (int *)calloc(ncomp,sizeof(int));
  acyclic = (int *)calloc(ncomp,sizeof(int));
  ccarray = (int *)calloc(ncomp,sizeof(int));

  for(i=0;i<ncomp;i++) {
    if (scanints(file,1,&(nvarray[i])) != 1) {
      linklib_error_num = 83;
      sprintf(linklib_error_str,"linklib_link_read: Couldn't parse number"
              "of vertices in component %d.",i);    
      return NULL;
    }
    if (nvarray[i] < 0) {
      /* A negative number of vertices indicates a CLOSED component. */
      acyclic[i] = FALSE;  
      nvarray[i] *= -1;
    } else {
      acyclic[i] = TRUE;
    }
  }

  /* We have set nvarray and acyclic and are ready to read the color data.  */

  for(i=0;i<ncomp;i++) {
    if (scanints(file,1,&(ccarray[i])) != 1) {
      linklib_error_num = 84;
      sprintf(linklib_error_str,"linklib_link_read: Couldn't parse <ncolors>"
      "for component %d.", i); 
      return NULL;
    }
  }

  /* We now allocate the link data structure. */

  L = linklib_link_new(ncomp,nvarray,acyclic,ccarray);

  if (L == NULL) {   /* If we don't have this much memory, then return NULL. */
    linklib_error_num = 85;
    sprintf(linklib_error_str,"linklib_link_read: Couldn't allocate enough"
    " memory for link.");
    return NULL;
  }

  /* And get ready to read the actual data. */

  for(i = 0; i < ncomp; i++) {
    nv = L->cp[i].nv;
    for(j = 0; j < nv; j++) {
      if (scandoubles(file,3,&L->cp[i].vt[j].c[0],
                             &L->cp[i].vt[j].c[1],
                             &L->cp[i].vt[j].c[2]) != 3) {
        linklib_link_free(L);
        linklib_error_num = 86;
        sprintf(linklib_error_str,"linklib_link_read: Couldn't parse "
        " <x> <y> <z> data for vertex %d of component %d.",j,i);
        return NULL;
      }
    }
  }
  /* Now set the "wrap-around" vertices */
  linklib_link_fix_wrap(L);

    /* And finally the colors. Unfortunately, to really comply with 
        the Geomview standard here we have to be kind of careful. */

  /* And finally the colors. */
  for (i=0; i < ncomp; i++) {
    for (j=0; j < L->cp[i].cc; j++) {
      if (scandoubles(file,4, &L->cp[i].clr[j].r, &L->cp[i].clr[j].g,
           &L->cp[i].clr[j].b, &L->cp[i].clr[j].alpha) != 4) {

	linklib_error_num = 28;
	sprintf(linklib_error_str,
		"linklib_read: Couldn't parse color %d "
		"in component %d of link.\n",j,i);
        return NULL;


      }
    }
  }

  free(ccarray);
  free(acyclic);
  free(nvarray);

  return L;
}

#define pline_edges(P) (((P).acyclic) ? (P).nv-1 : (P).nv)
/* 
 *   Return the total number of edges in link. 
 */
int linklib_link_edges(const linklib_link *L) 
{
  int i, edges = 0;

  for (i=0;i<L->nc;i++) {
    edges += pline_edges(L->cp[i]);
  }
  return edges;
}

/* 
 * Duplicate a link and return the duplicate 
 *
 */
linklib_link *linklib_link_copy(const linklib_link *L) {
  linklib_link *nL;
  int *nv,*acyclic,*ccarray;
  int cnt,cnt2;

  if ((nv = (int *)malloc((L->nc)*sizeof(int))) == NULL ||
      (acyclic = (int *)malloc((L->nc)*sizeof(int))) == NULL ||
      (ccarray = (int *)malloc((L->nc)*sizeof(int))) == NULL) {
    linklib_error_num = 75;
    sprintf(linklib_error_str,"Unable to malloc space for alternate link.\n");
    return NULL;
  }
  for (cnt = 0; cnt < L->nc; cnt++) {
    nv[cnt] = L->cp[cnt].nv;
    acyclic[cnt] = L->cp[cnt].acyclic;
    ccarray[cnt] = L->cp[cnt].cc;
  }
  nL = linklib_link_new(L->nc,nv,acyclic,ccarray);

  for (cnt = 0; cnt < L->nc; cnt++) {
    for (cnt2 = 0; cnt2 < L->cp[cnt].nv; cnt2++) {
      nL->cp[cnt].vt[cnt2] = L->cp[cnt].vt[cnt2];
    }
    for (cnt2 = 0; cnt2 < L->cp[cnt].cc; cnt2++) {
      nL->cp[cnt].clr[cnt2] = L->cp[cnt].clr[cnt2];
    }
  }

  free(ccarray);
  free(acyclic);
  free(nv);

  return nL;
}

linklib_vector linklib_link_tangent_vector(linklib_link *link,int cp, int vt)

/* Procedure computes a (unit) tangent vector to <link> 
   at the given vertex of the given component. */

{
  linklib_vector in, out, tan;

  if (link->cp[cp].acyclic) {

    if (vt == 0) {

      tan = linklib_vminus(link->cp[cp].vt[1],
		  link->cp[cp].vt[0]);

      linklib_vector_normalize(&tan);

      return tan;

    } else if (vt == link->cp[cp].nv-1) {

       tan = linklib_vminus(link->cp[cp].vt[link->cp[cp].nv-1],
		   link->cp[cp].vt[link->cp[cp].nv-2]);

       linklib_vector_normalize(&tan);

       return tan;

    }

  }

  /* We now know that either we are on a closed 
     component, or we are not at an endpoint.   */
  
   in = linklib_vminus(link->cp[cp].vt[vt+1],link->cp[cp].vt[vt]);
   linklib_vector_normalize(&in);

   out = linklib_vminus(link->cp[cp].vt[vt],link->cp[cp].vt[vt-1]);
   linklib_vector_normalize(&out);

   linklib_vweighted(tan,0.5,in,out);
   linklib_vector_normalize(&tan);
   
   return tan;

}


double linklib_link_length(linklib_link *L,double *component_lengths)

/* Procedure computes the length of each component of the link,
   and fills in the array of doubles "component_lengths", which 
   must be as long as L->nc. It returns the total length. We assume
   that fix_wrap has been called. */

{
  double tot_length;
  int cmp, nv, vert;
  linklib_pline *cp;

  tot_length = 0;
  for (cmp = 0; cmp < L->nc; cmp++) {

    component_lengths[cmp] = 0;
    cp = &L->cp[cmp];
    nv = (cp->acyclic) ? cp->nv-1 : cp->nv;

    for (vert = 0; vert < nv; vert++) {

      component_lengths[cmp] += linklib_vdist(cp->vt[vert+1],cp->vt[vert]);

    }

    tot_length += component_lengths[cmp];
  }

  return tot_length;
}

double linklib_link_parameter(linklib_link *L,int cmp,int vertnum)

/* Procedure reports the arclength distance from the given vertex */
/* to the 0th vertex of the given component of L. */

{
  double tot_length;
  int vert,nv;
  linklib_pline *cp;
  linklib_vector temp_vect;

  assert(L != NULL);
  assert(0 <= cmp && cmp <= L->nc);
  assert(0 <= vertnum && vertnum <= L->cp[cmp].nv-1);

  tot_length = 0;
  cp = &L->cp[cmp];
  nv = (cp->acyclic) ? cp->nv-1 : cp->nv;

  for (vert = 0; vert < vertnum; vert++) {
    temp_vect = cp->vt[vert+1];
    linklib_vsub(temp_vect,cp->vt[vert]);
    tot_length += linklib_norm(temp_vect);
  }

  return tot_length;
}

void linklib_link_force_closed(linklib_link *link)

     /* Procedure closes all open components of link by distributing a small
	change of all vertices of each such component. It also changes the 
	"acyclic" flag, calls fix_wrap. We lose one vertex in this process. */

{
  int i, this_cp;
  linklib_vector diff;

  for (this_cp=0;this_cp < link->nc;this_cp++) {

    if (link->cp[this_cp].acyclic == TRUE) {  
      /* Isolate the open components. */

      diff = link->cp[this_cp].vt[link->cp[this_cp].nv-1];   
      /* Compute the error in closure */

      linklib_vsub(diff,link->cp[this_cp].vt[0]);

      for(i=0;i<link->cp[this_cp].nv;i++) {

	linklib_vlincombine(link->cp[this_cp].vt[i],1.0,
			    diff,-((double)(i))/(double)(link->cp[this_cp].nv-1),
       			    link->cp[this_cp].vt[i]);

      }

      /* We claim to have moved the last vertex on top of the first. */

      diff = linklib_vminus(link->cp[this_cp].vt[0],
			    link->cp[this_cp].vt[link->cp[this_cp].nv-1]);
      assert(linklib_norm(diff) < 1e-10);

      /* Thus we eliminate the last vertex. */

      link->cp[this_cp].nv--;
      link->cp[this_cp].acyclic = FALSE;

    }

  }

  linklib_link_fix_wrap(link);

}
