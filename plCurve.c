/*
 *  Routines to create, destroy, read and write links (and plines)
 * 
 *  $Id: plCurve.c,v 1.29 2006-02-05 04:18:54 ashted Exp $
 *
 */

/* Copyright 2004 The University of Georgia. */

/* This file is part of plCurve.
   
plCurve is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

plCurve is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with plCurve; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif
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
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#ifdef HAVE_ASSERT_H
#include <assert.h>
#endif

#include <plCurve.h>

/*
 * Set up a new pline.  Pl should point to an *ALREADY ALLOCATED* pline (but
 * with an unallocated space for vertices).  The number of vertices is given in
 * nv and open is set to TRUE or FALSE depending on whether the pline is
 * open or closed.                                        
 *
 * We allocate two extra vertices, at -1 and nv to make "wrap-around" much 
 * simpler.
 */
static void pline_new(plCurve_pline *Pl,int nv, int open, int cc) {
 
  if (nv < 1) {
    plCurve_error_num = PLCL_E_TOO_FEW_VERTS;
    sprintf(plCurve_error_str,
      "pline_new: Can't create a pline with %d vertices.\n",nv);
    return;
  }

  Pl->open = open;
  Pl->nv = nv;
  if ((Pl->vt = 
       (plcl_vector *)calloc((nv+2),sizeof(plcl_vector))) == NULL) {
    plCurve_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plCurve_error_str,
      "pline_new: Can't allocate space for %d vertices in pline_new.\n",nv);
    return;
  }
  Pl->vt++; /* so that Pl->vt[-1] is a valid space */

  Pl->cc = cc;
  if ((Pl->clr = (plCurve_color *)calloc(cc,sizeof(plCurve_color))) == NULL) {
    plCurve_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plCurve_error_str,
      "pline_new: Can't allocate space for %d colors in pline_new.\n",cc);
    return;
  }
}

/*
 * Procedure allocates memory for a new link. The number of components is given
 * by components. The number of vertices in each component shows up in the
 * buffer pointed to be nv.  The closed/open nature of each pline is given in
 * the array pointed to by open.                           
 *
 */
plCurve *plCurve_new(int components, const int *nv, 
                     const int *open, const int *cc, 
                     const int ncst, const plCurve_constraint *cst)
{
  plCurve *L;
  int i;

  /* First, we check to see that the input values are reasonable. */

  plCurve_error_num = plCurve_error_str[0] = 0;
  if (components < 1) {
    plCurve_error_num = PLCL_E_TOO_FEW_COMPS;
    sprintf(plCurve_error_str,
      "plCurve_new: Can't create a link with %d components.\n",components);
    return NULL;
  }

  if (nv == NULL) {
    plCurve_error_num = PLCL_E_NULL_PTR;
#ifdef HAVE_STRLCPY
    strlcpy(plCurve_error_str,"plCurve_new: nv is NULL.\n",
      sizeof(plCurve_error_str));
#else
    strlncpy(plCurve_error_str,"plCurve_new: nv is NULL.\n",
      sizeof(plCurve_error_str)-1);
    plCurve_error_str[sizeof(plCurve_error_str)-1] = '\0';
#endif
    return NULL;
  }
  if (open == NULL) {
    plCurve_error_num = PLCL_E_NULL_PTR;
#ifdef HAVE_STRLCPY
    strlcpy(plCurve_error_str,"plCurve_new: open is NULL.\n",
      sizeof(plCurve_error_str));
#else
    strlncpy(plCurve_error_str,"plCurve_new: open is NULL.\n",
      sizeof(plCurve_error_str)-1);
    plCurve_error_str[sizeof(plCurve_error_str)-1] = '\0';
#endif
    return NULL;
  }
  if (cc == NULL) {
    plCurve_error_num = PLCL_E_NULL_PTR;
#ifdef HAVE_STRLCPY
    strlcpy(plCurve_error_str,"plCurve_new: cc is NULL.\n",
      sizeof(plCurve_error_str));
#else
    strlncpy(plCurve_error_str,"plCurve_new: cc is NULL.\n",
      sizeof(plCurve_error_str)-1);
    plCurve_error_str[sizeof(plCurve_error_str)-1] = '\0';
#endif
    return NULL;
  }

  if (ncst < 0) {
    plCurve_error_num = PLCL_E_NEG_CST;
    sprintf(plCurve_error_str,
      "plCurve_new: Can't have %d constraints.\n",ncst);
    return NULL;
  }

  if (ncst > 0 && cst == NULL) {
    plCurve_error_num = PLCL_E_NULL_PTR;
    sprintf(plCurve_error_str,
      "plCurve_new: Called with ncst=%d but cst NULL.\n",ncst);
    return NULL;
  }

  /* Now we attempt to allocate space for these components. */
  
  if ((L = (plCurve *)malloc(sizeof(plCurve))) == NULL) {
    plCurve_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plCurve_error_str,
      "plCurve_new: Could not allocate space for plCurve.\n");
    return NULL;
  }
  L->nc = components;
  if ((L->cp = (plCurve_pline *)malloc(L->nc*sizeof(plCurve_pline))) == NULL) {
    plCurve_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plCurve_error_str,
      "plCurve_new: Can't allocate array of pline ptrs.\n");
    return NULL;
  }
  L->ncst = ncst;
  if (ncst > 0) {
    if ((L->cst = (plCurve_constraint *)
      malloc(L->ncst*sizeof(plCurve_constraint))) == NULL) {
      plCurve_error_num = PLCL_E_CANT_ALLOC;
      sprintf(plCurve_error_str,
        "plCurve_new: Can't allocate array of constraints.\n");
      return NULL;
    }
  }

  for (i = 0; i < L->nc; i++) {
    pline_new(&L->cp[i],nv[i],open[i],cc[i]);
  }
  for (i = 0; i < L->ncst; i++) {
    L->cst[i] = cst[i];
  }

  return L;
}

/*
 * Find the closest point on the given line to the given point.
 *
 * coef should hold 6 doulbes, which we will call a,b,c,d,e,f
 * the components of the vector point we will call x,y,z
 *
 * The line is given parametrically by (at + b, ct + d, et + f) and the
 * closest point then has
 *  
 *       a(x - b) + c(y - d) + e(z - f)
 *   t = ------------------------------
 *              a^2 + c^2 + e^2
 *
 */
static plcl_vector Closest_line_point(const plcl_vector point, 
                                      const double *coef) {

  plcl_vector ret_vect;
  double a = coef[0];
  double b = coef[1];
  double c = coef[2];
  double d = coef[3];
  double e = coef[4];
  double f = coef[5];
  double x = point.c[0];
  double y = point.c[1];
  double z = point.c[2];

  double t = (a*(x-b) + c*(y-d) + e*(z-f))/(a*a + c*c + e*e);

  ret_vect.c[0] = a*t + b;
  ret_vect.c[1] = c*t + d;
  ret_vect.c[2] = e*t + f;

  return ret_vect;
}

/*
 * Find the closest point on the given plane to the given point.
 *
 * coef should hold 4 doubles, which we will call a,b,c,d
 * the components of the vector point we will call x1,y1,z1
 *
 * The plane can be given by ax + by + cz = d or parameterized as
 *
 *          d - ax - by
 *   (x, y, -----------)
 *               c
 *
 * The closest point then has
 *
 *              a x1 + b y1 + c z1 - d
 *   x = x1 - a ----------------------
 *                 a^2 + b^2 + c^2
 *
 *              a x1 + b y1 + c z1 - d
 *   y = y1 - b ----------------------
 *                 a^2 + b^2 + c^2
 *
 */
static plcl_vector Closest_plane_point(const plcl_vector point, 
                                       const double *coef) {

  plcl_vector ret_vect;
  double a = coef[0];
  double b = coef[1];
  double c = coef[2];
  double d = coef[3];
  double x1 = point.c[0];
  double y1 = point.c[1];
  double z1 = point.c[2];

  double frac = (a*x1 + b*y1 + c*z1 - d)/(a*a + b*b + c*c);
  double x = x1 - a*frac;
  double y = y1 - b*frac;

  ret_vect.c[0] = x;
  ret_vect.c[1] = y;
  ret_vect.c[2] = (d - a*x - b*y)/c;

  return ret_vect;
}

/*
 * Check to see if the given constraint is satisfied and return value for 
 * how far off it is.
 *
 */

double plCurve_cst_check(const plCurve L, const plCurve_constraint cst) {
  
  plcl_vector closest;

  plCurve_error_num = plCurve_error_str[0] = 0;
  if (cst.kind == PLCL_FIXED) {
    closest.c[0] = cst.coef[0];
    closest.c[1] = cst.coef[1];
    closest.c[2] = cst.coef[2];
  } else if (cst.kind == PLCL_ON_LINE) {
    closest = Closest_line_point(L.cp[cst.cp].vt[cst.vt],cst.coef);
  } else if (cst.kind == PLCL_IN_PLANE) {
    closest = Closest_plane_point(L.cp[cst.cp].vt[cst.vt],cst.coef);
  } else {
    plCurve_error_num = PLCL_E_BAD_CST_KIND;
    sprintf(plCurve_error_str,
      "plCurve_cst_check: Unknown constraint kind: %d.\n",cst.kind);
    return -1;
  }

  return linklib_vdist(L.cp[cst.cp].vt[cst.vt],closest);
}

/*
 * Move the vertex to the closest possible point which satisfies the
 * constraint.
 *
 */

double plCurve_cst_fix(const plCurve L, const plCurve_constraint cst) {
  
  plcl_vector closest;

  plCurve_error_num = plCurve_error_str[0] = 0;
  if (cst.kind == PLCL_FIXED) {
    closest.c[0] = cst.coef[0];
    closest.c[1] = cst.coef[1];
    closest.c[2] = cst.coef[2];
  } else if (cst.kind == PLCL_ON_LINE) {
    closest = Closest_line_point(L.cp[cst.cp].vt[cst.vt],cst.coef);
  } else if (cst.kind == PLCL_IN_PLANE) {
    closest = Closest_plane_point(L.cp[cst.cp].vt[cst.vt],cst.coef);
  } else {
    plCurve_error_num = PLCL_E_BAD_CST_KIND;
    sprintf(plCurve_error_str,
      "plCurve_cst_fix: Unknown constraint kind: %d.\n",cst.kind);
    return -1;
  }

  double dist = linklib_vdist(L.cp[cst.cp].vt[cst.vt],closest);
  L.cp[cst.cp].vt[cst.vt] = closest;

  return dist;
}

/*
 * Free the memory used to hold vertices in a given pline (not the memory of
 * the pline itself).  We then set all the values in the link data structure to
 * reflect the fact that the memory has been freed.  We can call pline_free
 * twice on the same pline without fear. 
 *
 */ 
static void pline_free(plCurve_pline *Pl) {
  
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
void plCurve_free(plCurve *L) {
  int i;

  plCurve_error_num = plCurve_error_str[0] = 0;

  /* First, we check the input. */
  if (L == NULL) {
    return; /* Move along, nothing to see here */
  }

  if (L->nc < 0) {
    plCurve_error_num = PLCL_E_TOO_FEW_COMPS;
    sprintf(plCurve_error_str,
      "plCurve_free: plCurve appears corrupted. L.nc = %d.\n",L->nc);
    return;
  }

  /* Now we can get to work. */
  for (i=0; i<L->nc; i++) {
    pline_free(&L->cp[i]);
  }

  free(L->cp);
  L->nc = 0;
  if (L->cst != NULL) {
    free(L->cst);
  }
  L->ncst = 0;

  free(L);
  L = NULL;
} /* plCurve_free */

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

int plCurve_write(FILE *file, const plCurve *L) {
  int i,j;              /* Counter for the for loops */ 
  int nverts = 0;       /* Total number of vertices of all components */
  int colors = 0;       /* Total number of colors of all components */

  plCurve_error_num = plCurve_error_str[0] = 0;

  /* First, do a little sanity checking. */
  if (L == NULL) {
    plCurve_error_num = PLCL_E_NULL_PTR;
    sprintf(plCurve_error_str,
      "plCurve_write: Passed NULL pointer as plCurve.\n");
    return -1;
  }

  if (file == NULL) {
    plCurve_error_num = PLCL_E_NULL_PTR;
    sprintf(plCurve_error_str,"plCurve_write: Passed NULL pointer as file.\n");
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
    if (L->cp[i].open) {
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
static int skip_whitespace_and_comments(FILE *infile)
{
  int thischar,commentflag = {FALSE};

  /* First, we check to make sure that infile looks legit. */
  if (infile == NULL) {
    plCurve_error_num = PLCL_E_NULL_PTR;
    sprintf(plCurve_error_str,
      "skip_whitespace_and_comments: infile is a null pointer.\n");
    return -1;
  }
  
  /* Now we start to work. */
  for(;;) {
    thischar = fgetc(infile);

    if (thischar == EOF) {  
      /* Reached end of file before a non-space, non-comment */
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
static int scandoubles(FILE *infile,int ndoubles, ...)
{
  int nconverted = 0,i;
  va_list ap;
  double *thisdouble;

  /* First, we check for overall sanity. */

  if (infile == NULL) {
    plCurve_error_num = PLCL_E_NULL_PTR;
    sprintf(plCurve_error_str,"scandoubles: infile is a null pointer.\n");
    return -1;
  }

  if (ndoubles < 1) {
    plCurve_error_num = PLCL_E_TOO_FEW_DBLS;
    sprintf(plCurve_error_str,
      "scandoubles: ndoubles (%d) is less than one.\n", ndoubles);
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
static int scanints(FILE *infile,int nints, ...)
{
  int nconverted = 0,i;
  va_list ap;
  int *thisint;

  /* First, we check for overall sanity. */

  if (infile == NULL) {
    plCurve_error_num = PLCL_E_NULL_PTR;
    sprintf(plCurve_error_str,"scanints: infile is a null pointer.\n");
    return -1;
  }

  if (nints < 1) {
    plCurve_error_num = PLCL_E_TOO_FEW_INTS;
    sprintf(plCurve_error_str,"scanints: nints (%d) is less than one.\n",nints);
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
static int skipwhitespace(FILE *infile)
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
void plCurve_fix_wrap(const plCurve *L) {
  int i,nv;

  plCurve_error_num = plCurve_error_str[0] = 0;

  for (i = 0; i < L->nc; i++) {
    nv = L->cp[i].nv;
    if (L->cp[i].open) {
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
plCurve *plCurve_read(FILE *file) 
{
  plCurve *L;
  int nverts, ncomp, ncolors;
  int *nvarray, *open, *ccarray;
  int i, j;
  int nv;
  
  plCurve_error_num = plCurve_error_str[0] = 0;

  /* First, we check for the 'VECT' keyword. */
  if (fscanf(file," VECT ") == EOF) {
  
    plCurve_error_num = PLCL_E_NO_VECT;
    sprintf(plCurve_error_str,"plCurve_read: Couldn't find VECT keyword.\n");
    return NULL;
  }

  /* Now we read the three integers giving vertices, components, and colors. */

  if (scanints(file,3,&ncomp,&nverts,&ncolors) != 3) {
  
    plCurve_error_num = PLCL_E_BAD_CVC_LINE;
    sprintf(plCurve_error_str,
      "plCurve_read: Couldn't parse <ncomp> <nverts> <ncolors> line.\n");
    return NULL;
  }

  /* We now try to read the array of numbers of vertices. */

  nvarray = (int *)calloc(ncomp,sizeof(int));
  open    = (int *)calloc(ncomp,sizeof(int));
  ccarray = (int *)calloc(ncomp,sizeof(int));

  for(i=0;i<ncomp;i++) {
    if (scanints(file,1,&(nvarray[i])) != 1) {
      plCurve_error_num = PLCL_E_BAD_CVRT_LINE;
      sprintf(plCurve_error_str,"plCurve_read: Couldn't parse number"
              "of vertices in component %d.\n",i);    
      return NULL;
    }
    if (nvarray[i] < 0) {
      /* A negative number of vertices indicates a CLOSED component. */
      open[i] = FALSE;  
      nvarray[i] *= -1;
    } else {
      open[i] = TRUE;
    }
  }

  /* We have set nvarray and open and are ready to read the color data.  */

  for(i=0;i<ncomp;i++) {
    if (scanints(file,1,&(ccarray[i])) != 1) {
      plCurve_error_num = PLCL_E_BAD_CLR_LINE;
      sprintf(plCurve_error_str,"plCurve_read: Couldn't parse <ncolors>"
        "for component %d.\n", i);   
      return NULL;
    }
  }

  /* We now allocate the link data structure. */

  L = plCurve_new(ncomp,nvarray,open,ccarray,0,NULL);

  if (L == NULL) {   /* If we don't have this much memory, then return NULL. */
    sprintf(plCurve_error_str,"plCurve_read: Error in plCurve_new: %d.\n",
      plCurve_error_num);
    plCurve_error_num = PLCL_E_NEW_FAILED;
    return NULL;
  }

  /* And get ready to read the actual data. */

  for(i = 0; i < ncomp; i++) {
    nv = L->cp[i].nv;
    for(j = 0; j < nv; j++) {
      if (scandoubles(file,3,&L->cp[i].vt[j].c[0],
                             &L->cp[i].vt[j].c[1],
                             &L->cp[i].vt[j].c[2]) != 3) {
        plCurve_free(L);
        plCurve_error_num = PLCL_E_BAD_VERT_LINE;
        sprintf(plCurve_error_str,"plCurve_read: Couldn't parse "
          " <x> <y> <z> data for vertex %d of component %d.\n",j,i);
        return NULL;
      }
    }
  }
  /* Now set the "wrap-around" vertices */
  plCurve_fix_wrap(L);

    /* And finally the colors. Unfortunately, to really comply with 
        the Geomview standard here we have to be kind of careful. */

  /* And finally the colors. */
  for (i=0; i < ncomp; i++) {
    for (j=0; j < L->cp[i].cc; j++) {
      if (scandoubles(file,4, &L->cp[i].clr[j].r, &L->cp[i].clr[j].g,
           &L->cp[i].clr[j].b, &L->cp[i].clr[j].alpha) != 4) {

        plCurve_error_num = PLCL_E_BAD_COLOR;
        sprintf(plCurve_error_str,"plCurve_read: Couldn't parse color %d "
          "in component %d of link.\n",j,i);
        return NULL;

      }
    }
  }

  free(ccarray);
  free(open);
  free(nvarray);

  return L;
}

#define pline_edges(P) (((P).open) ? (P).nv-1 : (P).nv)
/* 
 *   Return the total number of edges in link. 
 */
int plCurve_edges(const plCurve *L) 
{
  int i, edges = 0;

  plCurve_error_num = plCurve_error_str[0] = 0;

  for (i=0;i<L->nc;i++) {
    edges += pline_edges(L->cp[i]);
  }
  return edges;
}

/* Compute the curvature of L at vertex vt of component cp */

double plCurve_curvature(const plCurve *L, 
                              const int comp, 
                              const int vert)
{
  double         kappa;
  plcl_vector in,out;
  double normin, normout;
  double dot_prod,cross_prod_norm;

  /* We start with some initializations. */
  
  plCurve_error_num = plCurve_error_str[0] = 0;
  plCurve_fix_wrap(L);
  
  /* Now we work. */

  in = linklib_vminus(L->cp[comp].vt[vert],L->cp[comp].vt[vert-1]);
  normin = linklib_norm(in);
  out = linklib_vminus(L->cp[comp].vt[vert+1],L->cp[comp].vt[vert]);
  normout = linklib_norm(out);
      
  dot_prod = linklib_dot(in,out);
  cross_prod_norm = linklib_norm(linklib_cross(in,out));

  if (normin*normout + dot_prod < 1e-12) {
    plCurve_error_num = PLCL_E_INF_KAPPA;
    sprintf(plCurve_error_str,"plCurve_curvature: kappa not finite "
      "at vertex %d of component %d.\n",comp,vert);
    return -1.0;
  }
      
  kappa = (2*cross_prod_norm)/(normin*normout + dot_prod);
  kappa /= (normin < normout) ? normin : normout;

  return kappa;
}

/* 
 * Duplicate a link and return the duplicate 
 *
 */
plCurve *plCurve_copy(const plCurve *L) {
  plCurve *nL;
  int *nv,*open,*ccarray;
  int cnt,cnt2;

  plCurve_error_num = plCurve_error_str[0] = 0;

  if ((nv = (int *)malloc((L->nc)*sizeof(int))) == NULL ||
      (open = (int *)malloc((L->nc)*sizeof(int))) == NULL ||
      (ccarray = (int *)malloc((L->nc)*sizeof(int))) == NULL) {
    plCurve_error_num = PLCL_E_CANT_ALLOC;
    sprintf(plCurve_error_str,"Unable to malloc space for alternate link.\n");
    return NULL;
  }
  for (cnt = 0; cnt < L->nc; cnt++) {
    nv[cnt] = L->cp[cnt].nv;
    open[cnt] = L->cp[cnt].open;
    ccarray[cnt] = L->cp[cnt].cc;
  }
  nL = plCurve_new(L->nc,nv,open,ccarray,0,NULL);

  for (cnt = 0; cnt < L->nc; cnt++) {
    /*
     * Copy the vertices (including the "hidden" ones, so we don't have to call
     * plCurve_fix_wrap).
     */
    for (cnt2 = -1; cnt2 <= L->cp[cnt].nv; cnt2++) { 
      nL->cp[cnt].vt[cnt2] = L->cp[cnt].vt[cnt2];
    }
    for (cnt2 = 0; cnt2 < L->cp[cnt].cc; cnt2++) {
      nL->cp[cnt].clr[cnt2] = L->cp[cnt].clr[cnt2];
    }
  }

  free(ccarray);
  free(open);
  free(nv);

  return nL;
}

plcl_vector plCurve_tangent_vector(plCurve *link,int cp, int vt)

/* Procedure computes a (unit) tangent vector to <link> 
   at the given vertex of the given component. */

{
  plcl_vector in, out, tan;

  plCurve_error_num = plCurve_error_str[0] = 0;

  if (link->cp[cp].open) {

    if (vt == 0) {

      tan = linklib_vminus(link->cp[cp].vt[1],
                  link->cp[cp].vt[0]);

      plcl_vector_normalize(&tan);

      return tan;

    } else if (vt == link->cp[cp].nv-1) {

       tan = linklib_vminus(link->cp[cp].vt[link->cp[cp].nv-1],
                   link->cp[cp].vt[link->cp[cp].nv-2]);

       plcl_vector_normalize(&tan);

       return tan;

    }

  }

  /* We now know that either we are on a closed 
     component, or we are not at an endpoint.   */
  
   in = linklib_vminus(link->cp[cp].vt[vt+1],link->cp[cp].vt[vt]);
   plcl_vector_normalize(&in);

   out = linklib_vminus(link->cp[cp].vt[vt],link->cp[cp].vt[vt-1]);
   plcl_vector_normalize(&out);

   linklib_vweighted(tan,0.5,in,out);
   plcl_vector_normalize(&tan);
   
   return tan;

}


double plCurve_length(plCurve *L,double *component_lengths)

/* Procedure computes the length of each component of the link,
   and fills in the array of doubles "component_lengths", which 
   must be as long as L->nc. It returns the total length. We assume
   that fix_wrap has been called. */

{
  double tot_length;
  int cmp, nv, vert;
  plCurve_pline *cp;

  plCurve_error_num = plCurve_error_str[0] = 0;

  tot_length = 0;
  for (cmp = 0; cmp < L->nc; cmp++) {

    component_lengths[cmp] = 0;
    cp = &L->cp[cmp];
    nv = (cp->open) ? cp->nv-1 : cp->nv;

    for (vert = 0; vert < nv; vert++) {

      component_lengths[cmp] += linklib_vdist(cp->vt[vert+1],cp->vt[vert]);

    }

    tot_length += component_lengths[cmp];
  }

  return tot_length;
}


/* Procedure reports the arclength distance from the given vertex */
/* to the 0th vertex of the given component of L. */
double plCurve_parameter(plCurve *L,int cmp,int vertnum) {

  double tot_length;
  int vert,nv;
  plCurve_pline *cp;
  plcl_vector temp_vect;

  plCurve_error_num = plCurve_error_str[0] = 0;

  if (L == NULL) {
    plCurve_error_num = PLCL_E_NULL_PTR;
    sprintf(plCurve_error_str,
      "plCurve_parameter: Passed NULL pointer as plCurve.\n");
    return -1;
  }

  if (cmp < 0 || cmp >= L->nc) {
    plCurve_error_num = PLCL_E_BAD_COMPONENT;
    sprintf(plCurve_error_str,
      "plCurve_parameter: Component value out of range (0..%d): %d.\n",
      L->nc-1, cmp);
    return -1;
  }

  if (vertnum < 0 || vertnum >= L->cp[cmp].nv) {
    plCurve_error_num = PLCL_E_BAD_VERTEX;
    sprintf(plCurve_error_str,
      "plCurve_parameter: Vertex value out of range (0..%d): %d.\n",
      L->cp[cmp].nv-1, vertnum);
    return -1;
  }

  tot_length = 0;
  cp = &L->cp[cmp];
  nv = (cp->open) ? cp->nv-1 : cp->nv;

  for (vert = 0; vert < vertnum; vert++) {
    temp_vect = cp->vt[vert+1];
    linklib_vsub(temp_vect,cp->vt[vert]);
    tot_length += linklib_norm(temp_vect);
  }

  return tot_length;
}

/*
 * This procedure closes all open components of link by distributing a small
 * change of all vertices of each such component. It also changes the "open"
 * flag and calls fix_wrap. We lose one vertex in this process. 
 */
void plCurve_force_closed(plCurve *link)
{
  int i, cmp;
  plcl_vector diff;

  plCurve_error_num = plCurve_error_str[0] = 0;

  for (cmp=0;cmp < link->nc;cmp++) {
    if (link->cp[cmp].open == TRUE) {  /* Isolate the open components. */

      /* Compute the error in closure */
      diff = link->cp[cmp].vt[link->cp[cmp].nv-1];   
      linklib_vsub(diff,link->cp[cmp].vt[0]);

      for(i=0;i<link->cp[cmp].nv;i++) {
        linklib_vlincombine(link->cp[cmp].vt[i],1.0,
                            diff,-((double)(i))/(double)(link->cp[cmp].nv-1),
                            link->cp[cmp].vt[i]);
      }

      /* We claim to have moved the last vertex on top of the first. */
      diff = linklib_vminus(link->cp[cmp].vt[0],
                            link->cp[cmp].vt[link->cp[cmp].nv-1]);
      assert(linklib_norm(diff) < 1e-10);

      /* Thus we eliminate the last vertex. */
      link->cp[cmp].nv--;
      link->cp[cmp].open = FALSE;
    }
  }

  plCurve_fix_wrap(link);
}
