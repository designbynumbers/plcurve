/*
 *
 * Data structures and prototypes for linklib_links
 *
 *  $Id: plCurve.h,v 1.3 2005-07-01 01:56:33 cantarel Exp $
 *
 */

/* Copyright 2004 The University of Georgia */

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

#ifndef __LINKLIB_LINK_H
#define __LINKLIB_LINK_H

#if (__cplusplus || c_plusplus)
extern "C" {
#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif 

#ifdef HAVE_ASSERT_H
#include <assert.h>
#endif

#include <truefalse.h>
#include <vector.h>

int  linklib_error_num;
char linklib_error_str[80];

typedef struct linklib_color_type {
  double r;
  double g;
  double b;
  double alpha;
} linklib_color;

typedef struct linklib_pline_type {
  int             acyclic;   /* This is an "open" pline (with distinct ends) */
  int             nv;        /* Number of vertices */
  int             cc;        /* Color count (number of colors) */
  linklib_vector *vt;        /* Actual vertices */
  linklib_color  *clr;       /* Colors */
  /***** Need a way to specify constraints on endpoints here *****/
} linklib_pline;

typedef struct linklib_link_type {	
  int nc;			/* Number of components */
  linklib_pline *cp;            /* Components */
} linklib_link;

/* 
 * Prototypes for routines to deal with links.  More in-depth documentation is
 * available in link.c.
 *
 */

/* Build a new link (with associated plines) */
linklib_link *linklib_link_new(int components, 
                               const int *nv, 
                               const int *acyclic,
                               const int *cc);

/* Free the link (and plines) */
void          linklib_link_free(linklib_link *L);

/* Read link data from a file */
linklib_link *linklib_link_read(FILE *infile);

/* Write link data to a file */
int           linklib_link_write(FILE *outfile, const linklib_link *L);

/* Fix the "hidden vertices" for easy handling of closed components */
void          linklib_link_fix_wrap(const linklib_link *L);

/* Count the edges in a link (correctly handling open/closed) */
int           linklib_link_edges(const linklib_link *L);

/* Copy a link */
linklib_link *linklib_link_copy(const linklib_link *L);

/* Compute tangent vector */
linklib_vector linklib_link_tangent_vector(linklib_link *link,int cp, int vt);

/* Find the arclength of a link. */
double linklib_link_length(linklib_link *L,double *component_lengths);

/* Find the arclength position of a point on a link. */
double linklib_link_parameter(linklib_link *L,int cmp,int vertnum);

/* Force a linklib_link to close as gently as possible */
void linklib_link_force_closed(linklib_link *link);

#if (__cplusplus || c_plusplus)
};
#endif
#endif
