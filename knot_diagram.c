#include <config.h>
#include <plCurve.h>

#ifdef HAVE_STDIO_H
  #include <stdio.h>
#endif
#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif
#ifdef HAVE_FLOAT_H
  #include <float.h>
#endif
#ifdef HAVE_MATH_H
  #include <math.h>
#endif
#ifdef HAVE_ASSERT_H
  #include <assert.h>
#endif
#ifdef HAVE_MALLOC_H
  #include <malloc.h>
#endif
#ifdef HAVE_STRING_H
  #include <string.h>
#endif
#ifdef HAVE_UNISTD_H
  #include <unistd.h>
#endif
#ifdef HAVE_TIME_H
  #include <time.h>
#endif
#ifdef HAVE_ARGTABLE2_H
  #include <argtable2.h>
#endif

/* Remove a vertex from a plCurve by breaking components into pieces */
static void remove_vertex(plCurve * const L, const int cmp, const int vert,
                          bool *start_over, plc_vector N) {
  plc_vector *vt;
  plc_color *clr;

  assert(L != NULL);
  assert(cmp >= 0);
  if (cmp >= L->nc) {
    /* For some reason, we are asking to deal with too large a component.
       Print data so that the problem can be reproduced */
    fprintf(stderr,"Error: %d = cmp >= L->nc = %d",cmp,L->nc);
    fprintf(stderr,"  Angle: %g,%g,%g\n",plc_M_clist(N));
  }
  assert(cmp < L->nc);
  assert(vert >= 0);
  assert(vert < L->cp[cmp].nv);

  *start_over = false;
  
  if (vert == L->cp[cmp].nv-1 || 
      (vert == L->cp[cmp].nv-2 && L->cp[cmp].open)) { 
    /* ultimate (or penultimate on an open component) vertex, simply trim */
    L->cp[cmp].nv = vert;
    if (L->cp[cmp].cc > L->cp[cmp].nv) { 
      L->cp[cmp].cc = L->cp[cmp].nv;
    }
    L->cp[cmp].open = true;
  } else if (vert == 0 || (vert == 1 && L->cp[cmp].open)) {
    /* similarly for first or second-on-open vertex */
    L->cp[cmp].nv -= (vert + 1);
    memmove(L->cp[cmp].vt, &(L->cp[cmp].vt[vert+1]),
        L->cp[cmp].nv*sizeof(plc_vector));
    if (L->cp[cmp].cc > L->cp[cmp].nv) { 
      L->cp[cmp].cc = L->cp[cmp].nv;
      memmove(L->cp[cmp].clr, &(L->cp[cmp].clr[vert+1]),
          L->cp[cmp].nv*sizeof(plc_color));
    }
    L->cp[cmp].open = true;
  } else {
    /* Middle vertex */
    if (L->cp[cmp].open) {
      plc_add_component(L,cmp+1,L->cp[cmp].nv-(vert+1),true,
          (L->cp[cmp].cc < L->cp[cmp].nv) ? L->cp[cmp].cc :
          L->cp[cmp].nv-(vert+1),&(L->cp[cmp].vt[vert+1]),
          (L->cp[cmp].cc < L->cp[cmp].nv) ? L->cp[cmp].clr :
          &(L->cp[cmp].clr[vert+1]));
      L->cp[cmp].nv = vert;
      if (L->cp[cmp].cc > L->cp[cmp].nv) { 
        L->cp[cmp].cc = L->cp[cmp].nv;
      }
    } else {
      /* We're eliminating a middle vertex on a closed component.  We need to
         open the component, but do so *at that point*.  This will require
         moving the vertex after that one to the begining and the vertex before
         that one to the end */
      *start_over = true;
      L->cp[cmp].open = true;
      L->cp[cmp].nv--;
      vt = malloc(L->cp[cmp].nv*sizeof(plc_vector));
      assert(vt != NULL);
      memcpy(vt,&(L->cp[cmp].vt[vert+1]),
          (L->cp[cmp].nv-vert)*sizeof(plc_vector));
      memcpy(&(vt[L->cp[cmp].nv-vert]),L->cp[cmp].vt,vert*sizeof(plc_vector));
      memcpy(L->cp[cmp].vt,vt,L->cp[cmp].nv*sizeof(plc_vector));
      if (L->cp[cmp].cc > L->cp[cmp].nv) {
        clr = malloc(L->cp[cmp].nv*sizeof(plc_color));
        assert(clr != NULL);
        L->cp[cmp].cc = L->cp[cmp].nv;
        memcpy(clr,&(L->cp[cmp].clr[vert+1]),
            (L->cp[cmp].nv-vert)*sizeof(plc_color));
        memcpy(&(clr[L->cp[cmp].nv-vert]),L->cp[cmp].clr,
            vert*sizeof(plc_color));
        memcpy(L->cp[cmp].clr,clr,L->cp[cmp].nv*sizeof(plc_color));
        free(clr);
      }
      free(vt);
      /*@-branchstate@*/
    }
    /*@=branchstate@*/
  }
  if (L->cp[cmp].nv <= 1) {
    /* Ran out of edges here, get rid of the component */
    plc_drop_component(L,cmp);
    *start_over = true;
  }
  plc_fix_wrap(L);
}

#define pi 3.1415926
#define showCurve(F) \
  if (show_work >= 5) { \
    fprintf(geomview,"(geometry Curve "); \
    plc_write(geomview,F); \
    fprintf(geomview,") (look-recenter Curve) (look-encompass Curve)\n"); \
    fflush(NULL); \
  }

/* Project the vertices onto the plane to which N is normal (returning them
 * as (x,y,0) vectors in a plCurve framework. */
static plCurve *flatten(const plCurve * const L, 
                        const plc_vector N, 
                        const double gap,
                        const double neighbor_gap,
                              FILE *geomview,
                        const int delay,
                        const int show_work) {
  plCurve *F,*H;  /* Flattened knot, heights of flattened points */
  int cmp,cmp2,vert,vert2;
  plc_vector n;
  plc_vector *vt;
  plc_vector e3;
  double cos_theta,sin_theta;
  plc_vector axle,u_para,u_perp,v_cross_u;
  bool ok,far_enough;
  double gap_sq = gap*gap;
  double n_gap_sq = neighbor_gap*neighbor_gap;
  double dist;
  int last_to_check;
  bool start_over;
  int wherecnt;

  assert(L != NULL);

  F = plc_copy(L);

  /* First rotate so that the given normal vector lies along the z axis. */
  n = plc_normalize_vect(N,NULL); /* Bail out if we give it a bad vector */
  e3 = plc_build_vect(0.0,0.0,1.0);
  /* What should we rotate around? */
  axle = plc_cross_prod(n,e3);
  ok = true;
  axle = plc_normalize_vect(axle,&ok);
  if (ok) { 
    /* We actually do need to do some rotation */
    for (cmp = 0; cmp < F->nc; cmp++) {
      for (vert = 0; vert < F->cp[cmp].nv; vert++) {
        vt = &F->cp[cmp].vt[vert];
        /* Rotate space so that n is pointing along the z axis and just take 
           the first two components.  Formulas from
             http://en.wikipedia.org/wiki/Coordinate_rotation#Three_dimensions
         */
        cos_theta = plc_dot_prod(n,e3);
        sin_theta = sqrt(1-cos_theta*cos_theta);
        u_para = plc_scale_vect(plc_dot_prod(*vt,axle),axle);
        u_perp = plc_vect_diff(*vt,u_para);
        v_cross_u = plc_cross_prod(*vt,axle);
        *vt = plc_vect_sum(u_para,
            plc_vlincomb(cos_theta,u_perp,sin_theta,v_cross_u));
      }
    }
  }
  showCurve(F);
  if (show_work >=5) { usleep(delay*100); }

  H = plc_copy(F);
  assert(H != NULL);
  /* Flatten F.  Retain heights in H. */
  for (cmp = 0; cmp < F->nc; cmp++) {
    for (vert = 0; vert < F->cp[cmp].nv; vert++) {
      F->cp[cmp].vt[vert].c[2] = 0.0;
    }
  }
  showCurve(F);
  if (show_work >= 5) { usleep(delay*100); }

  /* Now eliminate vertices to show crossings. */
  for (cmp = 0; cmp < F->nc; cmp++) {
    for (vert = 0; vert < F->cp[cmp].nv; vert++) {
      if (show_work >= 10) {
        fprintf(geomview,"(geometry Where VECT 1 16 1 -16 1 ");
        for (wherecnt = 0; wherecnt < 16; wherecnt++) {
          fprintf(geomview,"%g %g 0.0 ",
              F->cp[cmp].vt[vert].c[0]+gap*cos(wherecnt*2*pi/16),
              F->cp[cmp].vt[vert].c[1]+gap*sin(wherecnt*2*pi/16));
        }
        fprintf(geomview,"1 0 0 1)\n");
        fflush(NULL);
        if (show_work >= 15) {
          usleep(delay);
        }
      }
      far_enough = false;
      vt = &(F->cp[cmp].vt[vert]);
      last_to_check = F->cp[cmp].nv-1;
      if (!F->cp[cmp].open && plc_sq_dist(*vt,F->cp[cmp].vt[0]) <= n_gap_sq) {
        /* We are close to the beginning of a closed strand, perhaps we
         * should avoid checking some of the end vertices. */
        while (last_to_check > vert &&
            plc_sq_dist(*vt,F->cp[cmp].vt[last_to_check]) <= n_gap_sq) {
          last_to_check--;
        }
      }
      for (cmp2 = (vert == F->cp[cmp].nv-1) ? cmp+1 : cmp;
           cmp2 < F->nc; cmp2++) {
        /* Another component is always "far enough away" to form a crossing. */
        far_enough = far_enough || (cmp != cmp2);
        for (vert2 = ((cmp2 == cmp) ? vert+1 : 0);
             vert2 < ((cmp2 == cmp) ? last_to_check+1 : F->cp[cmp2].nv); 
             vert2++) {
          if (show_work >= 20) {
            fprintf(geomview,"(geometry Versus VECT 2 4 2 2 2 1 1 "
                "%g %g 0 %g %g 0 %g %g 0 %g %g 0 0 0 1 1 0 0 1 1)\n",
                F->cp[cmp2].vt[vert2].c[0]-gap/2.0,
                F->cp[cmp2].vt[vert2].c[1]-gap/2.0,
                F->cp[cmp2].vt[vert2].c[0]+gap/2.0,
                F->cp[cmp2].vt[vert2].c[1]+gap/2.0,
                F->cp[cmp2].vt[vert2].c[0]-gap/2.0,
                F->cp[cmp2].vt[vert2].c[1]+gap/2.0,
                F->cp[cmp2].vt[vert2].c[0]+gap/2.0,
                F->cp[cmp2].vt[vert2].c[1]-gap/2.0);
            fflush(NULL);
            if (show_work >= 25) {
              usleep(delay);
            }
          }
          /* Wander along, looking for any time the other point comes within
           * gap of cmp:vert (after it first gets at least gap away from it).
           * That marks a crossing and we remove one of the vertices in
           * question -- the lower one. */
          dist = plc_sq_dist(F->cp[cmp].vt[vert],F->cp[cmp2].vt[vert2]);
          far_enough = far_enough || (dist >= n_gap_sq);
          if ((dist <= gap_sq) && far_enough) {
            /* Found a crossing */
            start_over = false;
            if (H->cp[cmp].vt[vert].c[2] >= H->cp[cmp2].vt[vert2].c[2]) {
              remove_vertex(F,cmp2,vert2,&start_over,N);
              remove_vertex(H,cmp2,vert2,&start_over,N);
              showCurve(F);
              vert2--;
              if (cmp2 == cmp) {
                /* This is now an open strand, check all vertices */
                last_to_check = F->cp[cmp].nv-1;
              }
              if (start_over) {
                if ((cmp2 != cmp) && (cmp2 < F->nc)) {
                  far_enough = true;
                  vert2 = -1;
                } else {
                  /* Could be we just renumbered every vertex in this
                   * component.  We'd better start clear over from scratch. */
                  far_enough = false;
                  vert = -1;
                  cmp2 = F->nc-1;
                  vert2 = F->cp[cmp2].nv-1;
                }
              }
            } else {
              remove_vertex(F,cmp,vert,&start_over,N);
              remove_vertex(H,cmp,vert,&start_over,N);
              showCurve(F);
              if (start_over) {
                /* Run this component again */
                vert = -1;
              } else {
                /* Just run this vertex again */
                vert--;
              }
              cmp2 = F->nc-1;
              vert2 = F->cp[cmp2].nv-1;
              if (vert2 < last_to_check) { vert2 = last_to_check; }
            }
          }
        }
      }
    }
  }
            
  plc_free(H);
  return F;
}

/* Add a component which contains the 2d convex hull of the other components
 * (looking only at the first 2 dimensions, of course).  We'll be using the
 * gift-wrapping algorithm. */
static void add_convex_hull(plCurve *L) {
  plc_vector *vt;
  int nv,verts;
  int cmp,vert;
  plc_vector first,next;
  plc_vector x_vect,y_vect,cross;
  plc_color *clr;

  for (nv = 0, cmp = 0; cmp < L->nc; cmp++) {
    nv += L->cp[cmp].nv;
  }
  vt = malloc(nv*sizeof(plc_vector));
  assert(vt != NULL);
  clr = malloc(nv*sizeof(plc_vector));
  assert(clr != NULL);
  verts = 0;

  /* Look for the lowest point */
  first = L->cp[0].vt[0];
  for (cmp = 0; cmp < L->nc; cmp++) {
    for (vert = 0; vert < L->cp[cmp].nv; vert++) {
      if (first.c[1] - L->cp[cmp].vt[vert].c[1] > DBL_EPSILON) {
        first = L->cp[cmp].vt[vert];
      }
    }
  }
  vt[verts++] = first;

  /* Area of the triangle made by three vertices.  Formula adapted from 
       http://www.cse.unsw.edu.au/~lambert/java/3d/triangle.html
     which says that the area of the triangle defined by 
       p1 = (x1,y1), p2 = (x2,y2), p3 = (x3,y3) is

          1  | x1 x2 x3 |
         --- | y1 y2 y3 |
          2  |  1  1  1 |

       which looks to me a lot like

         1
         - (x1,x2,x3) x (y1,y2,y3)
         2 
   */

  next = plc_build_vect(0.0,0.0,0.0);
  /*@-compdef@*/
  while (!plc_vecteq(vt[0],next) && verts < nv) {
  /*@=compdef@*/
    /* There is one other vertex out there which, when it is the other end of a
     * line segment with "first" all of the other vectors are to the left of
     * the line.  Find it in as pedantic a way as possible, ye olde brute
     * search (and a particularly stupid one at that). */
    next = (plc_vecteq(L->cp[0].vt[0],first)) ?
      L->cp[0].vt[1] : L->cp[0].vt[0];
    x_vect.c[0] = first.c[0];
    y_vect.c[0] = first.c[1];
    x_vect.c[1] = next.c[0];
    y_vect.c[1] = next.c[1];
    for (cmp = 0; cmp < L->nc; cmp++) {
      for (vert = 0; vert < L->cp[cmp].nv; vert++) {
        if (!plc_vecteq(L->cp[cmp].vt[vert],first) &&
            !plc_vecteq(L->cp[cmp].vt[vert],next)) {
          x_vect.c[2] = L->cp[cmp].vt[vert].c[0];
          y_vect.c[2] = L->cp[cmp].vt[vert].c[1];
          cross = plc_cross_prod(x_vect,y_vect);
          if (- cross.c[0] - cross.c[1] - cross.c[2] > DBL_EPSILON) {
            /* cmp:vert is to the right of (or same as) next, take it instead */
            next = L->cp[cmp].vt[vert];
            x_vect.c[1] = x_vect.c[2];
            y_vect.c[1] = y_vect.c[2];
          }
        }
      }
    }
    vt[verts++] = next;
    first = next;
  } /* while */
  verts--;
  /*@+loopexec@*/
  for (vert = 0; vert < verts; vert++) {
    clr[vert].r = 1.0 - ((double)vert/(double)verts);
    clr[vert].g = 0.0;
    clr[vert].b = (double)vert/(double)verts;
    clr[vert].alpha = 1.0;
  }
  /*@=loopexec@*/
  plc_add_component(L,L->nc,verts,false,verts,vt,clr);
  free(vt);
  free(clr);
}

/* Rotate (in 2-dimensions) the first n-1 components based on the convex hull
 * in the nth component, so that the shortest "diameter" lies in the y
 * direction.  */
static void rotate_2pic(plCurve *L,plc_vector *dir,bool dir_given) {
  plc_strand *cp;
  int this,that;
  double dist = DBL_MAX;
  double first,second,third;
  plc_vector this_edge,edge1,edge2,edge3;
  double cur_dist;
  double cos_theta,sin_theta;
  int cmp,vert;
  bool ok;

  if (!dir_given) { 
    this = 0;
    that = 2;
    cp = &(L->cp[L->nc-1]);
    while (that < cp->nv-1) {
      this_edge = plc_normalize_vect(
          plc_vect_diff(cp->vt[this+1],cp->vt[this]),&ok);
      edge1 = plc_normalize_vect(
          plc_vect_diff(cp->vt[that],cp->vt[that-1]),&ok);
      edge2 = plc_normalize_vect(
          plc_vect_diff(cp->vt[that+1],cp->vt[that]),&ok);
      edge3 = plc_normalize_vect(
          plc_vect_diff(cp->vt[that+2],cp->vt[that+1]),&ok);
      first = plc_M_dot(this_edge,edge1);
      second = plc_M_dot(this_edge,edge2);
      third = plc_M_dot(this_edge,edge3);
      /* Go looking for the most-negative dot product */
      while (second - first > DBL_EPSILON ||
          (that < cp->nv-1 && second - third > DBL_EPSILON)) {
        first = second;
        edge1 = edge2;
        second = third;
        edge2 = edge3;
        that++;
        edge3 = plc_normalize_vect(
            plc_vect_diff(cp->vt[that+2],cp->vt[that+1]),&ok);
        third = plc_M_dot(this_edge,edge3);
      }
      cur_dist = plc_sq_dist(cp->vt[this],cp->vt[that]);
      if (cur_dist <= dist) {
        dist = cur_dist;
        *dir = plc_normalize_vect(plc_vect_diff(cp->vt[that],cp->vt[this]),&ok);
      }
      this++;
    }
  }
  cos_theta = dir->c[1];
  sin_theta = dir->c[0];

  /* Drop the convex hull */
  plc_drop_component(L,L->nc-1);

  /* Rotate */
  for (cmp = 0; cmp < L->nc; cmp++) {
    for (vert = 0; vert < L->cp[cmp].nv; vert++) {
      L->cp[cmp].vt[vert] = 
        plc_build_vect(cos_theta*L->cp[cmp].vt[vert].c[0] -
                        sin_theta*L->cp[cmp].vt[vert].c[1],
                        sin_theta*L->cp[cmp].vt[vert].c[0] +
                        cos_theta*L->cp[cmp].vt[vert].c[1],0.0);
    }
  }
}

static void writeout(FILE *out, plCurve *L, 
                     struct arg_str *format,
                     struct arg_str *osize,
                     const char * const revision,
                     plc_vector F_dir, plc_vector rot_dir) {
  bool pgf = false;
  bool eepic = false;
  bool vect = true;
  double min_x,max_x,min_y,max_y;
  int cmp,vert;
  plc_color cur_col = {0.0,0.0,0.0,1.0};
  plc_color new_col;
  char colorspec[80];
  char fmttmp[80],*fmt;

  if (format->count > 0) {
    fmt = fmttmp;
    fmt = strncpy(fmttmp,format->sval[0],79);
    fmttmp[79] = '\0';
    if (*fmt == '=') { fmt++; }
#ifdef HAVE_STRCASECMP
    pgf = (strcasecmp(fmt,"pgf") == 0) || (strcasecmp(fmt,"tikz") == 0);
    eepic = (strcasecmp(fmt,"eepic") == 0) || (strcasecmp(fmt,"epic") == 0);
    vect = (strcasecmp(fmt,"vect") == 0);
#elif HAVE_STRCASESTR
    pgf = ((strlen(fmt) == 3) && (strcasestr(fmt,"pgf") == 0)) ||
          ((strlen(fmt) == 4) && (strcasestr(fmt,"tikz") == 0));
    eepic = ((strlen(fmt) == 5) && (strcasestr(fmt,"eepic") == 0)) ||
            ((strlen(fmt) == 4) && (strcasestr(fmt,"epic") == 0));
    vect = ((strlen(fmt) == 4) && (strcasestr(fmt,"vect") == 0));
#else
    pgf = (strcmp(fmt,"pgf") == 0) || (strcmp(fmt,"tikz") == 0);
    eepic = (strcmp(fmt,"eepic") == 0) || (strcmp(fmt,"epic") == 0);
    vect = (strcmp(fmt,"vect") == 0);
#endif
  }
  if (!pgf && !eepic && !vect) {
    fprintf(stderr,"Error: Unknown output format '%s'.\n",fmttmp);
    exit(EXIT_FAILURE);
  }
  min_x = DBL_MAX;
  min_y = DBL_MAX;
  max_x = DBL_MIN;
  max_y = DBL_MIN;
  if (pgf || eepic) {
    /* Find the bounding box */
    for (cmp = 0; cmp < L->nc; cmp++) {
      for (vert = 0; vert < L->cp[cmp].nv; vert++) {
        if (L->cp[cmp].vt[vert].c[0] > max_x) { 
          max_x = L->cp[cmp].vt[vert].c[0];
        }
        if (L->cp[cmp].vt[vert].c[0] < min_x) { 
          min_x = L->cp[cmp].vt[vert].c[0];
        }
        if (L->cp[cmp].vt[vert].c[1] > max_y) { 
          max_y = L->cp[cmp].vt[vert].c[1];
        }
        if (L->cp[cmp].vt[vert].c[1] < min_y) { 
          min_y = L->cp[cmp].vt[vert].c[1];
        }
      }
    }
    fprintf(out,"%% Knot code produced using knot_diagram version %s.\n",
      &revision[11]);
    fprintf(out,"%% Projection: %g,%g,%g  Rotation: %g,%g,%g\n",
      plc_M_clist(F_dir),plc_M_clist(rot_dir));
    fprintf(out, 
      "%% To change the size of the picture, adjust the ``%s'' below.\n",
      osize->sval[0]);
  }
  if (pgf) {
    fprintf(out,
      "%% Remember to \\usepackage{tikz} in the preamble\n");
    fprintf(out,"%%\n%% We have to define our own arrow head since pgf's "
      "arrows all \n%%   shorten their lines somewhat, which causes them to "
      "point\n%%   in the wrong direction for short edge lengths.\n");
    fprintf(out,"\\ifdefined\\knotDiagramArrowDimen\\else\\newdimen{\\knotDiagramArrowDimen}\n");
    fprintf(out,"\\pgfarrowsdeclare{toKnotDiag}{toKnotDiag}\n");
    fprintf(out,"{\n");
    fprintf(out,"  \\pgfarrowsleftextend{0pt}\n");
    fprintf(out,"  \\pgfarrowsrightextend{0pt}\n");
    fprintf(out,"}\n");
    fprintf(out,"{\n");
    fprintf(out,"  \\knotDiagramArrowDimen=0.28pt%%\n");
    fprintf(out,"  \\advance\\knotDiagramArrowDimen by.3\\pgflinewidth%%\n");
    fprintf(out,"  \\pgfsetlinewidth{0.8\\pgflinewidth}\n");
    fprintf(out,"  \\pgfsetdash{}{0pt}\n");
    fprintf(out,"  \\pgfsetroundcap\n");
    fprintf(out,"  \\pgfsetroundjoin\n");
    fprintf(out,"  \\pgfpathmoveto{\\pgfpoint{-3\\knotDiagramArrowDimen}{4\\knotDiagramArrowDimen}}\n");
    fprintf(out,"  \\pgfpathcurveto\n");
    fprintf(out,"  {\\pgfpoint{-2.75\\knotDiagramArrowDimen}{2.5\\knotDiagramArrowDimen}}\n");
    fprintf(out,"  {\\pgfpoint{0pt}{0.25\\knotDiagramArrowDimen}}\n");
    fprintf(out,"  {\\pgfpoint{0.75\\knotDiagramArrowDimen}{0pt}}\n");
    fprintf(out,"  \\pgfpathcurveto\n");
    fprintf(out,"  {\\pgfpoint{0pt}{-0.25\\knotDiagramArrowDimen}}\n");
    fprintf(out,"  {\\pgfpoint{-2.75\\knotDiagramArrowDimen}{-2.5\\knotDiagramArrowDimen}}\n");
    fprintf(out,"  {\\pgfpoint{-3\\knotDiagramArrowDimen}{-4\\knotDiagramArrowDimen}}\n");
    fprintf(out,"  \\pgfusepathqstroke\n");
    fprintf(out,"}\n");
    fprintf(out,"\\fi\n");
    fprintf(out,"\\ifdefined\\knotDiagramLength\\else"
                "\\newlength{\\knotDiagramLength}\\fi\n");
    fprintf(out,"\\setlength{\\knotDiagramLength}{%s / \\real{%3.3f}}\n",
      osize->sval[0],
      (max_x - min_x > max_y - min_y) ? max_x - min_x : max_y - min_y);
    fprintf(out,
      "\\begin{tikzpicture}[x=\\knotDiagramLength,y=\\knotDiagramLength]\n");
    fprintf(out,"\\tikzstyle{firsthalf}=[thick,-toKnotDiag]\n");
    fprintf(out,"\\tikzstyle{secondhalf}=[thick]\n");
    for (cmp = 0; cmp < L->nc; cmp++) {
      if (L->cp[cmp].cc == 1 &&
          (L->cp[cmp].clr[0].r != 0.0 ||
           L->cp[cmp].clr[0].g != 0.0 ||
           L->cp[cmp].clr[0].b != 0.0)) {
        sprintf(colorspec,",color={rgb,1:red,%3.3f;green,%3.3f;blue,%3.3f}",
          L->cp[cmp].clr[0].r,L->cp[cmp].clr[0].g,L->cp[cmp].clr[0].b);
      } else {
        colorspec[0] = '\0';
      }
      fprintf(out,"\\draw[firsthalf%s] (%3.3f,%3.3f)",colorspec,
        L->cp[cmp].vt[0].c[0],L->cp[cmp].vt[0].c[1]);
      for (vert = 1; vert < L->cp[cmp].nv / 2; vert++) { 
        fprintf(out," -- (%3.3f,%3.3f)",L->cp[cmp].vt[vert].c[0],
          L->cp[cmp].vt[vert].c[1]);
      }
      fprintf(out,";\n\\draw[secondhalf%s] (%3.3f,%3.3f)",colorspec,
        L->cp[cmp].vt[vert-1].c[0],L->cp[cmp].vt[vert-1].c[1]);
      /* We overwrite the previous line segment so as to get continuity */
      for (; vert < L->cp[cmp].nv; vert++) { 
        fprintf(out," -- (%3.3f,%3.3f)",L->cp[cmp].vt[vert].c[0],
          L->cp[cmp].vt[vert].c[1]);
      }
      fprintf(out,";\n");
    }
    fprintf(out,"\\end{tikzpicture}");
  } else if (eepic) {
    fprintf(out,
      "%% Remember to \\usepackage{epic,eepic,color,calc} in the preamble\n");
    /* Default to a picture with maximum dimension 2 inches. */
    fprintf(out,"{ %% To keep the effects of \\setlength contained \n");
    fprintf(out,"\\setlength{\\unitlength}{%s / \\real{%3.3f}}\n", 
      osize->sval[0],
      (max_x - min_x > max_y - min_y) ? max_x - min_x : max_y - min_y);
    fprintf(out,"\\begin{picture}(%3.3f,%3.3f)(%3.3f,%3.3f)\n", max_x - min_x,
      max_y - min_y, min_x, min_y);
    for (cmp = 0; cmp < L->nc; cmp++) {
      if (L->cp[cmp].cc != 1) {
        /* In case of no color give or of different colors for each vertex,
           set the color for this strand to black */
        new_col.r = 0.0; new_col.g = 0.0; new_col.b = 0.0;
      } else {
        new_col = L->cp[cmp].clr[0];
      }
      if (new_col.r != cur_col.r ||
          new_col.g != cur_col.g ||
          new_col.b != cur_col.b) {
        if (cur_col.r != 0.0 || cur_col.g != 0.0 || cur_col.b != 0.0) {
          fprintf(out,"\\end{color}%%\n");
        }
        if (new_col.r != 0.0 || new_col.g != 0.0 || new_col.b != 0.0) {
          fprintf(out,"\\begin{color}[rgb]{%3.3f,%3.3f,%3.3f}\n",
              new_col.r,new_col.g,new_col.b);
        }
        cur_col = new_col;
      }
      fprintf(out,"\\drawline ");
      for (vert = 0; vert < L->cp[cmp].nv; vert++) {
        fprintf(out,"(%3.3f,%3.3f) ",L->cp[cmp].vt[vert].c[0],
                                     L->cp[cmp].vt[vert].c[1]);
      }
      fprintf(out,"\n");
    }
    if (cur_col.r != 0.0 || cur_col.g != 0.0 || cur_col.b != 0.0) {
      fprintf(out,"\\end{color}\n");
    }
    fprintf(out,"\\end{picture}\n}\n");
  } else { 
    plc_write(out,L);
    fprintf(out,"# Projection: %g,%g,%g  Rotation: %g,%g,%g\n",
      plc_M_clist(F_dir),plc_M_clist(rot_dir));
  }
}
/*
Error: 3 = cmp >= L->nc = 3  Angle: 0.963819,-0.175632,0.200513
Error: 7 = cmp >= L->nc = 7  Angle: -0.352259,-0.609383,0.710328
*/

int main(int argc, char *argv[]) {
  FILE *vectfile;
  plCurve *L, *G;
  plCurve *F = NULL; 
  char err_str[80];
  int err_num;
  double max_edge,cur_edge;
  int cmp,vert;
  int cnt;
  int verts_left, max_verts_left;
  double total_curvature, min_total_curvature;

  int show_work = 0;
  int tries = 20;
  int delay = 16000;
  FILE *geomview = NULL;
  char revision[20] = "$Revision: 1.29 $";
  char *dollar;
  char plcv[80];

  plc_vector direction;
  plc_vector F_dir;
  plc_vector rot_dir;

#ifdef HAVE_ARGTABLE2_H
  struct arg_file *filename = arg_file1(NULL,NULL,"<file>",
      "The input VECT file");
  struct arg_file *outname = arg_file0("o","outfile","<file>",
      "The output VECT file, by default stdout");
  struct arg_str *format = arg_str0("O","format","<str>",
      "Output format (vect,pgf,eepic)");
  struct arg_str *osize = arg_str0("S","size","<tex length>",
      "Size of largest dimension for pgf/eepic [2in]");
  struct arg_lit *help = arg_lit0("h","help","Print this help and exit");
  struct arg_int *view = arg_int0("v","view","<int>",
      "search view level (0,5,10,15,20,25) [0]");
  struct arg_int *tryopt = arg_int0("t","tries","<int>",
      "Number of different views to choose amongst [20]");
  struct arg_int *delayopt = arg_int0(NULL,"delay","<int>",
      "Delay in ms to use with -v 15 and -v 25 [16]");
  struct arg_dbl *rot_x = arg_dbl0("p","rot_x","<dbl>",
      "x comp. of rot. vect. (forces -t 1)");
  struct arg_dbl *rot_y = arg_dbl0("q","rot_y","<dbl>",
      "y comp. of rot. vect.");
  struct arg_dbl *proj_x = arg_dbl0("x","proj_x","<dbl>",
      "x comp. of proj. vect. (forces -t 1)");
  struct arg_dbl *proj_y = arg_dbl0("y","proj_y","<dbl>",
      "y comp. of proj. vect.");
  struct arg_dbl *proj_z = arg_dbl0("z","proj_z","<dbl>",
      "z comp. of proj. vect.");
  struct arg_end *end = arg_end(20);

  void *argtable[] = {help,view,tryopt,delayopt,
      proj_x,proj_y,proj_z,rot_x,rot_y,outname,format,osize,filename,end};
  int nerrors;

#else
  if (argc < 2 || strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0) {
    fprintf(stderr,"usage: %s <file>\n",argv[0]);
    fprintf(stderr,"  (to get more options, rebuild using argtable2)\n");
    exit(EXIT_SUCCESS);
  }
  vectfile = fopen(argv[1],"r");
#endif
        
  direction = F_dir = rot_dir = plc_build_vect(0,0,0);

  dollar = strchr(&revision[1],'$');
  if (dollar != NULL) {
    dollar[-1] = '\0';
  }
  plc_version(plcv,80);
  fprintf(stderr,"knot_diagram v%s (plCurve %s)\n",&revision[11],plcv);
  fprintf(stderr,"  Produce a knot diagram from a VECT file.\n");

#ifdef HAVE_ARGTABLE2_H
  if (arg_nullcheck(argtable) != 0) {
    /* NULL entries detected, allocations must have failed. */
    fprintf(stderr,"%s: insufficient memory\n",argv[0]);
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    exit(EXIT_FAILURE);
  }

  /* Set the default output size */
  osize->sval[0] = "2in";

  /* Parse the command line as defined by argtable[] */
  nerrors = arg_parse(argc,argv,argtable);

  /* special case: '--help' takes preceence over error reporting */
  if (help->count > 0) {
    fprintf(stderr,"Usage: %s ",argv[0]);
    arg_print_syntax(stderr,argtable,"\n");
    arg_print_glossary(stderr,argtable,"  %-25s %s\n");
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return 0;
  }

  /* If the parser returned any errors then display them and exit */
  if (nerrors > 0) {
    /* Display the error details contained in the arg_end struct.*/
    fprintf(stderr,"\n");
    arg_print_errors(stderr,end,argv[0]);
    fprintf(stderr,"\nUsage: %s ",argv[0]);
    arg_print_syntax(stderr,argtable,"\n");
    arg_print_glossary(stderr,argtable,"  %-25s %s\n");
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return 1;
  }

  if (delayopt->count > 0) {
    delay = delayopt->ival[delayopt->count-1]*1000;
  }
  if (tryopt->count > 0) {
    tries = tryopt->ival[tryopt->count-1];
  }
  if (proj_x->count+proj_y->count+proj_z->count + 
      rot_x->count+rot_y->count > 0) {
    if (proj_x->count*proj_y->count*proj_z->count * 
        rot_x->count*rot_y->count==0) {
      fprintf(stderr,
        "Must specify complete projection and rotation vectors.\n");
      arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
      return 1;
    } 
    tries = 1;
    direction = plc_normalize_vect(plc_build_vect(
      proj_x->dval[proj_x->count-1],
      proj_y->dval[proj_y->count-1],
      proj_z->dval[proj_z->count-1]),NULL);
    rot_dir = plc_normalize_vect(plc_build_vect(
      rot_x->dval[rot_x->count-1],rot_y->dval[rot_y->count-1],0.0),NULL);
  }
  if (view->count > 0) {
    show_work = view->ival[view->count-1];
  }

  vectfile = fopen(filename->filename[0],"r");
#endif

  if (vectfile == NULL) {
    fprintf(stderr,"Unable to open file '%s' for reading.\n",
      filename->filename[0]);
    exit(EXIT_FAILURE);
  }
  L = plc_read(vectfile,&err_num,err_str,80);
  (void)fclose(vectfile);
  if (err_num != 0) {
#ifdef HAVE_ARGTABLE2_H
    fprintf(stderr,"Error reading file %s:\n%s\n",
        filename->filename[0],err_str);
#else
    fprintf(stderr,"Error reading file %s:\n%s\n",
        argv[1],err_str);
#endif
    exit(EXIT_FAILURE);
  }
  assert(L != NULL);

  srand((unsigned int)time((time_t *)NULL));

  /* Fire up geomview, if requested */
  if (show_work > 0) {
    geomview = popen("geomview -c -","w");
    if (geomview == NULL) {
      fprintf(stderr,"Unable to start geomview.");
      show_work = 0;
    } else {
      fprintf(geomview,"(normalization World none) (bbox-draw World no)\n");
      showCurve(L);
      if (show_work >= 5) { usleep(delay*100); }
    }
  }
  max_edge = 0.0;
  for (cmp = 0; cmp < L->nc; cmp++) {
    for (vert = 0; vert < L->cp[cmp].nv; vert++) {
      cur_edge = plc_distance(L->cp[cmp].vt[vert],L->cp[cmp].vt[vert+1]);
      max_edge = (max_edge >= cur_edge) ? max_edge : cur_edge;
    }
  }
  max_verts_left = 0;
  min_total_curvature = DBL_MAX;

  /* Main loop.  Try different directions */
  for (cnt = 0; cnt < tries; cnt++) {
    if (proj_x->count == 0) {
      direction = plc_random_vect();
    }
    G = flatten(L,direction,max_edge/1.2,2.0*max_edge,geomview,delay,show_work);
    assert(G != NULL);
    showCurve(G);
    if (show_work >= 5) { usleep(delay*100); }
    verts_left = 0;
    total_curvature = 0.0;
    for (cmp = 0; cmp < G->nc; cmp++) {
      verts_left += G->cp[cmp].nv;
      for (vert = 0; vert < G->cp[cmp].nv; vert++) {
        total_curvature += plc_MR_curvature(G,cmp,vert);
      }
    }
    if (show_work > 0) {
      fprintf(stderr,"Verts: %d  Total Curvature: %g",
          verts_left,total_curvature);
    }
    if ((verts_left > max_verts_left && 
        total_curvature/1.21 < min_total_curvature) ||
        (verts_left*1.0 >= max_verts_left*0.99 && 
         total_curvature < min_total_curvature)) {
      if (show_work > 0) { fprintf(stderr,"*\n"); }
      if (F != NULL) { plc_free(F); }
      F = G;
      F_dir = direction;
      max_verts_left = verts_left;
      min_total_curvature = total_curvature;
    } else {
      if (show_work > 0) { fprintf(stderr,"\n"); }
      plc_free(G);
    }
  }
  if (show_work > 0) {
    fprintf(stderr,
      "Final selection -- Verts: %d  Total Curvature: %g\n"
      "  Direction: %g,%g,%g\n",
        max_verts_left,min_total_curvature,plc_M_clist(F_dir));
  }
  if (show_work >= 10) {
    fprintf(geomview,"(delete Where)");
    if (show_work >= 20) {
      fprintf(geomview,"(delete Versus)");
    }
  }
  assert(F != NULL);
  add_convex_hull(F);
  showCurve(F);
  if (show_work >= 5 && rot_x->count == 0) { usleep(delay*250); }
  rotate_2pic(F,&rot_dir,rot_x->count > 0);
  showCurve(F);
#ifdef HAVE_ARGTABLE2_H
  if (outname->count > 0) {
    vectfile = fopen(outname->filename[0],"w");
    assert(vectfile != NULL);
    writeout(vectfile,F,format,osize,revision,F_dir,rot_dir);
    (void)fclose(vectfile);
  } else {
    writeout(stdout,F,format,osize,revision,F_dir,rot_dir);
  }
#else
  plc_write(stdout,F);
  fprintf(stdout,"# Projection: %g,%g,%g  Rotation: %g,%g,%g\n",
    plc_M_clist(F_dir),plc_M_clist(rot_dir));
#endif

  plc_free(F);
  plc_free(L);
#ifdef HAVE_ARGTABLE2_H
  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
#endif

  exit(EXIT_SUCCESS);
}
