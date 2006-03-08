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
                          bool *start_over) {
  int ctr;
  plCurve_strand temp_strand;
  plcl_vector *vt;
  plCurve_color *clr;

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
        L->cp[cmp].nv*sizeof(plcl_vector));
    if (L->cp[cmp].cc > L->cp[cmp].nv) { 
      L->cp[cmp].cc = L->cp[cmp].nv;
      memmove(L->cp[cmp].clr, &(L->cp[cmp].clr[vert+1]),
          L->cp[cmp].nv*sizeof(plCurve_color));
    }
    L->cp[cmp].open = true;
  } else {
    /* Middle vertex */
    if (L->cp[cmp].open) {
      plCurve_add_component(L,L->cp[cmp].nv-(vert+1),true,
          (L->cp[cmp].cc < L->cp[cmp].nv) ? L->cp[cmp].cc :
          L->cp[cmp].nv-(vert+1),&(L->cp[cmp].vt[vert+1]),
          (L->cp[cmp].cc < L->cp[cmp].nv) ? L->cp[cmp].clr :
          &(L->cp[cmp].clr[vert+1]));
      L->cp[cmp].nv = vert;
      if (L->cp[cmp].cc > L->cp[cmp].nv) { 
        L->cp[cmp].cc = L->cp[cmp].nv;
      }
      temp_strand = L->cp[L->nc-1];
      for (ctr = L->nc-1; ctr > cmp+1; ctr--) {
        L->cp[ctr] = L->cp[ctr-1];
        /*@-branchstate@*/
      }
      /*@=branchstate@*/
      L->cp[cmp+1] = temp_strand;
    } else {
      /* We're eliminating a middle vertex on a closed component.  We need to
         open the component, but do so *at that point*.  This will require
         moving the vertex after that one to the begining and the vertex before
         that one to the end */
      *start_over = true;
      L->cp[cmp].open = true;
      L->cp[cmp].nv--;
      vt = malloc(L->cp[cmp].nv*sizeof(plcl_vector));
      assert(vt != NULL);
      memcpy(vt,&(L->cp[cmp].vt[vert+1]),
          (L->cp[cmp].nv-vert)*sizeof(plcl_vector));
      memcpy(&(vt[L->cp[cmp].nv-vert]),L->cp[cmp].vt,vert*sizeof(plcl_vector));
      memcpy(L->cp[cmp].vt,vt,L->cp[cmp].nv*sizeof(plcl_vector));
      if (L->cp[cmp].cc > L->cp[cmp].nv) {
        clr = malloc(L->cp[cmp].nv*sizeof(plCurve_color));
        assert(clr != NULL);
        L->cp[cmp].cc = L->cp[cmp].nv;
        memcpy(clr,&(L->cp[cmp].clr[vert+1]),
            (L->cp[cmp].nv-vert)*sizeof(plCurve_color));
        memcpy(&(clr[L->cp[cmp].nv-vert]),L->cp[cmp].clr,
            vert*sizeof(plCurve_color));
        memcpy(L->cp[cmp].clr,clr,L->cp[cmp].nv*sizeof(plCurve_color));
        free(clr);
      }
      free(vt);
      /*@-branchstate@*/
    }
    /*@=branchstate@*/
  }
  if (L->cp[cmp].nv <= 1) {
    /* Ran out of edges here, get rid of the component */
    plCurve_drop_component(L,cmp);
    *start_over = true;
  }
  plCurve_fix_wrap(L);
}

#define pi 3.1415926
#define showCurve(F) \
  if (show_work >= 5) { \
    fprintf(geomview,"(geometry Curve "); \
    plCurve_write(geomview,F); \
    fprintf(geomview,") (look-recenter Curve) (look-encompass Curve)\n"); \
  }
/* Project the vertices onto the plane to which N is normal (returning them
 * as (x,y,0) vectors in a plCurve framework. */
static plCurve *flatten(const plCurve * const L, 
                        const plcl_vector N, 
                        const double gap,
                        const double neighbor_gap,
                              FILE *geomview,
                        const int show_work) {
  plCurve *F,*H;  /* Flattened knot, heights of flattened points */
  int cmp,cmp2,vert,vert2;
  plcl_vector n;
  plcl_vector *vt;
  plcl_vector e3;
  double cos_theta,sin_theta;
  plcl_vector axle,u_para,u_perp,v_cross_u;
  bool ok,far_enough;
  double gap_sq = gap*gap;
  double n_gap_sq = neighbor_gap*neighbor_gap;
  double dist;
  int last_to_check;
  bool start_over;
  int wherecnt;

  F = plCurve_copy(L);

  printf("Flattening using vector (%g,%g,%g) . . .\n",plcl_M_clist(N));
  /* First rotate so that the given normal vector lies along the z axis. */
  n = plcl_normalize_vect(N,NULL); /* Bail out if we give it a bad vector */
  e3 = plcl_build_vect(0.0,0.0,1.0);
  /* What should we rotate around? */
  axle = plcl_cross_prod(n,e3);
  ok = true;
  axle = plcl_normalize_vect(axle,&ok);
  if (ok) { 
    /* We actually do need to do some rotation */
    for (cmp = 0; cmp < F->nc; cmp++) {
      for (vert = 0; vert < F->cp[cmp].nv; vert++) {
        vt = &F->cp[cmp].vt[vert];
        /* Rotate space so that n is pointing along the z axis and just take 
           the first two components.  Formulas from
             http://en.wikipedia.org/wiki/Coordinate_rotation#Three_dimensions
         */
        cos_theta = plcl_dot_prod(n,e3);
        sin_theta = sqrt(1-cos_theta*cos_theta);
        u_para = plcl_scale_vect(plcl_dot_prod(*vt,axle),axle);
        u_perp = plcl_vect_diff(*vt,u_para);
        v_cross_u = plcl_cross_prod(*vt,axle);
        *vt = plcl_vect_sum(u_para,
            plcl_vlincomb(cos_theta,u_perp,sin_theta,v_cross_u));
      }
    }
  }
  showCurve(F);

  H = plCurve_copy(F);
  assert(H != NULL);
  /* Flatten F.  Retain heights in H. */
  for (cmp = 0; cmp < F->nc; cmp++) {
    for (vert = 0; vert < F->cp[cmp].nv; vert++) {
      F->cp[cmp].vt[vert].c[2] = 0.0;
    }
  }
  showCurve(F);

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
        fflush(geomview);
        if (show_work >= 15) {
          usleep(16000);
        }
      }
      far_enough = false;
      vt = &(F->cp[cmp].vt[vert]);
      last_to_check = F->cp[cmp].nv-1;
      if (!F->cp[cmp].open && plcl_sq_dist(*vt,F->cp[cmp].vt[0]) <= n_gap_sq) {
        /* We are close to the beginning of a closed strand, perhaps we
         * should avoid checking some of the end vertices. */
        while (last_to_check > vert &&
            plcl_sq_dist(*vt,F->cp[cmp].vt[last_to_check]) <= n_gap_sq) {
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
            fflush(geomview);
            if (show_work >= 25) {
              usleep(16000);
            }
          }
          /* Wander along, looking for any time the other point comes within
           * gap of cmp:vert (after it first gets at least gap away from it).
           * That marks a crossing and we remove one of the vertices in
           * question -- the lower one. */
          dist = plcl_sq_dist(F->cp[cmp].vt[vert],F->cp[cmp2].vt[vert2]);
          far_enough = far_enough || (dist >= n_gap_sq);
          if ((dist <= gap_sq) && far_enough) {
            /* Found a crossing */
            start_over = false;
            if (H->cp[cmp].vt[vert].c[2] >= H->cp[cmp2].vt[vert2].c[2]) {
              remove_vertex(F,cmp2,vert2,&start_over);
              remove_vertex(H,cmp2,vert2,&start_over);
              showCurve(F);
              vert2--;
              if (cmp2 == cmp) {
                /* This is now an open strand, check all vertices */
                last_to_check = F->cp[cmp].nv-1;
              }
              if (start_over) {
                if (cmp2 != cmp) {
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
              remove_vertex(F,cmp,vert,&start_over);
              remove_vertex(H,cmp,vert,&start_over);
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
            
  plCurve_free(H);
  return F;
}

/* Add a component which contains the 2d convex hull of the other components
 * (looking only at the first 2 dimensions, of course).  We'll be using the
 * gift-wrapping algorithm. */
static void add_convex_hull(plCurve *L) {
  plcl_vector *vt;
  int nv,verts;
  int cmp,vert;
  plcl_vector first,next;
  plcl_vector x_vect,y_vect,cross;
  plCurve_color *clr;

  for (nv = 0, cmp = 0; cmp < L->nc; cmp++) {
    nv += L->cp[cmp].nv;
  }
  vt = malloc(nv*sizeof(plcl_vector));
  assert(vt != NULL);
  clr = malloc(nv*sizeof(plcl_vector));
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

  next = plcl_build_vect(0.0,0.0,0.0);
  /*@-compdef@*/
  while (!plcl_vecteq(vt[0],next) && verts < nv) {
  /*@=compdef@*/
    /* There is one other vertex out there which, when it is the other end of a
     * line segment with "first" all of the other vectors are to the left of
     * the line.  Find it in as pedantic a way as possible, ye olde brute
     * search (and a particularly stupid one at that). */
    next = (plcl_vecteq(L->cp[0].vt[0],first)) ?
      L->cp[0].vt[1] : L->cp[0].vt[0];
    x_vect.c[0] = first.c[0];
    y_vect.c[0] = first.c[1];
    x_vect.c[1] = next.c[0];
    y_vect.c[1] = next.c[1];
    for (cmp = 0; cmp < L->nc; cmp++) {
      for (vert = 0; vert < L->cp[cmp].nv; vert++) {
        if (!plcl_vecteq(L->cp[cmp].vt[vert],first) &&
            !plcl_vecteq(L->cp[cmp].vt[vert],next)) {
          x_vect.c[2] = L->cp[cmp].vt[vert].c[0];
          y_vect.c[2] = L->cp[cmp].vt[vert].c[1];
          cross = plcl_cross_prod(x_vect,y_vect);
          if (- cross.c[0] - cross.c[1] - cross.c[2] > DBL_EPSILON) {
            /* cmp:vert is to the right of (or same as) next, take it instead */
            next = L->cp[cmp].vt[vert];
            x_vect.c[1] = x_vect.c[2];
            y_vect.c[1] = y_vect.c[2];
          }
        }
      }
    }
    if (plcl_vecteq(next,vt[verts-1])) {
      fprintf(stderr,"Looping...\n");
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
  plCurve_add_component(L,verts,false,verts,vt,clr);
  free(vt);
  free(clr);
}

/* Rotate (in 2-dimensions) the first n-1 components based on the convex hull
 * in the nth component, so that the shortest "diameter" lies in the y
 * direction.  */
static void rotate_2pic(plCurve *L) {
  plCurve_strand *cp;
  int this,that;
  plcl_vector dir;
  double dist = DBL_MAX;
  double first,second,third;
  plcl_vector this_edge,edge1,edge2,edge3;
  double cur_dist;
  double cos_theta,sin_theta;
  int cmp,vert;
  bool ok;

  this = 0;
  that = 2;
  cp = &(L->cp[L->nc-1]);
  while (that < cp->nv-1) {
    this_edge = plcl_normalize_vect(
        plcl_vect_diff(cp->vt[this+1],cp->vt[this]),&ok);
    edge1 = plcl_normalize_vect(
        plcl_vect_diff(cp->vt[that],cp->vt[that-1]),&ok);
    edge2 = plcl_normalize_vect(
        plcl_vect_diff(cp->vt[that+1],cp->vt[that]),&ok);
    edge3 = plcl_normalize_vect(
        plcl_vect_diff(cp->vt[that+2],cp->vt[that+1]),&ok);
    first = plcl_M_dot(this_edge,edge1);
    second = plcl_M_dot(this_edge,edge2);
    third = plcl_M_dot(this_edge,edge3);
    /* Go looking for the most-negative dot product */
    while (second - first > DBL_EPSILON ||
        (that < cp->nv-1 && second - third > DBL_EPSILON)) {
      first = second;
      edge1 = edge2;
      second = third;
      edge2 = edge3;
      that++;
      edge3 = plcl_normalize_vect(
          plcl_vect_diff(cp->vt[that+2],cp->vt[that+1]),&ok);
      third = plcl_M_dot(this_edge,edge3);
    }
    cur_dist = plcl_sq_dist(cp->vt[this],cp->vt[that]);
    if (cur_dist <= dist) {
      dist = cur_dist;
      dir = plcl_normalize_vect(plcl_vect_diff(cp->vt[that],cp->vt[this]),&ok);
    }
    this++;
  }
  cos_theta = dir.c[1];
  sin_theta = dir.c[0];

  /* Drop the convex hull */
  plCurve_drop_component(L,L->nc-1);

  /* Rotate */
  for (cmp = 0; cmp < L->nc; cmp++) {
    for (vert = 0; vert < L->cp[cmp].nv; vert++) {
      L->cp[cmp].vt[vert] = 
        plcl_build_vect(cos_theta*L->cp[cmp].vt[vert].c[0] -
                        sin_theta*L->cp[cmp].vt[vert].c[1],
                        sin_theta*L->cp[cmp].vt[vert].c[0] +
                        cos_theta*L->cp[cmp].vt[vert].c[1],0.0);
    }
  }
}

int main(int argc, char *argv[]) {
  FILE *vectfile;
  plCurve *L, *F, *G;
  char err_str[80];
  int err_num;
  double max_edge,cur_edge;
  int cmp,vert;
  int cnt;
  int verts_left, max_verts_left;

  int show_work = 10;
  FILE *geomview;

#ifdef HAVE_ARGTABLE2_H
  if (argc < 2 || strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0) {
    printf("usage: %s <file>\n",argv[0]);
    exit(EXIT_SUCCESS);
  }
  vectfile = fopen(argv[1],"r");
#endif
        
  assert(vectfile != NULL);
  L = plCurve_read(vectfile,&err_num,err_str,80);
  assert(L != NULL);
  (void)fclose(vectfile);
  if (err_num != 0) {
    fprintf(stderr,"Error reading file %s:\n%s\n","kl_2_2_1_I.vect",err_str);
    exit(EXIT_FAILURE);
  }

  srand((unsigned int)time((time_t *)NULL));

  /* Fire up geomview, if requested */
  if (show_work > 0) {
    geomview = popen("geomview -c -","w");
    if (geomview == NULL) {
      printf("Unable to start geomview.");
      show_work = 0;
    } else {
      fprintf(geomview,"(normalization World none) (bbox-draw World no)\n");
      showCurve(L);
    }
  }
  max_edge = 0.0;
  for (cmp = 0; cmp < L->nc; cmp++) {
    for (vert = 0; vert < L->cp[cmp].nv; vert++) {
      cur_edge = plcl_distance(L->cp[cmp].vt[vert],L->cp[cmp].vt[vert+1]);
      max_edge = (max_edge >= cur_edge) ? max_edge : cur_edge;
    }
  }
  max_verts_left = 0;
  for (cnt = 0; cnt < 20; cnt++) {
    G = flatten(L,plcl_random_vect(),max_edge/1.4,2.0*max_edge,
          geomview,show_work);
    showCurve(G);
    verts_left = 0;
    for (cmp = 0; cmp < G->nc; cmp++) {
      verts_left += G->cp[cmp].nv;
    }
    if (verts_left > max_verts_left) {
      if (F != NULL) { plCurve_free(F); }
      F = G;
      max_verts_left = verts_left;
    } else {
      plCurve_free(G);
    }
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
  rotate_2pic(F);
  showCurve(F);
#ifdef HAVE_ARGTABLE2_H
  plCurve_write(stdout,F);
#else
  vectfile = fopen("kl_2_2_1_flat.vect","w");
  assert(vectfile != NULL);
  plCurve_write(vectfile,F);
  (void)fclose(vectfile);
#endif
  exit(EXIT_SUCCESS);
}
