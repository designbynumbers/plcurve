/* 
 * Take the given link, replacing each corner with its minrad circle arc (not
 * using the global MinRad value).  Place enough points along the resulting
 * curve(s) so that using an n^2 algorithm to check distances, we can find
 * its ropelength to within the desired accuracy.
 *
 * $Id: roundout_rl.c,v 1.16 2007-03-16 02:39:35 ashted Exp $
 *
 */

/* Copyright 2004 The University of Georgia. */

/* This file is part of liboctrope.
   
liboctrope is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

liboctrope is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with liboctrope; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif

#include <string.h>

#include"octrope_internal.h"
#include"octrope.h"

#include"argtable2.h"

#define PI 3.141592653589793238462643383279502884197169399375

typedef struct point_n_curvature_type {
  plc_vector point;
  double curvature;
} point_n_curvature;

typedef struct pair_of_interest_type {
  int comp1;
  int point1;
  int comp2;
  int point2;
  double dist;
} pair_of_interest;

/*inline*/ void add_arc_points(point_n_curvature *points, int *num_points,
                    plc_strand P, int vert, double accuracy,
                    double *total_len, double *total_curvature,
                    double *MinRad) {
  /* Add enough points from this arc to keep us within accuracy */
  double normin,normout;
  plc_vector in, out, cross_prod, center, first_rad, second_rad;
  double dot_prod,cross_prod_norm,rad;
  double radius, scale;
  double alpha, arc_len, angle;
  int steps,cnt;
  double halflength, position;
  plc_vector half_shortside;

  in = plc_vect_diff(P.vt[vert],P.vt[vert-1]);
  out = plc_vect_diff(P.vt[vert+1],P.vt[vert]);
  normin = plc_norm(in);
  normout = plc_norm(out);
  
  dot_prod = plc_M_dot(in,out);
  cross_prod = plc_cross_prod(in,out);
  cross_prod_norm = plc_norm(cross_prod);
  alpha = atan2(cross_prod_norm,dot_prod); 
  if (fabs(cross_prod_norm) < 1E-12) {
    /* This must be a straight line */
    halflength = (normin < normout) ? normin/2 : normout/2;
    for (position = normin - halflength; position < normin; 
         position += accuracy) {
      plc_M_vweighted(points[*num_points].point,position/normin,
                        P.vt[vert-1],P.vt[vert]);
      points[*num_points].curvature = *total_curvature;
      (*num_points)++;
    }
    for (position = 0; position < halflength; position += accuracy) {
      plc_M_vweighted(points[*num_points].point,position/normout,
                        P.vt[vert],P.vt[vert+1]);
      points[*num_points].curvature = *total_curvature;
      (*num_points)++;
    }
    *total_len += 2*halflength;
    *total_curvature += alpha;
    return;
  }
  /* Angle subtended by the minrad arc */
  rad = (normin*normout + dot_prod)/(2*cross_prod_norm);

  radius = rad*((normin < normout) ? normin : normout);
  if (radius < *MinRad) { 
    *MinRad = radius;
  }
  first_rad = plc_cross_prod(in,cross_prod);
  scale = radius/plc_norm(first_rad); /* For the macro's sake */
  plc_M_scale_vect(scale,first_rad);
  second_rad = plc_cross_prod(out,cross_prod);
  scale = radius/plc_norm(second_rad); /* For the macro's sake */
  plc_M_scale_vect(scale,second_rad);
  if (normin < normout) {
    center = plc_scale_vect(-1,first_rad);
    plc_M_vmadd(center,1/2,in);
    plc_M_vweighted(half_shortside,0.5,P.vt[vert-1],P.vt[vert]);
    plc_M_add_vect(center,half_shortside);
  } else {
    center = plc_scale_vect(-1,second_rad);
    plc_M_vmadd(center,1/2,out);
    plc_M_vweighted(half_shortside,0.5,P.vt[vert],P.vt[vert+1]);
    plc_M_add_vect(center,half_shortside);
  }
  arc_len = radius*alpha; /* Length of minrad arc */
  *total_len += arc_len;
  steps = arc_len/accuracy;
  for (cnt = 0; cnt < steps; cnt++) {
    angle = cnt*alpha/steps;
    points[*num_points].point = center;
    plc_M_vmadd(points[*num_points].point,cos(angle),first_rad);
    plc_M_vmadd(points[*num_points].point,sin(angle)*radius/normin,in);
    points[*num_points].curvature = *total_curvature+angle;
    (*num_points)++;
  }
  *total_curvature += alpha;
}

/*inline*/ int pair_compare(const void *A,const void *B) {
  pair_of_interest *a,*b;
  
  a = (pair_of_interest *)(A);
  b = (pair_of_interest *)(B);
  if (a->comp1 > b->comp1) { return 3; }
  if (a->comp1 < b->comp1) { return -3; }
  if (a->comp2 > b->comp2) { return 3; }
  if (a->comp2 < b->comp2) { return -3; }
  if (a->point1 > b->point1) { return 2; }
  if (a->point1 < b->point1) { return -2; }
  if (a->point2 > b->point2) { return 2; }
  if (a->point2 < b->point2) { return -2; }
  if (a->dist > b->dist) { return 1; }
  if (a->dist < b->dist) { return -1; }
  return 0;
}

/*inline*/ int wrap(plCurve*L,int *np_array,const int comp, int point) {
  if (point < 0) {
    point = (L->cp[comp].open) ?  0 : point + np_array[comp];
  } else if (point >= np_array[comp]) {
    point = (L->cp[comp].open) ?  np_array[comp] : point - np_array[comp];
  }
  return point;
}
  
/*inline*/ void add_three_points(plCurve*L,int *np_array,
                      point_n_curvature **points,double *curvatures,
                      char *argv[],pair_of_interest **poi,
                      const int poi_ptr,int *poi_cnt,
                      const int skip1,const int skip2,
                      int *poi_max) {
  int comp1,comp2,point1,point2,curve1,curve2,point1_skip,point2_skip,temp_int;

  comp1 = (*poi)[poi_ptr].comp1;
  comp2 = (*poi)[poi_ptr].comp2;
  /* if skip1  > 0 and skip2 == 0, just do the (-skip1,0) point.
   * if skip1 == 0 and skip2  > 0, just do the (0,-skip2) point.
   * otherwise, do points (0,-skip2), (-skip1,0) and (-skip1,-skip2) in that
   * order */
  for (point1_skip = (skip2) ? 0 : skip1; 
       point1_skip <= skip1; 
       point1_skip += (skip1) ? skip1 : 1) {
    for (point2_skip = (point1_skip) ? 0 : skip2; 
         point2_skip <= skip2;
         point2_skip += (skip2) ? skip2 : 1) {
      point1 = wrap(L,np_array,comp1,(*poi)[poi_ptr].point1 - point1_skip);
      point2 = wrap(L,np_array,comp2,(*poi)[poi_ptr].point2 - point2_skip);
      if (comp1 == comp2) { 
        if (point1 > point2) {
          temp_int = point1;
          point1 = point2;
          point2 = temp_int;
        }
        curve1 = points[comp1][point1].curvature;
        /* remember, comp1 == comp2 at this point */
        curve2 = points[comp1][point2].curvature;
      }
      if ((comp1 != comp2) || ((curve2 - curve1 >= PI) && 
           (L->cp[comp1].open ||
            (curvatures[comp1] + curve1 - curve2 >= PI)))) {
        (*poi)[*poi_cnt].comp1 = comp1;
        (*poi)[*poi_cnt].point1 = point1;
        (*poi)[*poi_cnt].comp2 = comp2;
        (*poi)[*poi_cnt].point2 = point2;
        (*poi)[*poi_cnt].dist = -1;
        (*poi_cnt)++;
        if (*poi_cnt >= *poi_max) {
          (*poi_max) += 10000;
          fprintf(stderr,"Increasing poi_max to %d (1) . . .\n",*poi_max);
          if ((*poi = 
            realloc(*poi,*poi_max*sizeof(pair_of_interest))) == NULL) {
            fprintf(stderr,
              "%s:Error -- Pairs of interest table overflow (%d).\n",
              argv[0],*poi_max);
            exit(3);
          }
        }
      }
    }
  }
}
  
#define min(A,B) ((A < B) ? A : B)
#define max(A,B) ((A > B) ? A : B)
int main(int argc,char *argv[]) {

  plCurve*L;
  FILE *infile;
  FILE *outfile;

  struct arg_file *knotfiles = 
    arg_filen(NULL,NULL,"<file>",1,1024,"input .VECT file");
  struct arg_file *outputfile =
    arg_file0("o","output","<file>","output .VECT file (for debugging)");
  struct arg_dbl *acc = arg_dbl0("a","accuracy","<dbl>",
    "Tolerance for length calculation (default=1e-6");
  struct arg_lit *help = arg_lit0("h","help","print this help and exit");
  struct arg_end *end = arg_end(20);

  void *argtable[] = {help,acc,outputfile,knotfiles,end};
  int nerrors;
  int filecnt;
  int *np_array;
  double *curvatures;

  point_n_curvature **points;
  int num_points, comp1, comp2, vert, cnt, scale_factor;
  int start_vert, end_vert;
  double total_len, total_curvature;
  double inlen, outlen, position;
  double side_len, longest_side;
  int point1, point2, point1w, point2w; 
  double shortest_distance = DBL_MAX;
  double MinRad = DBL_MAX;
  double thickness,thickbnd,esquared;
  pair_of_interest temp_pair;
  pair_of_interest *poi;
  int poi_cnt = 0; 
  int poi_cnt2 = 0;
  int poi_skip, pass_cnt, poi_ptr, top, bot, mid, temp_int, poi_max;
  int comparison;
  double curve1,curve2;
  double tolerance;
  char revision[20] = "$Revision: 1.16 $";
  char *dollar;

  dollar = strchr(&revision[1],'$');
  dollar[0] = '\0';
  printf("Roundout_RL v%s\n",&revision[11]);
  printf("  Replace the corners of a polygonal knot with (local) MinRad\n");
  printf(
    "  circles and estimate the ropelength of the resulting C^{1,1} knot.\n"
  );

  if (arg_nullcheck(argtable) != 0) {
    /* NULL entries detected, allocations must have failed. */
    fprintf(stderr,"%s: insufficient memory\n",argv[0]);
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return 1;
  }

  /* Set defaults */
  acc->dval[0] = 1e-6;
  outputfile->filename[0]="-";

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

  for (filecnt = 0; filecnt < knotfiles->count; filecnt++) {
    infile = fopen(knotfiles->filename[filecnt],"r");
    if (infile == NULL) {
      fprintf(stderr,"%s: Couldn't find input file %s.\n",argv[0],
        knotfiles->filename[filecnt]);
      return 1;
    }

    octrope_error_num = 0;
    L = plc_read(infile,&octrope_error_num,octrope_error_str,80);
    fclose(infile);
    if (octrope_error_num != 0) { 
      fprintf(stderr,"%s:%s",argv[0],octrope_error_str); 
      arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
      return 1;
    }

    /* Find the longest side so as to know how many points we need */

    longest_side = 0;
    num_points = 0;
    if ((points = malloc(L->nc*sizeof(point_n_curvature *))) == NULL) {
      fprintf(stderr,"Error, unable to allocate space for %d lists of points.\n", L->nc);
      exit(1);
    }
    if ((np_array = malloc(L->nc*sizeof(int))) == NULL) {
      fprintf(stderr,"Error, unable to allocate space for %d point counts.\n", L->nc);
      exit(1);
    }
    if ((curvatures = malloc(L->nc*sizeof(double))) == NULL) {
      fprintf(stderr,"Error, unable to allocate space for %d curvatures.\n", L->nc);
      exit(1);
    }
    for (comp1 = 0; comp1 < L->nc; comp1++) {
      num_points += L->cp[comp1].nv;
      if (L->cp[comp1].open) {
        start_vert = 1;
        end_vert = L->cp[comp1].nv;
      } else {
        start_vert = 0;
        end_vert = L->cp[comp1].nv;
      }
      for (vert = start_vert; vert < end_vert; vert++) {
        side_len =
          plc_norm(plc_vect_diff(L->cp[comp1].vt[vert-1],
                                      L->cp[comp1].vt[vert]));
        if (side_len > longest_side) {
          longest_side = side_len;
        }
      }
    }
    scale_factor = ceil(longest_side/acc->dval[0]);
    if (num_points*scale_factor >= INT32_MAX) {
      fprintf(stderr,"Requested accuracy results in too many points.\n");
      fprintf(stderr,"Use a number greater than %f.\n",
        num_points*longest_side/INT32_MAX);
      exit(1);
    }
    for (comp1 = 0; comp1 < L->nc; comp1++) {
      printf("Allocating %d.\n",L->cp[comp1].nv*scale_factor);
      if ((points[comp1] = 
           malloc(L->cp[comp1].nv*scale_factor*sizeof(point_n_curvature))) ==
           NULL) {
        fprintf(stderr,"Error, unable to allocate space for %d vectors.\n",
          L->cp[comp1].nv*scale_factor);
        fprintf(stderr,"You may wish to reduce the accuracy.\n");
        exit(1);
      }
      np_array[comp1] = 0;
    }
    poi_max = 10000;
    printf("Allocating space for %d pairs.\n",poi_max);
    if ((poi = calloc(poi_max,sizeof(pair_of_interest))) == NULL) {
      fprintf(stderr,"%s:Unable to allocate space for %d pairs.\n",argv[0],
        poi_max);
      fprintf(stderr,"You may wish to reduce the accuracy.\n");
      exit(1);
    }

    printf("Creating point list.\n");
    num_points = 0;
    total_len = 0;
    for (comp1 = 0; comp1 < L->nc; comp1++) {
      total_curvature = 0;
      if (L->cp[comp1].open) {
        start_vert = 1;
        end_vert = L->cp[comp1].nv-1;
        /* Do the first half-edge */
        inlen = plc_norm(plc_vect_diff(L->cp[comp1].vt[1],
                                            L->cp[comp1].vt[0]));
        for (position = 0; position < inlen/2; position += acc->dval[0]) {
          plc_M_vweighted(points[comp1][np_array[comp1]].point,
                            position/inlen,L->cp[comp1].vt[0],
                                           L->cp[comp1].vt[1]);
          points[comp1][np_array[comp1]].curvature = total_curvature;
          np_array[comp1]++;
        }
        total_len += inlen/2;
      } else {
        start_vert = 0;
        end_vert = L->cp[comp1].nv;
      }
      for (vert = start_vert; vert < end_vert; vert++) {
        inlen = plc_norm(plc_vect_diff(L->cp[comp1].vt[vert],
                                            L->cp[comp1].vt[vert-1]));
        outlen = plc_norm(plc_vect_diff(L->cp[comp1].vt[vert+1],
                                             L->cp[comp1].vt[vert]));
        /* Do the half straight-segment preceeding the curve */
        if ((inlen - outlen)/2 > acc->dval[0]) {
          for (position = inlen/2; 
               position < inlen - (outlen/2); 
               position += acc->dval[0]) {
            plc_M_vweighted(points[comp1][np_array[comp1]].point,
            position/inlen,L->cp[comp1].vt[vert-1],L->cp[comp1].vt[vert]);
            points[comp1][np_array[comp1]].curvature = total_curvature;
            np_array[comp1]++;
          }
        }
        /* Do the curve */
        add_arc_points(points[comp1],&np_array[comp1],L->cp[comp1],vert,
                       acc->dval[0],&total_len, &total_curvature,&MinRad);
        /* Do the half straight-segment following the curve */
        if ((outlen - inlen)/2 > acc->dval[0]) {
          for (position = inlen/2; 
               position < outlen/2;
               position += acc->dval[0]) {
            plc_M_vweighted(points[comp1][np_array[comp1]].point,
                              position/outlen,L->cp[comp1].vt[vert],
                                              L->cp[comp1].vt[vert+1]);
            points[comp1][np_array[comp1]].curvature = total_curvature;
            np_array[comp1]++;
          }
        }
        total_len += fabs((outlen - inlen)/2);
      }
      if (L->cp[comp1].open) {
        /* Do the last half-edge */
        outlen = plc_norm(
          plc_vect_diff(L->cp[comp1].vt[L->cp[comp1].nv],
                         L->cp[comp1].vt[L->cp[comp1].nv-1]));
        for (position = outlen/2; position < outlen; position += acc->dval[0]) {
          plc_M_vweighted(points[comp1][np_array[comp1]].point,
                            position/outlen,L->cp[comp1].vt[L->cp[comp1].nv-2],
                                            L->cp[comp1].vt[L->cp[comp1].nv-1]);
          points[comp1][np_array[comp1]].curvature = total_curvature;
          np_array[comp1]++;
        }
        /* Including the final point */
        points[comp1][np_array[comp1]].point = 
          L->cp[comp1].vt[L->cp[comp1].nv-1];
        points[comp1][np_array[comp1]].curvature = total_curvature;
        np_array[comp1]++;
        total_len += outlen/2;
      }  
      curvatures[comp1] = total_curvature;
      num_points += np_array[comp1];
    }

    /* Produce VECT file */
    if (outputfile->filename[0][0] != '-') {
      printf("Opening output file: %s\n",outputfile->filename[0]);
      outfile = fopen(outputfile->filename[0],"w");
      if (outfile == NULL) {
        fprintf(stderr,"%s: Couldn't create output file %s.\n",argv[0],
          outputfile->filename[0]);
      } else {
        printf("Writing output file: %s\n",outputfile->filename[0]);
        fprintf(outfile,"VECT\n");
        fprintf(outfile,"%d %d %d\n",L->nc,num_points,L->nc);
        for (comp1 = 0; comp1 < L->nc; comp1++) {
          if (L->cp[comp1].open) {
            fprintf(outfile,"%d ",np_array[comp1]);
          } else {
            fprintf(outfile,"%d ",-np_array[comp1]);
          }
        }
        fprintf(outfile,"\n");
        for (comp1 = 0; comp1 < L->nc; comp1++) {
          fprintf(outfile,"1 ");
        }
        fprintf(outfile,"\n");
        for (comp1 = 0; comp1 < L->nc; comp1++) {
          for (cnt = 0; cnt < np_array[comp1]; cnt++) {
            fprintf(outfile,"%f %f %f\n",points[comp1][cnt].point.c[0],
                                points[comp1][cnt].point.c[1],
                                points[comp1][cnt].point.c[2]);
          }
        }
        for (comp1 = 0; comp1 < L->nc; comp1++) {
          fprintf(outfile,"0 0 1 1\n");
        } 
        fclose(outfile);
      }
      printf("File written.\n");
    }

    /* Measure distances */
    printf("Doing distance checks (n is %d)\n",num_points);
    /* Use an adaptive algorithm, since the smallest actual distance must be
     * near the smallest measured distance. */
    /* Pick some points */
    poi_skip = floor(pow(2.0,ceil(log(num_points/100.0)/log(2.0)))+0.5);
    poi_skip = (poi_skip) ? poi_skip : 1;
    printf("Poi_skip: %d\n",poi_skip);
    shortest_distance = DBL_MAX;
    for (comp1 = 0; comp1 < L->nc; comp1++) {
      for (comp2 = comp1; comp2 < L->nc; comp2++) {
        for (point1 = 0; point1 < np_array[comp1] - poi_skip; 
             point1 += poi_skip) {
          for (point2 = (comp1 == comp2) ? point1 + poi_skip : 0; 
               point2 < np_array[comp2]; point2 += poi_skip) {
            curve1 = points[comp1][point1].curvature;
            curve2 = points[comp2][point2].curvature;
            if ((comp1 != comp2) ||
                ((curve2 - curve1 >= PI) && 
                 (L->cp[comp1].open ||
                  (curvatures[comp1] + curve1 - curve2 >= PI)))) {
              poi[poi_cnt].comp1 = comp1;
              poi[poi_cnt].point1 = point1;
              poi[poi_cnt].comp2 = comp2;
              poi[poi_cnt].point2 = point2;
              poi[poi_cnt].dist = plc_norm(
                plc_vect_diff(points[comp1][point1].point,
                               points[comp2][point2].point));
              if (poi[poi_cnt].dist < shortest_distance) {
                shortest_distance = poi[poi_cnt].dist;
              } 
              poi_cnt++;
              if (poi_cnt >= poi_max) {
                poi_max += 10000;
                fprintf(stderr,"Increasing poi_max to %d (2) . . .\n",poi_max);
                if ((poi = 
                  realloc(poi,poi_max*sizeof(pair_of_interest))) == NULL) {
                  fprintf(stderr,
                    "%s:Error -- Pairs of interest table overflow (%d).\n",
                    argv[0],poi_max);
                  exit(1);
                }
              }
            }
          }
        }
      }
    }
    pass_cnt = 1;
    printf("  Pass %d: %d pairs\n",pass_cnt++,poi_cnt);
    while (poi_skip > 1) {
      /* Build new point list */
      poi_cnt2 = poi_cnt;
/* x is a point in the s,t plane which has distance within the current epsilon.
 * We check points 1-8.  If they have been measured and also are within 
 * epsilon, we let them provide points b-h respectively.  If they have been
 * measured but are not within epsilon, we put in the respective points.  If
 * they have not been measured, we put them in as well.  x is responsible for
 * the points marked +.
 * 1 b 2 c 3
 * d + + e e
 * 4 + x e 5
 * f g g h h
 * 6 g 7 h 8
 */
      /* Tolerance should be the smaller of (1 + K/4) e^2 and just e. */
      tolerance = (1 + 1/(4*MinRad))*poi_skip*acc->dval[0]*
                                     poi_skip*acc->dval[0];
      tolerance = min(tolerance,poi_skip*acc->dval[0]);
      for (poi_ptr = 0; poi_ptr < poi_cnt; poi_ptr++) {
        if (poi[poi_ptr].dist <= shortest_distance + tolerance) {
          comp1 = poi[poi_ptr].comp1;
          comp2 = poi[poi_ptr].comp2;
          for (point1 = poi[poi_ptr].point1 - poi_skip;
               point1 <= poi[poi_ptr].point1 + poi_skip;
               point1 += poi_skip) {
            for (point2 = poi[poi_ptr].point2 - poi_skip;
                 point2 <= poi[poi_ptr].point2 + poi_skip;
                 point2 += poi_skip) {
              if (point1 == poi[poi_ptr].point1 &&
                  point2 == poi[poi_ptr].point2) {
                /* The original point, add the three points */
                add_three_points(L,np_array,points,curvatures,argv,
                                 &poi,poi_ptr,&poi_cnt2,poi_skip >> 1,
                                 poi_skip >> 1,&poi_max);
              } else {
                /* wrap if necessary */
                point1w = wrap(L,np_array,comp1,point1);
                point2w = wrap(L,np_array,comp2,point2);
                if (comp1 == comp2) {
                  if (point1w > point2w) {
                    temp_int = point1w;
                    point1w = point2w;
                    point2w = temp_int;
                  }
                  curve1 = points[comp1][point1w].curvature;
                  curve2 = points[comp2][point2w].curvature;
                }
                if ((comp1 != comp2) ||
                    ((curve2 - curve1 >= PI) && 
                     (L->cp[comp1].open ||
                      (curvatures[comp1] + curve1 - curve2 >= PI)))) {
                  /* This is a point which might already have been checked,
                   * seek it. */
                  temp_pair = poi[poi_ptr];
                  temp_pair.point1 = point1w;
                  temp_pair.point2 = point2w;
                  temp_pair.dist = -1;
                  /* binary search */
                  comparison = pair_compare(&temp_pair,&poi[poi_ptr]);
                  if (comparison > 1) {
                    bot = poi_ptr;
                    top = poi_cnt - 1;
                  } else {
                    bot = 0;
                    top = poi_ptr;
                  }
                  mid = (top + bot)/2;
                  comparison = pair_compare(&temp_pair,&poi[mid]);
                  while ((comparison > 1 || comparison < -1) && top - bot > 1) {
                    if (comparison > 1) {
                      bot = mid;
                    } else {
                      top = mid;
                    }
                    mid = (top+bot)/2;
                    comparison = pair_compare(&temp_pair,&poi[mid]);
                  }
                  if (comparison > 1) {
                    mid = top;
                    comparison = pair_compare(&temp_pair,&poi[mid]);
                  }
                  if (comparison <= -2 || comparison >= 2) {
                    /* The point hasn't been checked, add it and its 
                     * (up to) 3 points */
                    poi[poi_cnt2] = temp_pair;
                    poi_cnt2++;
                    if (poi_cnt2 >= poi_max) {
                      poi_max += 10000;
                      fprintf(stderr,"Increasing poi_max to %d (3) . . .\n",
                        poi_max);
                      if ((poi = 
                        realloc(poi,poi_max*sizeof(pair_of_interest))) 
                        == NULL) {
                        fprintf(stderr,
                          "%s:Error -- Pairs of interest table overflow (%d).\n",
                          argv[0],poi_max);
                        exit(4);
                      }
                    }
                    /* Yes, these should be point1/2 and not point1w/2w, we're
                       checking to see where we are in the loop */
                    if (point1 >= poi[poi_ptr].point1 ||
                        point2 >= poi[poi_ptr].point2) {
                      add_three_points(L,np_array,points,curvatures,argv,
                                       &poi,poi_cnt2-1,&poi_cnt2,
                        (point1 < poi[poi_ptr].point1) ? 0 : poi_skip >> 1,
                        (point2 < poi[poi_ptr].point2) ? 0 : poi_skip >> 1,
                        &poi_max);
                    }
                  } else {
                    /* We found the pair */
                    if (poi[mid].dist > 
                        shortest_distance + tolerance) {
                      /* This one won't add its points, step in */
                      if (point1 >= poi[poi_ptr].point1 ||
                          point2 >= poi[poi_ptr].point2) {
                        add_three_points(L,np_array,points,curvatures,argv,
                                         &poi,mid,&poi_cnt2,
                          (point1 < poi[poi_ptr].point1) ? 0 : poi_skip >> 1,
                          (point2 < poi[poi_ptr].point2) ? 0 : poi_skip >> 1,
                          &poi_max);
                      } /* if */
                    } /* if */
                  } /* else */
                } /* if */
              } /* else */
            } /* for */
          } /* for */
        } /* If distance is reasonable */
      } /* for */
      /* Sort list */
      printf("Used %d of %d (%5.5f%%) in pair table.\n",poi_cnt2,poi_max,
        100.0*poi_cnt2/poi_max);
      qsort(poi,poi_cnt2,sizeof(pair_of_interest),pair_compare);
      /* Remove duplicates */
      poi_cnt = 0;
      for (poi_ptr = 0; poi_ptr < poi_cnt2; poi_ptr++) {
        if (poi_cnt == 0 ||
            poi[poi_ptr].comp1  != poi[poi_cnt-1].comp1 ||
            poi[poi_ptr].point1 != poi[poi_cnt-1].point1 ||
            poi[poi_ptr].comp2  != poi[poi_cnt-1].comp2 ||
            poi[poi_ptr].point2 != poi[poi_cnt-1].point2) {
          poi[poi_cnt] = poi[poi_ptr];
          poi_cnt++;
          if (poi_cnt >= poi_max) {
            poi_max += 10000;
            fprintf(stderr,"Increasing poi_max to %d (4) . . .\n",poi_max);
            if ((poi = realloc(poi,poi_max*sizeof(pair_of_interest))) == NULL) {
              fprintf(stderr,
                "%s:Error -- Pairs of interest table overflow (%d).\n",
                argv[0],poi_max);
              exit(1);
            }
          }
        }
      }
      poi_skip >>= 1;
      printf("  Pass %d: %d pairs (skip == %d)\n",pass_cnt++,poi_cnt,poi_skip);
      /* Find new lengths */
      for (poi_ptr = 0; poi_ptr < poi_cnt; poi_ptr++) {
        if (poi[poi_ptr].dist < 0) {
          poi[poi_ptr].dist = plc_norm(plc_vect_diff(
            points[poi[poi_ptr].comp1][poi[poi_ptr].point1].point,
            points[poi[poi_ptr].comp2][poi[poi_ptr].point2].point));
          if (poi[poi_ptr].dist < shortest_distance) {
            shortest_distance = poi[poi_ptr].dist;
          } 
        }
      }
    }
    thickness = (2*MinRad < shortest_distance) ? 2*MinRad : shortest_distance;
    thickbnd  = (2*MinRad < shortest_distance - acc->dval[0]) ? 
      2*MinRad : shortest_distance - acc->dval[0];
    esquared = shortest_distance - (1+1/(4*MinRad))*acc->dval[0]*acc->dval[0];
    esquared = (2*MinRad < esquared) ? 2*MinRad : esquared;
    printf("Shortest pairwise distance found: %15.15f\n",shortest_distance);
    printf("MinRad for the curve: %15.15f\n",MinRad);
    printf("Thickness used: measured: %15.15f\n",thickness);
    printf("               e-bounded: %15.15f\n",thickbnd);
    printf("             e^2-bounded: %15.15f\n",esquared);
    printf("Curve length: %15.15f\n",total_len);
    printf("Diameter Ropelength: measured: %15.15f\n",total_len/thickness);
    printf("                e-Upper Bound: %15.15f\n",total_len/thickbnd);
    printf("              e^2-Upper Bound: %15.15f\n",total_len/esquared);
    printf("Radius Ropelength: measured: %15.15f\n",2*total_len/thickness);
    printf("              e-Upper Bound: %15.15f\n",2*total_len/thickbnd);
    printf("            e^2-Upper Bound: %15.15f\n",2*total_len/esquared);

    free(poi);
    for (comp1 = 0; comp1 < L->nc; comp1++) {
      free(points[comp1]);
    } 
    free(points);
    free(np_array);
    plc_free(L);
  }
  return(0);
}
