/*
 * Check all the various ways of calling octrope.
 *
 * $Id: call_check.c,v 1.11 2006-04-18 19:14:36 ashted Exp $
 *
 */

/* Copyright 2004 The University of Georgia */

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

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <math.h>

#include "octrope.h"

/* Free and return */

#define fandr\
  if (L != NULL) { plc_free(L); }\
  if (L2 != NULL) { plc_free(L2); }\
  if (L3 != NULL) { plc_free(L3); }\
  return -1;

#define check_err\
  if (octrope_error_num != 0) {\
    fprintf(stderr,"%s:Octrope error %d:%s\n",\
      argv[0],octrope_error_num,octrope_error_str);\
    fandr;\
  }

int main(int argc,char *argv[]) {

  int num_struts;
  plCurve *L = NULL, *L2 = NULL, *L3 = NULL;
  octrope_strut strutlist[8];
  double sl_size = 8;
  double min_strut;
  void *mem;
  int memsize;
  double min_rad;
  octrope_mrloc min_rad_locs[25];
  int mr_size = 25;
  int num_min_rad_locs;
  double ropelength;
  double thickness;
  double curve_length;
  int nv,cc;
  bool open;
  FILE *test_file;
  char *filename = "call_check_testfile";
  int call_cnt,cnt;
  plc_vector strut[2];
  int i;
  
  /* We will attempt to make all the possible calls to the octrope library.
     These are:
       octrope_struts
       octrope_minrad
       octrope_curvelength
       octrope_poca
       octrope_minradval
       octrope_thickness
       octrope_ropelength
       octrope
       octrope_est_mem
       octrope_set_levels
       octrope_set_debug
       octrope_debug_level
       octrope_strut_ends
       plc_new
       plc_free
       plc_read
       plc_write
       plc_fix_wrap
       plc_num_edges
       plc_copy
    */


  /* Test setting debug levels */
  octrope_set_debug(97);
  check_err;
  if (octrope_debug_level() != 97) {
    fprintf(stderr,"%s:Error setting or checking octrope_debug_level.\n",
       argv[0]);
    fprintf(stderr,"%s:Attempted to set to 97, currently %d.\n",
      argv[0],octrope_debug_level());
  }
  check_err;
  octrope_set_debug(0);
  check_err;
  if (octrope_debug_level() != 0) {
    fprintf(stderr,"%s:Error setting or checking octrope_debug_level.\n",
      argv[0]);
    fprintf(stderr,"%s:Attempted to set to 97, currently %d.\n",
      argv[0],octrope_debug_level());
  }
  check_err;
  octrope_set_levels(1);
  check_err;
  octrope_set_levels(0);
  check_err;

  /* Now set up the link (a four-pointed star, folded in the middle)*/
  nv = 8;
  open = false;
  cc = 1;
  L = plc_new(1,&nv,&open,&cc);
  check_err;
  L->cp[0].vt[0].c[0] = 0;
  L->cp[0].vt[0].c[1] = 1;  /* (0,1,0) */
  L->cp[0].vt[0].c[2] = 0;
  
  L->cp[0].vt[1].c[0] = 4;
  L->cp[0].vt[1].c[1] = 4;  /* (4,4,0) */
  L->cp[0].vt[1].c[2] = 0;
  
  L->cp[0].vt[2].c[0] = 1;
  L->cp[0].vt[2].c[1] = 0;  /* (1,0,0) */
  L->cp[0].vt[2].c[2] = 0;
  
  L->cp[0].vt[3].c[0] = 4;
  L->cp[0].vt[3].c[1] = 0;  /* (4,0,4) */
  L->cp[0].vt[3].c[2] = 4;
  
  L->cp[0].vt[4].c[0] = 0;
  L->cp[0].vt[4].c[1] = 0;  /* (0,0,1) */
  L->cp[0].vt[4].c[2] = 1;
  
  L->cp[0].vt[5].c[0] = -4;
  L->cp[0].vt[5].c[1] = 0;  /* (-4,0,4) */
  L->cp[0].vt[5].c[2] = 4;
  
  L->cp[0].vt[6].c[0] = -1;
  L->cp[0].vt[6].c[1] = 0;  /* (-1,0,0) */
  L->cp[0].vt[6].c[2] = 0;
  
  L->cp[0].vt[7].c[0] = -4;
  L->cp[0].vt[7].c[1] = 4;  /* (-4,4,0) */
  L->cp[0].vt[7].c[2] = 0;

  L->cp[0].clr[0].r = 0;
  L->cp[0].clr[0].g = 0;
  L->cp[0].clr[0].b = 1;
  L->cp[0].clr[0].alpha = 1;

  plc_fix_wrap(L);
  check_err;

  if (plc_num_edges(L) != 8) {
    fprintf(stderr,"%s:Incorrect number of edges in link:%d\n",
      argv[0],plc_num_edges(L));
    fandr;
  }
  check_err;
  L2 = plc_copy(L);
  if (plc_num_edges(L2) != 8) {
    fprintf(stderr,"%s:Incorrect number of edges in copied link:%d\n",
      argv[0],plc_num_edges(L2));
    fandr;
  }
  if (L2->cp[0].vt[7].c[1] != 4) {
    fprintf(stderr,"%s:Link copied incorrectly.\n", argv[0]);
    fandr;
  }
  test_file = fopen(filename,"w");
  if (test_file == NULL) {
    fprintf(stderr,"%s: Unable to open %s for writing.\n",argv[0],filename);
    fandr;
  }
  plc_write(test_file,L);
  check_err;
  fclose(test_file);
  test_file = fopen(filename,"r");
  if (test_file == NULL) {
    fprintf(stderr,"%s: Unable to open %s for reading.\n",argv[0],filename);
    fandr;
  }
  octrope_error_num = 0;
  L3 = plc_read(test_file,&octrope_error_num,octrope_error_str,80);
  check_err;
  fclose(test_file);
  if (plc_num_edges(L3) != 8) {
    fprintf(stderr,"%s:Incorrect number of edges in written link:%d\n",
      argv[0],plc_num_edges(L2));
    fandr;
  }
  if (L3->cp[0].vt[7].c[1] != 4) {
    fprintf(stderr,"%s:Link written or read incorrectly.\n", argv[0]);
    fandr;
  }
  unlink(filename);
  
  memsize = octrope_est_mem(plc_num_edges(L));
  check_err;
  if ((mem = (void *)malloc(memsize)) == NULL) {
    fprintf(stderr,"%s:Unable to allocate %d bytes.\n",argv[0],memsize);
    fandr;
  }

  curve_length = octrope_curvelength(L);
  check_err;
  if (fabs(curve_length - 40) > 1E-12) {
    fprintf(stderr,"%s:Expected curve length to be 40, got %17.17f.\n",
      argv[0],curve_length);
    fandr;
  }
  min_strut = octrope_poca(L,mem,memsize);
  check_err;
  if (fabs(min_strut*min_strut - 2.0) > 1E-12) {
    fprintf(stderr,"%s:Expected POCA to be sqrt(2), got %17.17f.\n",
      argv[0],min_strut);
    fandr;
  }
  min_rad = octrope_minradval(L);
  check_err;
  if (fabs(min_rad - 5.0/14.0) > 1E-12) {
    fprintf(stderr,"%s:Expected MinRad to be 5/14, got %17.17f.\n",
      argv[0],min_strut);
    fandr;
  }
  /* We check the lambda = 1 case below */
  thickness = octrope_thickness(L,mem,memsize,0);
  check_err;
  if (fabs(thickness - min_strut/2.0) > 1E-12) {
    fprintf(stderr,"%s:Expected thickness to be %17.17f, got %17.17f.\n",
      argv[0],min_strut,thickness);
    fandr;
  }
  ropelength = octrope_ropelength(L,mem,memsize,0);
  check_err;
  if (fabs(ropelength - curve_length/thickness) > 1E-12) {
    fprintf(stderr,"%s:Expected ropelength to be %17.17f, got %17.17f.\n",
      argv[0],curve_length/thickness,ropelength);
    fandr;
  }

  min_strut = 25;
  strutlist[0].component[0] = -1;
  num_struts = octrope_struts(L,0,0,strutlist,sl_size,&min_strut,NULL,0);
  check_err;
  if (num_struts != 1) {
    fprintf(stderr,"%s:Expected a single strut, got %d.\n",
      argv[0],num_struts);
    fandr;
  }
  if (fabs(min_strut*min_strut - 2.0) > 1E-12) {
    fprintf(stderr,"%s:Expected strut length to be sqrt(2), got %17.17f.\n",
      argv[0],min_strut);
    fandr;
  }
  octrope_strut_ends(L,strutlist,strut);
  check_err;
  if ((strut[0].c[0] != 0 || strut[0].c[1] != 0 || strut[0].c[2] != 1 ||
       strut[1].c[0] != 0 || strut[1].c[1] != 1 || strut[1].c[2] != 0) &&
      (strut[0].c[0] != 0 || strut[0].c[1] != 1 || strut[0].c[2] != 0 ||
       strut[1].c[0] != 0 || strut[1].c[1] != 0 || strut[1].c[2] != 1)) {
    fprintf(stderr,
      "%s:Expected a strut (0,0,1)-(0,1,0), got (%f,%f,%f)-(%f,%f,%f).\n",
      argv[0],strut[0].c[0],strut[0].c[1],strut[0].c[2],
              strut[1].c[0],strut[1].c[1],strut[1].c[2]);
    fandr;
  }
  min_strut = 25;
  strutlist[0].component[0] = -1;
  num_struts = octrope_struts(L,2,0,strutlist,sl_size,&min_strut,NULL,0);
  check_err;
  if (num_struts != 2) {
    fprintf(stderr,"%s:Expected two struts, got %d.\n",
      argv[0],num_struts);
    fandr;
  }
  if (fabs(min_strut*min_strut - 2.0) > 1E-12) {
    fprintf(stderr,"%s:Expected strut length to be sqrt(2), got %17.17f.\n",
      argv[0],min_strut);
    fandr;
  }
  octrope_strut_ends(L,strutlist,strut);
  check_err;
  if ((strut[0].c[0] != 1 || strut[0].c[1] != 0 || strut[0].c[2] != 0 ||
       strut[1].c[0] != -1|| strut[1].c[1] != 0 || strut[1].c[2] != 0) &&
      (strut[0].c[0] != -1|| strut[0].c[1] != 0 || strut[0].c[2] != 0 ||
       strut[1].c[0] != 1 || strut[1].c[1] != 0 || strut[1].c[2] != 0)) {
    fprintf(stderr,
      "%s:Expected a strut (1,0,0)-(-1,0,0), got (%f,%f,%f)-(%f,%f,%f).\n",
      argv[0],strut[0].c[0],strut[0].c[1],strut[0].c[2],
              strut[1].c[0],strut[1].c[1],strut[1].c[2]);
    fandr;
  }
  octrope_strut_ends(L,&strutlist[1],strut);
  check_err;
  if ((strut[0].c[0] != 0 || strut[0].c[1] != 0 || strut[0].c[2] != 1 ||
       strut[1].c[0] != 0 || strut[1].c[1] != 1 || strut[1].c[2] != 0) &&
      (strut[0].c[0] != 0 || strut[0].c[1] != 1 || strut[0].c[2] != 0 ||
       strut[1].c[0] != 0 || strut[1].c[1] != 0 || strut[1].c[2] != 1)) {
    fprintf(stderr,
      "%s:Expected a strut (0,0,1)-(0,1,0), got (%f,%f,%f)-(%f,%f,%f).\n",
      argv[0],strut[0].c[0],strut[0].c[1],strut[0].c[2],
              strut[1].c[0],strut[1].c[1],strut[1].c[2]);
    fandr;
  }
  min_strut = 25;
  num_struts = octrope_struts(L,0,0,NULL,0,&min_strut,NULL,0);
  check_err;
  if (num_struts != 0) {
    fprintf(stderr,"%s:Expected no struts, got %d.\n",
      argv[0],num_struts);
    fandr;
  }
  if (fabs(min_strut*min_strut - 2.0) > 1E-12) {
    fprintf(stderr,"%s:Expected min_strut to be sqrt(2), got %17.17f.\n",
      argv[0],min_strut);
    fandr;
  }
  strutlist[0].component[0] = -1;
  num_struts = octrope_struts(L,0,0,strutlist,sl_size,NULL,NULL,0);
  check_err;
  octrope_strut_ends(L,strutlist,strut);
  check_err;
  if ((strut[0].c[0] != 0 || strut[0].c[1] != 0 || strut[0].c[2] != 1 ||
       strut[1].c[0] != 0 || strut[1].c[1] != 1 || strut[1].c[2] != 0) &&
      (strut[0].c[0] != 0 || strut[0].c[1] != 1 || strut[0].c[2] != 0 ||
       strut[1].c[0] != 0 || strut[1].c[1] != 0 || strut[1].c[2] != 1)) {
    fprintf(stderr,
      "%s:Expected a strut (0,0,1)-(0,1,0), got (%f,%f,%f)-(%f,%f,%f).\n",
      argv[0],strut[0].c[0],strut[0].c[1],strut[0].c[2],
              strut[1].c[0],strut[1].c[1],strut[1].c[2]);
    fandr;
  }

  min_rad = 25;
  num_min_rad_locs = 27;
  min_rad = octrope_minrad(L,0,0,min_rad_locs,mr_size,&num_min_rad_locs);
  check_err;
  if (fabs(min_rad - 5.0/14.0) > 1E-12) {
    fprintf(stderr,"%s:Expected MinRad of 5/14, got %17.17f\n",
      argv[0],min_rad);
    fandr;
  }
  if (num_min_rad_locs != 8) {
    fprintf(stderr,"%s:Expected 8 MinRad locations, but got %d.\n",
      argv[0],num_min_rad_locs);
    printf("  Cutoff: 0   minrad: %f  epsilon: 0\n",min_rad);
    for (i=0; i < num_min_rad_locs; i++) {
      printf("  %d : %d : %f\n",min_rad_locs[i].component,min_rad_locs[i].vert,
        min_rad_locs[i].mr);
    }
    fandr;
  }
  min_rad = 25;
  num_min_rad_locs = 27;
  min_rad = octrope_minrad(L,2,0,min_rad_locs,mr_size,&num_min_rad_locs);
  check_err;
  if (fabs(min_rad - 5.0/14.0) > 1E-12) {
    fprintf(stderr,"%s:Expected MinRad of 5/14, got %17.17f\n",
      argv[0],min_rad);
    fandr;
  }
  if (num_min_rad_locs != 12) {
    fprintf(stderr,"%s:Expected 12 MinRad locations, but got %d.\n",
      argv[0],num_min_rad_locs);
    for (i=0; i < num_min_rad_locs; i++) {
      printf("  %d : %d : %f\n",min_rad_locs[i].component,min_rad_locs[i].vert,
        min_rad_locs[i].mr);
    }
    fandr;
  }
  min_rad = 22;
  min_rad = octrope_minrad(L,0,0,NULL,0,NULL);
  check_err;
  if (fabs(min_rad - 5.0/14.0) > 1E-12) {
    fprintf(stderr,"%s:Expected MinRad of 5/14, got %17.17f\n",
      argv[0],min_rad);
    fandr;
  }
  
  for (call_cnt =  0; call_cnt < 1<<7; call_cnt++) {
    /* Set bad values so we know the call worked */
    ropelength = -1;
    thickness = -1;
    curve_length = -1;
    min_rad_locs[0].component = -1;
    min_rad_locs[0].vert = -1;
    min_rad = -1;
    num_min_rad_locs = -1;
    strutlist[0].component[0] = -1;
    strutlist[0].lead_vert[0] = -1;
    min_strut = -1;
    num_struts = -1;
    octrope(L,
            (call_cnt & 1<<6) ? &ropelength       : NULL,
            (call_cnt & 1<<5) ? &thickness        : NULL,
            (call_cnt & 1<<4) ? &curve_length     : NULL,
            (call_cnt & 1<<3) ? &min_rad          : NULL,
            (call_cnt & 1<<1) ? &min_strut        : NULL,
            0,     /* mr_cutoff */
            0.7,   /* mr_epsilon */
            (call_cnt & 1<<2) ? min_rad_locs      : NULL,
              (call_cnt & 1<<2) ? mr_size : 0,
              (call_cnt & 1<<2) ? &num_min_rad_locs : NULL,
            0,     /* strut_cutoff */
            0.7,   /* strut_epsilon */
            (call_cnt & 1<<0) ? strutlist         : NULL,
              (call_cnt & 1<<0) ? sl_size : 0,
              (call_cnt & 1<<0) ? &num_struts       : NULL,
            mem,memsize,
            1      /* lambda        */
            );
    check_err;
    if (call_cnt & 1<<6 && (fabs(ropelength - 112) > 1E-12)) {
      fprintf(stderr,"%s:Expected ropelength 112, found %17.17f.\n",
        argv[0],ropelength);
      fandr;
    }
    if (call_cnt & 1<<5 && (fabs(thickness - 5.0/14.0) > 1E-12)) {
      fprintf(stderr,"%s:Expected thickness 5/14, found %17.17f.\n",
        argv[0],thickness);
      fandr;
    }
    if (call_cnt & 1<<4 && (fabs(curve_length - 40) > 1E-12)) {
      fprintf(stderr,"%s:Expected curve length 40, found %17.17f.\n",
        argv[0],curve_length);
      fandr;
    }
    if (call_cnt & 1<<3 && (fabs(min_rad - 5.0/14.0) > 1E-12)) {
      fprintf(stderr,"%s:Expected MinRad 5/14, found %17.17f.\n",
        argv[0],min_rad);
      fandr;
    }
    if (call_cnt & 1<<2) {
      if (num_min_rad_locs != 8) {
        fprintf(stderr,"%s:Expected 8 minrad locations, found %d:\n",
          argv[0],num_min_rad_locs);
        for (cnt=0; cnt < num_min_rad_locs; cnt++) {
          fprintf(stderr,"  (%d:%d)",
            min_rad_locs[cnt].component, min_rad_locs[cnt].vert);
        }
        fprintf(stderr,"\n");
        fandr;
      }
      if (min_rad_locs[0].component != 0 || min_rad_locs[0].vert != 7 || min_rad_locs[0].svert != 6 ||
	  min_rad_locs[1].component != 0 || min_rad_locs[1].vert != 7 || min_rad_locs[1].svert != 8 ||
	  min_rad_locs[2].component != 0 || min_rad_locs[2].vert != 1 || min_rad_locs[2].svert != 0 ||
	  min_rad_locs[3].component != 0 || min_rad_locs[3].vert != 1 || min_rad_locs[3].svert != 2 ||
	  min_rad_locs[4].component != 0 || min_rad_locs[4].vert != 3 || min_rad_locs[4].svert != 2 ||
	  min_rad_locs[5].component != 0 || min_rad_locs[5].vert != 3 || min_rad_locs[5].svert != 4 ||
	  min_rad_locs[6].component != 0 || min_rad_locs[6].vert != 5 || min_rad_locs[6].svert != 4 ||
	  min_rad_locs[7].component != 0 || min_rad_locs[7].vert != 5 || min_rad_locs[7].svert != 6 )
	{
        fprintf(stderr,"%s:Expected minrads at 0:7(6), 0:7(8), 0:1(0), 0:1(2), 0:3(2), 0:3(4), and 0:5(4), 0:5(6) found\n",
          argv[0]);
        for (cnt=0; cnt < num_min_rad_locs; cnt++) {
          fprintf(stderr,"  (%d:%d(%d))",
		  min_rad_locs[cnt].component, min_rad_locs[cnt].vert,min_rad_locs[cnt].svert);
        }
        fprintf(stderr,"\n");
        fandr;
      }
    }
    if (call_cnt & 1<<1 && (fabs(min_strut*min_strut - 2.0) > 1E-12)) {
      fprintf(stderr,"%s:Expected min strut length sqrt(2), found %17.17f.\n",
        argv[0],min_strut);
      fandr;
    }
    if (call_cnt & 1<<0) {
      if ((strutlist[0].component[0] != 0 || strutlist[0].lead_vert[0] != 2) &&
          (strutlist[0].component[1] != 0 || strutlist[0].lead_vert[1] != 2)) {
        fprintf(stderr,"%s:Expected strut at 0,2, found %d,%d.\n",
          argv[0],strutlist[0].component[0],strutlist[0].lead_vert[0]);
        fandr;
      } 
      if (num_struts != 2) {
        fprintf(stderr,"%s:Expected 2 struts, found %d.\n",
          argv[0],num_struts);
        fandr;
      }
    }
  }

  /* Clean up and quit */
  plc_free(L);
  check_err;
  plc_free(L2);
  check_err;
  plc_free(L3);
  check_err;
  return 0;
}
