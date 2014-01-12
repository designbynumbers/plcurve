/*
 * Test the effects of changing the number of levels in the octree
 *
 * $Id: test_levels.c,v 1.4 2006-04-18 19:14:37 ashted Exp $
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
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#ifdef HAVE_MATH_H
#include <math.h>
#endif
#ifdef HAVE_TIME_H
#include <time.h>
#endif
#include "octrope.h"

int main(int argc,char *argv[]) {

  plCurve *L;
  FILE *infile;
  int cnt;
  void *mem;
  int memsize;
  double rl;
  clock_t time_used;
  int num_edges;

  octrope_set_debug(4);
  memsize = 10000000;
  printf("Reserving 10,000,000 bytes.\n");
  fflush(stdout);
  if ((mem = (void *)malloc(memsize)) == NULL) {
    fprintf(stderr,"Unable to allocate 10,000,000 bytes of memory\n");
    exit(-1);
  }
  printf("Reserved.\n");
  fflush(stdout);

  infile = fopen("hopf.vect","r");
  if (infile == NULL) {
    fprintf(stderr,"test_levels: Couldn't find input file %s.\n",argv[1]);
    exit(-1);
  }

  octrope_error_num = 0;
  L = plc_read(infile,&octrope_error_num,octrope_error_str,80);
  if (octrope_error_num != 0) {
    fprintf(stderr,"test_levels: %s\n",octrope_error_str);
    exit(-1);
  }
  fclose(infile);

  num_edges = plc_num_edges(L);
  if (octrope_error_num != 0) {
    fprintf(stderr,"test_levels: %s\n",octrope_error_str);
    exit(-1);
  }
  printf("test_levels TEST: %d-edge Hopf\n", num_edges);
  fflush(stdout);
  if (L != NULL) {
    printf("test_levels: Loaded link from hopf.vect.\n");
    fflush(stdout);

    for (cnt = 16; cnt >= 1; cnt--) {
      octrope_set_levels(cnt);
      printf("---- Testing with levels = %d ----\n",cnt);
      fflush(stdout);
      time_used = -clock();
      rl = octrope_ropelength(L,mem,memsize,1);
      if (octrope_error_num != 0) {
        fprintf(stderr,"test_levels: %s\n",octrope_error_str);
        exit(-1);
      }
      time_used += clock();
      printf("Ropelength: %7.7f   Time Used:%7.7f\n",rl,
        1.0*time_used/CLOCKS_PER_SEC);
      fflush(stdout);
    }
  }
  plc_free(L);
  if (octrope_error_num != 0) {
    fprintf(stderr,"test_levels: %s\n",octrope_error_str);
    exit(-1);
  }

  printf("End of test.\n");
  fflush(stdout);

  infile = fopen("random.vect","r");
  if (infile == NULL) {
    fprintf(stderr,"test_levels: Couldn't open random.vect.\n");
    exit(-1);
  }

  octrope_error_num = 0;
  L = plc_read(infile,&octrope_error_num,octrope_error_str,80);
  if (octrope_error_num != 0) {
    fprintf(stderr,"test_levels: %s\n",octrope_error_str);
    exit(-1);
  }
  fclose(infile);

  num_edges = plc_num_edges(L);
  if (octrope_error_num != 0) {
    fprintf(stderr,"test_levels: %s\n",octrope_error_str);
    exit(-1);
  }
  printf("test_levels TEST: %d-edge random knot\n", num_edges);
  fflush(stdout);
  if (L != NULL) {
    printf("test_levels: Loaded link from random.vect.\n");
    fflush(stdout);

    for (cnt = 16; cnt >= 1; cnt--) {
      octrope_set_levels(cnt);
      printf("---- Testing with levels = %d ----\n",cnt);
      fflush(stdout);
      time_used = -clock();
      rl = octrope_ropelength(L,mem,memsize,1);
      if (octrope_error_num != 0) {
        fprintf(stderr,"test_levels: %s\n",octrope_error_str);
        exit(-1);
      }
      time_used += clock();
      printf("Ropelength: %7.7f   Time Used:%7.7f\n",rl,
        1.0*time_used/CLOCKS_PER_SEC);
      fflush(stdout);
    }
  }
  plc_free(L);
  if (octrope_error_num != 0) {
    fprintf(stderr,"test_levels: %s\n",octrope_error_str);
    exit(-1);
  }

  printf("End of test.\n");
  fflush(stdout);

  printf("End of run.\n");
  fflush(stdout);

  free(mem);
  exit(0);
}
