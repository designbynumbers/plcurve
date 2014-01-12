/*
 * Sample program to show the use of liboctrope.a
 *
 * $Id: test_minrad.c,v 1.8 2006-04-18 19:14:37 ashted Exp $
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
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#ifdef HAVE_MATH_H
#include <math.h>
#endif
#ifdef HAVE_TIME_H
#include <time.h>
#endif
#include "octrope.h"

#define check_err \
if (octrope_error_num != 0) {\
  fprintf(stderr,"test_minrad:Octrope error %d:%s\n",octrope_error_num,\
    octrope_error_str);\
  if (L != NULL) {\
    plc_free(L);\
    if (octrope_error_num != 0) {\
      fprintf(stderr,"test_minrad:Octrope error %d:%s\n",octrope_error_num,\
        octrope_error_str);\
    }\
  }\
  exit(-1);\
}

void performance_test() {
  int nv = {100};
  bool open = false;
  int cc = 0;
  clock_t    timeused;
  int NUM_TESTS = 1000;
  int i;
  double t_step = (2*3.1415926)/(double)(nv),theta = 0.0;
  double performance;

  plCurve *L;
  plc_vector *lv;

  fprintf(stderr,"Performance testing octrope_minrad.\n");

  L = plc_new(1,&nv,&open,&cc);
  check_err;

  for(i=0;i<nv;i++,theta += t_step) {
    lv = &L->cp[0].vt[i];
    lv->c[0] = cos(theta);
    lv->c[1] = sin(theta);
    lv->c[2] = 0;
  }

  /* Now test its minrad. */

  for(i=0;i<NUM_TESTS;i++) {
    octrope_minrad(L, 0, 0, NULL, 0, NULL);
    check_err;
  }

  timeused = clock();
  performance = (double)(NUM_TESTS)*(double)(nv)*(double)(CLOCKS_PER_SEC)/(double)(timeused);

  fprintf(stderr,"Approximately %g octrope_minrads per second on this machine.\n",performance);
  exit(0);
} 

int main(int argc,char *argv[]) {

  plCurve *L;
  FILE *infile;
  octrope_mrloc *mrl;
  int num_mrls,i;
  double minrad;

  if (argc < 2) {
    printf("usage: testminrad <file.vect>\n");
    exit(2);
  }

  if (!strcmp(argv[1],"-p")) {
    /* Enter performance testing mode. */
    performance_test();
    exit(0);
  }

  infile = fopen(argv[1],"r");

  if (infile == NULL) {
    printf("testminrad: Couldn't find input file %s.\n",argv[1]);
  }

  octrope_error_num = 0;
  L = plc_read(infile,&octrope_error_num,octrope_error_str,80);
  fclose(infile);
  check_err;

  if (L != NULL) {
    printf("testminrad: Loaded link from %s.\n",argv[1]);

    minrad = octrope_minrad(L,0,0,NULL,0,NULL);
    check_err;

    printf("testminrad: minrad = %g.\n",minrad);
    printf("testminrad: Now rerunning with this target * 1.01. \n");

    mrl = calloc(plc_num_edges(L),sizeof(octrope_mrloc));
    check_err;
    
    octrope_minrad(L, 0, 1E-12, mrl, plc_num_edges(L), &num_mrls);
    check_err;

    printf("testminrad: Found %d vertices within 1%% of minimum minrad.\n\n",
	   num_mrls);

    for(i=0;i<num_mrls;i++) {
      printf("Component: %d \t\t Vertex: %d \n", mrl[i].component, 
	     mrl[i].vert);
    }
    
    plc_free(L);
    check_err;

    exit(0);
  } else {
    printf("testminrad: Couldn't load link from %s. \n",argv[1]);
    exit(1);
  }
}
