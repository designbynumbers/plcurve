/*
 * timing.c : This test program generates 

                  Hopf links,
		  trefoils,
		  random walks

	      at various numbers of edges and tests the speed of the
	      program on these test cases for our algorithm and for
	      the "standard" algorithm.
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
#include <time.h>

#include "octrope.h"

/* Generates a carefully-built Hopf link with <verts_per_comp> vertices in *
 * each ring.                                                              */

plCurve *hopf_link(int total_verts)
{
  int i, nv[2], ccarray[2] = {0,0};
  bool open[2] = {false,false}; 
  double theta, t_step;
  double pi = 3.14159265358979;
  plCurve   *L;
  int verts_per_comp;

  nv[0] = nv[1] = verts_per_comp = total_verts/2;
  
  L = plc_new(2,nv,open,ccarray);
  if (octrope_error_num != 0) {
    fprintf(stderr,"timing: %s\n",octrope_error_str);
    exit(-1);
  }

  for(i=0,t_step = 2.0*pi/(double)(verts_per_comp),theta = t_step/2.0;
      i<verts_per_comp;
      i++,theta += t_step) {

    L->cp[0].vt[i].c[1] = (1/cos(t_step/2.0))*cos(theta);
    L->cp[0].vt[i].c[0] = (1/cos(t_step/2.0))*sin(theta);
    L->cp[0].vt[i].c[2] = 0;

    L->cp[1].vt[i].c[0] = 0;
    L->cp[1].vt[i].c[1] = -(1/cos(t_step/2.0))*cos(theta) + 1.0;
    L->cp[1].vt[i].c[2] = (1/cos(t_step/2.0))*sin(theta);
  }
  return L;
}

/* Generates a carefully-built trefoil knot with N vertices. */

plCurve *trefoil(int verts_per_comp)
{
  int i, nv, cc = 0; 
  bool open = {false};
  double theta, t_step;
  double pi = 3.14159265358979;
  plCurve   *L;

  nv = verts_per_comp;
  
  L = plc_new(1,&nv,&open,&cc);
  if (octrope_error_num != 0) {
    fprintf(stderr,"timing: %s\n",octrope_error_str);
    exit(-1);
  }

  for(i=0,t_step = 4.0*pi/(double)(verts_per_comp),theta = t_step/2.0;
      i<verts_per_comp;
      i++,theta += t_step) {

    L->cp[0].vt[i].c[1] = (1 + 0.66*cos(1.5*theta))*cos(theta);
    L->cp[0].vt[i].c[0] = (1 + 0.66*cos(1.5*theta))*sin(theta);
    L->cp[0].vt[i].c[2] = 0.66*sin(1.5*theta);

  }
  return L;
}

plCurve *random_walk(int verts_per_comp)
{
  int i, nv, cc = 0; 
  bool open = {true};
  plCurve *L;
  plc_vector v;

  nv = verts_per_comp;
  
  L = plc_new(1,&nv,&open,&cc);
  if (octrope_error_num != 0) {
    fprintf(stderr,"timing: %s\n",octrope_error_str);
    exit(-1);
  }

  L->cp[0].vt[0].c[0] = L->cp[0].vt[0].c[1] = L->cp[0].vt[0].c[2] = 0;

  for(i=1;i<nv;i++) {

    v.c[0] = (double)(random() % 1000 - 500)/1000.0;
    v.c[1] = (double)(random() % 1000 - 500)/1000.0;
    v.c[2] = (double)(random() % 1000 - 500)/1000.0;

    L->cp[0].vt[i] = L->cp[0].vt[i-1];
    plc_M_add_vect(L->cp[0].vt[i],v);

  }

  return L;
}
    

double test_time( plCurve *L, int repeat)

     /* Procedure tests L and returns the average run time 
	over "repeat" trials. */

{

  time_t start, end;
  int    i;

  start = clock();
   
  for(i=0;i<repeat;i++) {

    octrope_ropelength(L, NULL, 0, 1.0);

    if (octrope_error_num != 0) {
      fprintf(stderr,"timing: %s\n",octrope_error_str);
      exit(-1);
    }
  }
  
  end = clock();

  return (double)(end - start)/(double)(CLOCKS_PER_SEC*repeat);
}
 
void test_sequence( plCurve *(makelink)(int), 
		    int start_v, int end_v, int step, 
		    int levels, FILE *outfile)

     /* Procedure takes as input a pointer to a function which
	makes a link with a given number of vertices. It then 
	runs speed tests on that family of links, outputting the
	results to outfile. */

{
 int MAX_SAMPLES = 5000;
 
 plCurve *L;
 int n_verts;

 int    n_samples = 0;
 int    verts[MAX_SAMPLES];
 double runtime[MAX_SAMPLES]; 
 
 if (levels) {
   octrope_set_levels(levels);  
 }

 for(n_verts = start_v; n_verts <= end_v; n_verts += step ) {

    L = makelink(n_verts);
 
    if (n_verts < 500) {
      runtime[n_samples] = test_time(L,10);
    } else {
      runtime[n_samples] = test_time(L,5);
    }
    
    verts[n_samples++] = n_verts;
    plc_free(L);

    if (octrope_error_num != 0) {
      fprintf(stderr,"timing: %s\n",octrope_error_str);
      exit(-1);
    }

    fprintf(outfile,"%d %g\n",verts[n_samples-1],runtime[n_samples-1]);
    fflush(outfile);
  }
}


int main(int argc,char *argv[]) {

  /* First, the fast tests. */

  FILE *outfile;

  fprintf(stderr,"timing: Making directory ./timingresults for output\n");
  system("mkdir timingresults");

  /* Runs "fast sequence of tests". */

  fprintf(stderr,"timing: Running fast series of octrope tests.\n");

  outfile = fopen("timingresults/rwalk_fast.dat","w");
  if (outfile == NULL) {fprintf(stderr,"timing: Couldn't open file.\n"); exit(1);}  
  test_sequence(random_walk,10,500,40,0,outfile); 
  test_sequence(random_walk,600,1000,100,0,outfile);
  fclose(outfile);

  outfile = fopen("timingresults/hopflink_fast.dat","w");
  if (outfile == NULL) {fprintf(stderr,"timing: Couldn't open file.\n"); exit(1);}  
  test_sequence(hopf_link,10,500,40,0,outfile);
  test_sequence(hopf_link,500,1000,100,0,outfile);
  fclose(outfile);

  outfile = fopen("timingresults/trefoil_fast.dat","w");
  if (outfile == NULL) {fprintf(stderr,"timing: Couldn't open file.\n"); exit(1);}  
  test_sequence(trefoil,10,500,40,0,outfile);
  test_sequence(trefoil,600,1000,100,0,outfile);
  fclose(outfile);
  
  /* And the "slow sequence of tests". */

  fprintf(stderr,"timing: Running slow series of octope tests.\n");

  outfile = fopen("timingresults/rwalk_slow.dat","w");
  if (outfile == NULL) {fprintf(stderr,"timing: Couldn't open file.\n"); exit(1);}  
  test_sequence(random_walk,10,500,40,1,outfile); 
  test_sequence(random_walk,600,1000,100,1,outfile);
  fclose(outfile);

  outfile = fopen("timingresults/hopflink_slow.dat","w");
  if (outfile == NULL) {fprintf(stderr,"timing: Couldn't open file.\n"); exit(1);}  
  test_sequence(hopf_link,10,500,40,1,outfile);
  test_sequence(hopf_link,600,1000,100,1,outfile);
  fclose(outfile);

  outfile = fopen("timingresults/trefoil_slow.dat","w");
  if (outfile == NULL) {fprintf(stderr,"timing: Couldn't open file.\n"); exit(1);}  
  test_sequence(trefoil,10,500,40,1,outfile);
  test_sequence(trefoil,600,1000,100,1,outfile);
  fclose(outfile);

  /* Now we try to "gnuplot" the results. */

  fprintf(stderr,"timing: Plotting results with gnuplot.\n");

  outfile = fopen("timingresults/gpcommands","w");
  fprintf(outfile,
	  "set logscale xy;\n"
	  "set ylabel \"time (sec)\";\n"
	  "set term postscript;\n"

	  "set xlabel \"number of edges in random walk\";\n"
	  "set output \"timingresults/rwalk_plot.ps\";"
	  "plot \"timingresults/rwalk_slow.dat\" with linespoints "
	  " title \"random walk (levels = 1)\", "
	  " \"timingresults/rwalk_fast.dat\" with linespoints "
	  " title \"random walk (levels = 0)\" ;\n"

	  "set xlabel \"number of edges in hopf link\";\n"
	  "set output \"timingresults/hopflink_plot.ps\";"
	  "plot \"timingresults/hopflink_slow.dat\" with linespoints "
	  " title \"hopf link (levels = 1)\", "
	  " \"timingresults/hopflink_fast.dat\" with linespoints "
	  " title \"hopf link (levels = 0)\";\n"

	  "set xlabel \"number of edges in trefoil\";\n"
	  "set output \"timingresults/trefoil_plot.ps\";"
	  "plot \"timingresults/trefoil_slow.dat\" with linespoints "
	  " title \"trefoil (levels = 1)\", "
	  " \"timingresults/trefoil_fast.dat\" with linespoints "
          " title \"trefoil (levels = 0)\";\n");

  fclose(outfile);
  system("gnuplot timingresults/gpcommands");

  fprintf(stderr,"timing: Results plotted in timing/*.ps. Raw data in timing/*.dat. \n");

  return 0;
}
  
