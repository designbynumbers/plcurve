/*
 * Sample program to show the use of liboctrope.a
 *
 * $Id: ropelength.c,v 1.17 2007-12-10 21:19:01 cantarel Exp $
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
#ifdef HAVE_FLOAT_H
#include <float.h>
#endif
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#ifdef HAVE_ASSERT_H
#include <assert.h>
#endif

#include "octrope.h"
#include "octrope_internal.h"

#include"argtable2.h"

int VERBOSE = 0;
int PD_VERBOSE = 0;


int main(int argc,char *argv[]) {

  plCurve *L;
  FILE *infile;
  double curve_length;
  double thickness;
  double ropelength;
  double min_rad;
  double min_strut;

  struct arg_file *knotfiles = arg_filen(NULL,NULL,"<file>",1,65536,
    "input file (up to 65536)");
  struct arg_int *d = arg_int0("dD",NULL,"<int>","debug level (>=0)");
  struct arg_lit *dbg = arg_lit0(NULL,"debug","same as -d 4");
  struct arg_int *l = arg_int0("l","levels","<int>","levels in octree (>=0)");
  struct arg_dbl *lam = arg_dbl0("L","lambda","<dbl>","lambda (default=1)");
  struct arg_lit *help = arg_lit0("h","help","print this help and exit");
  struct arg_lit *excesslength = arg_lit0(NULL,"excesslength","excess length for clasp or open knot");
  struct arg_lit *diameter = arg_lit0(NULL,"diameter","computes pointset diameter (slow)");
  struct arg_lit *q = arg_lit0("q","quiet","prints ropelength only");
  struct arg_int *verbosity = arg_int0("v","verbosity","<int>","outputs debugging information");
  struct arg_end *end = arg_end(20);

  void *argtable[] = {lam,l,d,dbg,help,excesslength,diameter,q,verbosity,knotfiles,end};

  int nerrors;
  int filecnt;
  int edges;
  char revision[20] = "$Revision: 1.17 $";
  char *dollar;

  char winning_file[1024] = "No file";
  double winning_ropelength = 1e128, winning_thickness = -1;
  /*int    winning_strutcount = -1, winning_mrstruts  = -1;*/

  dollar = strchr(&revision[1],'$');
  dollar[0] = '\0';
  
  if (arg_nullcheck(argtable) != 0) {
    /* NULL entries detected, allocations must have failed. */
    fprintf(stderr,"%s: insufficient memory\n",argv[0]);
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return 1;
  }
 
  /* Set default lambda */
  lam->dval[0] = 1;
  
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

  /* Now print banner. */

  if (q->count == 0) {
    printf("Ropelength v%s\n",&revision[11]);
    printf("  Measure the (polygonal) ropelength of a polygonal knot.\n");
  }

  if (dbg->count > 0) {
    fprintf(stderr,"%s:Setting debug level to 4.\n",argv[0]);
    octrope_set_debug(4);
  } 
  if (d->count > 0) {
    fprintf(stderr,"%s:Setting debug level to %d.\n",argv[0],
      d->ival[d->count-1]);
    octrope_set_debug(d->ival[d->count-1]);
  }
  if (l->count > 0) {
    if (octrope_debug_level() > 2) {
      fprintf(stderr,"%s:Setting levels to %d\n",argv[0],
        l->ival[l->count-1]);
    }
    octrope_set_levels(l->ival[l->count-1]); 
  }

  if (verbosity->count > 0) {
    VERBOSE = verbosity->ival[0];
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
      fprintf(stderr,"%s:%s\n",argv[0],octrope_error_str); 
      arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
      return 1;
    }

    curve_length = octrope_curvelength(L);
    if (octrope_error_num != 0) {
      fprintf(stderr,"%s:%s\n",argv[0],octrope_error_str); 
      plc_free(L);
      if (octrope_error_num != 0) {
        fprintf(stderr,"%s:%s",argv[0],octrope_error_str); 
      }
      arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
      return 1;
    }
    edges = plc_num_edges(L);
    if (octrope_error_num != 0) {
      fprintf(stderr,"%s:%s\n",argv[0],octrope_error_str); 
      plc_free(L);
      if (octrope_error_num != 0) {
        fprintf(stderr,"%s:%s\n",argv[0],octrope_error_str); 
      }
      arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
      return 1;
    }

    /*************************************************************/

    octrope(L,&ropelength,&thickness,NULL,&min_rad,&min_strut,0,0,NULL,
            0,NULL,0,0,NULL,0,NULL,NULL,0,lam->dval[0]);

    /************************************************************/

    if (octrope_error_num != 0) {
      fprintf(stderr,"%s:%s\n",argv[0],octrope_error_str); 
      plc_free(L);
      if (octrope_error_num != 0) {
        fprintf(stderr,"%s:%s\n",argv[0],octrope_error_str); 
      }
      arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
      return 1;
    }
    if (thickness == 0) {
      if (q->count > 0) {
	printf("inf");
      } else {	
	printf("%s:Ropelength: infinity  Thickness: %f\n",
	       knotfiles->filename[filecnt], thickness);
	printf("%s:Diameter Ropelength: infinity  Diameter Thickness: %f\n",
	       knotfiles->filename[filecnt], thickness);
      }
    } else {
      ropelength = curve_length/thickness;
      if (q->count > 0) {
	printf("%.16g\n",ropelength);
      } else {
	printf("%s:Ropelength: %f  Thickness: %f\n",
	       knotfiles->filename[filecnt], ropelength, thickness);
	printf("%s:Diameter Ropelength: %f  Diameter Thickness: %f\n",
	       knotfiles->filename[filecnt], ropelength/2.0, 2.0*thickness);
      }
    }

    /* If we're checking more than one file, keep track of the best ropelength */

    if (knotfiles->count > 1) {

      if (thickness > 0 && ropelength < winning_ropelength) {

	strncpy(winning_file,knotfiles->filename[filecnt],sizeof(winning_file));
	winning_ropelength = ropelength;
	winning_thickness = thickness;

      }

    }

    /* And now a little extra information, gratis */

    if (excesslength->count > 0) {
      
      /* First collect some information */

      int nopen=0,openguy,n;

      for(n=0;n<L->nc;n++) { 

	if (L->cp[n].open) {nopen++; openguy = n;}

      }

      /* Now break into cases. */

      if (L->nc == 2 && L->cp[0].open && L->cp[1].open) { /* This looks like a clasp */

	/* To determine the excess length for a clasp, there are several steps. */
	
	/* Compute ending tangents */

	plc_vector tans[4];
	bool ok;

	tans[0] = plc_mean_tangent(L,0,1,&ok);
	tans[1] = plc_mean_tangent(L,0,L->cp[0].nv-1,&ok);
	tans[2] = plc_mean_tangent(L,1,1,&ok);
	tans[3] = plc_mean_tangent(L,1,L->cp[1].nv-1,&ok);

	//z = plc_build_vect(0,0,1);
	
	double taus[2],tavg;

	taus[0] = sin(plc_angle(tans[0],tans[1],&ok)/2.0);
	taus[1] = sin(plc_angle(tans[2],tans[3],&ok)/2.0);
	tavg = (taus[0] + taus[1])/2.0;

	if (fabs(taus[0] - taus[1]) > 1e-4) { 

	  printf("ropelength: Tau values for clasp computed from each component (%g, %g)\n"
		 "            don't match. Proceeding with computation, but BEWARE.\n",
		 taus[0],taus[1]);

	}


	double tau;  /* tau = sin theta */

	tau = tavg;

	printf("%s:tau %g ",
	       knotfiles->filename[filecnt], tau);

	/* 3. Compute the inradius. */

	double inradius;

	if (tau > 0.99) { /* This is vertical, but we don't know which way is up. We will assume that z is up. */ 
	
	  plc_vector up,across[2];

	  across[0] = plc_vect_diff(L->cp[0].vt[0],L->cp[0].vt[L->cp[0].nv-1]);
	  across[1] = plc_vect_diff(L->cp[1].vt[1],L->cp[1].vt[L->cp[1].nv-1]);

	  up = plc_normalize_vect(plc_cross_prod(across[0],across[1]),&ok);

	  plc_vector k = {{0,0,1}};

	  if (plc_distance(k,up) < 1e-1) {

	    up = k;

	  } 

	  inradius 
	    = fabs(plc_dot_prod(up,L->cp[0].vt[L->cp[0].nv-1]) + plc_dot_prod(up,L->cp[0].vt[0]) 
		   - plc_dot_prod(up,L->cp[1].vt[0]) - plc_dot_prod(up,L->cp[1].vt[L->cp[1].nv-1]))/4.0;
	   	    	  
	} else { /* We need to actually do the computation. */
	  
	  plc_vector A,B,C,D;
	  plc_vector N1,N2,N3,N4;
	  plc_vector P1,P2,P3,P4;
	  bool Aok,Bok,Cok,Dok;

	  N1 = tans[0]; P1 = L->cp[0].vt[0];
	  N2 = tans[1]; P2 = L->cp[0].vt[L->cp[0].nv-1];
	  N3 = tans[2]; P3 = L->cp[1].vt[0];
	  N4 = tans[3]; P4 = L->cp[1].vt[L->cp[1].nv-1];
	  
	  A = plc_3plane_intersection(N2,P2,N3,P3,N4,P4,&Aok); 
	  B = plc_3plane_intersection(N1,P1,N3,P3,N4,P4,&Bok); 
	  C = plc_3plane_intersection(N1,P1,N2,P2,N4,P4,&Cok); 
	  D = plc_3plane_intersection(N1,P1,N2,P2,N3,P3,&Dok); 

	  assert((Aok && Bok) && (Cok && Dok));

	  inradius = plc_tetrahedron_inradius(A,B,C,D);
	
	  if (VERBOSE > 5) { 

	    printf("A: (%g,%g,%g)\nB: (%g,%g,%g)\nC: (%g,%g,%g)\nD: (%g,%g,%g)\n",
		    plc_M_clist(A), plc_M_clist(B), plc_M_clist(C), plc_M_clist(D));
	  }

	}

	/* We have now computed inradius. */

	printf(" inradius: %g. excesslength: %g\n",inradius,
	       ropelength/2.0 - 4*inradius);
	  

      } else if (nopen == 1) {  /* Looks like an open knot. */
	
	printf("%s:excesslength: %f\n",
	       knotfiles->filename[filecnt], 
	       ropelength - plc_distance(L->cp[openguy].vt[0],
					 L->cp[openguy].vt[L->cp[openguy].nv-1]));
	
      } else {

	printf("ERROR: Excess length requested, but this file appears to be\n"
	       "       a tangle (more than one component is open) which is not\n"
	       "       a simple clasp (there are not two components). Will not\n"
	       "       compute excess length.\n\n");

      }
    }

    if (q->count == 0) {
      
      if (lam->dval[0] > 0) {

	printf("%s:minRad: %f    minStrut: %f\n",
	       knotfiles->filename[filecnt], min_rad, min_strut);

      } else {

	printf("%s:minRad: (not recorded)    minStrut: %f\n",
	       knotfiles->filename[filecnt], min_strut);

      }

      printf("%s:Edges: %d     Average Edgelength: %f\n",
	     knotfiles->filename[filecnt], edges, curve_length/edges);

      if (diameter->count > 0) {

	printf("%s:Pointset Diameter: %g\n",
	       knotfiles->filename[filecnt], plc_pointset_diameter(L));

      }

    }
      
    plc_free(L);
    if (octrope_error_num != 0) {
      fprintf(stderr,"%s:%s\n",argv[0],octrope_error_str); 
      return 1;
    }
  }

  if (knotfiles->count > 1 && winning_ropelength > 0) {

    printf("----------------------------------------------------------\n");
    printf("      Least Ropelength among %d files checked             \n",knotfiles->count);
    printf("----------------------------------------------------------\n");

    printf("%s Rop: %g Thi: %g \n",
	   winning_file,winning_ropelength,winning_thickness);

    printf("\n\n");

  }

  return 0;
}
