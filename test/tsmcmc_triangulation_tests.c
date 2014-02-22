/*

  This file runs test code for the tsmcmc_triangulation and plCurve building code.. 

*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#ifdef HAVE_STDIO_H
   #include<stdio.h>
#endif 

#ifdef HAVE_STRING_H
   #include<string.h>
#endif

#ifdef HAVE_STDINT_H
   #include<stdint.h>
#endif

#ifdef HAVE_STDLIB_H
   #include<stdlib.h>
#endif

#ifdef HAVE_GSL_GSL_RNG_H
   #include<gsl/gsl_rng.h>
#endif

#ifdef HAVE_GSL_GSL_RANDIST_H
   #include<gsl/gsl_randist.h>
#endif

#ifdef HAVE_TIME_H
  #include<time.h>
#endif

#ifdef HAVE_MATH_H
  #include<math.h>
#endif

#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

#ifdef HAVE_ASSERT_H
  #include<assert.h>
#endif

#ifdef HAVE_ARGTABLE2_H
  #include<argtable2.h>
#endif

#include<plCurve.h>
//#include <Accelerate/Accelerate.h>
#include <tsmcmc.h>

/* Global variables live here. */

struct arg_int  *seed;        // specify a seed for the random number generator (for debugging)
struct arg_lit  *verbose;
struct arg_lit  *help;
struct arg_lit  *quiet;

struct arg_end  *end;
struct arg_end  *helpend;

gsl_rng *rng; /* The global random number generator */

bool fan_triangulation_tests() {

  printf("\nFan Triangulation Test Suite \n"
	 "-----------------------------\n");

  printf("\tgenerating triangulations for 3-1001 edges...");

  int i;
  tsmcmc_triangulation_t T;

  for(i=3;i<1002;i++) {
    
    T = tsmcmc_fan_triangulation(i);
   
    if (!tsmcmc_triangulation_ok(T)) { 

      printf("%d edge triangulation FAIL\n",i);
      return false;

    }
    
    tsmcmc_triangulation_free(T);

  }

  printf("pass\n");
  printf("-----------------------------\n");

  return true;

}

bool vector_ok(plc_vector a) 
  
/* Tests to see if the vector is ok in floating point terms. */
  
{
  return (isfinite(a.c[0]) && isfinite(a.c[1]) && isfinite(a.c[2]));
}

bool triangle_edgelengths_ok(plc_vector a,plc_vector b,plc_vector c,double A,double B,double C)
{
  if (!vector_ok(a)) { 

    printf("vector a = (%g,%g,%g) is not legitimate floating point\n",plc_M_clist(a));
    return false;

  }

  if (!vector_ok(b)) { 

    printf("vector b = (%g,%g,%g) is not legitimate floating point\n",plc_M_clist(b));
    return false;

  }

  if (!vector_ok(c)) { 

    printf("vector c = (%g,%g,%g) is not legitimate floating point\n",plc_M_clist(c));
    return false;

  }
    
  if (fabs(plc_distance(a,b) - C) > 1e-8) {

    printf("distance from (%g,%g,%g) <-> (%g,%g,%g) = %g != %g\n",
	   plc_M_clist(a), plc_M_clist(b), plc_distance(a,b), C);
    return false;

  }

  if (fabs(plc_distance(b,c) - A) > 1e-8) {

    printf("distance from (%g,%g,%g) <-> (%g,%g,%g) = %g != %g\n",
	   plc_M_clist(b), plc_M_clist(c), plc_distance(b,c), A);
    return false;

  }

  if (fabs(plc_distance(c,a) - B) > 1e-8) {

    printf("distance from (%g,%g,%g) <-> (%g,%g,%g) = %g != %g\n",
	   plc_M_clist(c), plc_M_clist(a), plc_distance(c,a), B);
    return false;

  }

  return true;
}

bool triangle_from_edgelengths_tests() {

  printf("\nTriangle From Edgelengths Test Suite \n"
	 "-------------------------------------\n");

  printf("\tequilateral triangle...");
  plc_vector test[3];
  
  tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],1.0,1.0,1.0);

  if (triangle_edgelengths_ok(test[0],test[1],test[2],1.0,1.0,1.0)) { 

    printf("pass\n");

  } else {

    printf("FAIL\n");
    return false;

  }

  printf("\t3-4-5 triangle...");
    
  tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],3.0,4.0,5.0);

  if (triangle_edgelengths_ok(test[0],test[1],test[2],3.0,4.0,5.0)) { 

    printf("pass\n");

  } else {

    printf("FAIL\n");
    return false;

  }

  printf("\t4-5-3 triangle...");
    
  tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],4.0,5.0,3.0);

  if (triangle_edgelengths_ok(test[0],test[1],test[2],4.0,5.0,3.0)) { 

    printf("pass\n");

  } else {

    printf("FAIL\n");
    return false;

  }

  printf("\t5-3-4 triangle...");
    
  tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],5.0,3.0,4.0);

  if (triangle_edgelengths_ok(test[0],test[1],test[2],5.0,3.0,4.0)) { 

    printf("pass\n");

  } else {

    printf("FAIL\n");
    return false;

  }

  printf("\t1-2-3 triangle...");
 
  tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],1.0,2.0,3.0);

  if (triangle_edgelengths_ok(test[0],test[1],test[2],1.0,2.0,3.0)) { 

    printf("pass\n");

  } else {

    printf("FAIL\n");
    return false;

  }

  printf("\t2-1-3 triangle...");
 
  tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],2.0,1.0,3.0);

  if (triangle_edgelengths_ok(test[0],test[1],test[2],2.0,1.0,3.0)) { 

    printf("pass\n");

  } else {

    printf("FAIL\n");
    return false;

  }

  printf("\t1-3-2 triangle...");
 
  tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],1.0,3.0,2.0);

  if (triangle_edgelengths_ok(test[0],test[1],test[2],1.0,3.0,2.0)) { 

    printf("pass\n");

  } else {

    printf("FAIL\n");
    return false;

  }

  printf("\t2-3-1 triangle...");
 
  tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],2.0,3.0,1.0);

  if (triangle_edgelengths_ok(test[0],test[1],test[2],2.0,3.0,1.0)) { 

    printf("pass\n");

  } else {

    printf("FAIL\n");
    return false;

  }

  printf("\t3-1-2 triangle...");
 
  tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],3.0,1.0,2.0);

  if (triangle_edgelengths_ok(test[0],test[1],test[2],3.0,1.0,2.0)) { 

    printf("pass\n");

  } else {

    printf("FAIL\n");
    return false;

  }

  printf("\t3-2-1 triangle...");
 
  tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],3.0,2.0,1.0);

  if (triangle_edgelengths_ok(test[0],test[1],test[2],3.0,2.0,1.0)) { 

    printf("pass\n");

  } else {

    printf("FAIL\n");
    return false;

  }

  printf("\t0-1-1 triangle...");
  
  tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],0.0,1.0,1.0);

  if (triangle_edgelengths_ok(test[0],test[1],test[2],0.0,1.0,1.0)) { 

    printf("pass\n");

  } else {

    printf("FAIL\n");
    return false;

  }

  printf("\t1-0-1 triangle...");
 
  tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],1.0,0.0,1.0);

  if (triangle_edgelengths_ok(test[0],test[1],test[2],1.0,0.0,1.0)) { 

    printf("pass\n");

  } else {

    printf("FAIL\n");
    return false;

  }

  printf("\t1-1-0 triangle...");
 
  tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],1.0,1.0,0.0);

  if (triangle_edgelengths_ok(test[0],test[1],test[2],1.0,1.0,0.0)) { 

    printf("pass\n");

  } else {

    printf("FAIL\n");
    return false;

  }

  printf("\t0-0-0 triangle...");
 
  tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],0.0,0.0,0.0);

  if (triangle_edgelengths_ok(test[0],test[1],test[2],0.0,0.0,0.0)) { 

    printf("pass\n");

  } else {

    printf("FAIL\n");
    return false;

  }



  printf("\t(-1)-2-3 triangle...");
  
  if (tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],-1.0,2.0,3.0)) { 

    printf("FAIL\n");
    return false; 

  } else {

    printf("pass\n");

  }

  printf("\t1-(-2)-3 triangle...");
  
  if (tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],1.0,-2.0,3.0)) { 

    printf("FAIL\n");
    return false; 

  } else {

    printf("pass\n");

  }

  printf("\t1-2-(-3) triangle...");
  
  if (tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],1.0,2.0,-3.0)) { 

    printf("FAIL\n");
    return false; 

  } else {

    printf("pass\n");

  }

  printf("\t17-2-3 triangle...");
  
  if (tsmcmc_triangle_from_edgelengths(&test[0],&test[1],&test[2],17.0,2.0,3.0)) { 

    printf("FAIL\n");
    return false; 

  } else {

    printf("pass\n");

  }

  printf("------------------------------------\n");

  return true;

}

bool embed_polygon_tests() {

  double PI = 3.141592653589793;

  printf("\nEmbed Polygon Test Suite \n"
	 "-----------------------------\n");

  /* We practice embedding some 4-gons with various dihedrals. */

  double edge_lengths[4] = {1.0,1.0,1.0,1.0};
  double diagonal_lengths[1] = {sqrt(2.0)};
  double dihedral_angles[1] = {PI};

  plCurve *L;
  tsmcmc_triangulation_t T = tsmcmc_fan_triangulation(4);

  printf("\tplanar square - edges (1,1,1,1) - diag (sqrt(2)) - dihedral PI test...");

  L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedral_angles);
  if (!tsmcmc_polygon_embedding_ok(L,T,edge_lengths,diagonal_lengths,dihedral_angles)) { 

    printf("FAIL\n");
    return false;

  }

  if (plc_distance(L->cp[0].vt[0],plc_build_vect(0,0,0)) > 1e-9) { 

    printf("FAIL\n");
    return false;

  }

  if (plc_distance(L->cp[0].vt[1],plc_build_vect(1,0,0)) > 1e-9) { 

    printf("FAIL\n");
    return false;

  }

  if (plc_distance(L->cp[0].vt[2],plc_build_vect(1,1,0)) > 1e-9) { 

    printf("FAIL\n");
    return false;

  }

  if (plc_distance(L->cp[0].vt[3],plc_build_vect(0,1,0)) > 1e-9) { 

    printf("FAIL\n");
    return false;

  }

  plc_free(L);
 
  printf("pass\n");

  printf("\tnondegenerate quadrilateral 1 - edges (1,1,1,1) - diag (sqrt(2)) - various dihedrals...");
  
  double theta = 0;
  int i;

  for(i=0,theta=0;i<100;i++,theta+=2*PI/99.0) {
  
    dihedral_angles[0] = theta;
    L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedral_angles);
    if (!tsmcmc_polygon_embedding_ok(L,T,edge_lengths,diagonal_lengths,dihedral_angles)) { 

      printf("FAIL\n");
      return false;

    }

    if (PI - theta > 0.01) {  /* theta < pi */

      if (L->cp[0].vt[3].c[2] < 0) {

	printf("last vertex (%g,%g,%g) at dihedral angle %g should have POSITIVE z coord\n",
	       plc_M_clist(L->cp[0].vt[3]),theta);
	printf("FAIL.\n");
  
      }

    } else if (theta - PI > 0.01) { /* theta > PI */

       if (L->cp[0].vt[3].c[2] > 0) {

	printf("last vertex (%g,%g,%g) at dihedral angle %g should have NEGATIVE z coord\n",
	       plc_M_clist(L->cp[0].vt[3]),theta);
	printf("FAIL.\n");
  
      }
       
    }

    plc_free(L);

  }
  
  printf("pass\n");

  printf("\tnondegenerate quadrilateral 2 - edges (1,1,1,1) - diagonal (1.0) - various dihedrals...");

  diagonal_lengths[0] = 1.0;

  for(i=0,theta=0;i<100;i++,theta+=2*PI/99.0) {
  
    dihedral_angles[0] = theta;
    L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedral_angles);
    if (!tsmcmc_polygon_embedding_ok(L,T,edge_lengths,diagonal_lengths,dihedral_angles)) { 

      printf("FAIL\n");
      return false;

    }
    plc_free(L);
  
  }
  
  printf("pass\n");

  printf("\tdegenerate quadrilateral 1 - edges (3,1,1,1) - diagonal (2.0) - various dihedrals...");

  edge_lengths[0] = 3.0; edge_lengths[1] = 1.0; edge_lengths[2] = 1.0; edge_lengths[3] = 1.0;
  diagonal_lengths[0] = 2.0;
  
  for(i=0,theta=0;i<100;i++,theta+=2*PI/99.0) {
  
    dihedral_angles[0] = theta;
    L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedral_angles);
    if (!tsmcmc_polygon_embedding_ok(L,T,edge_lengths,diagonal_lengths,dihedral_angles)) { 

      printf("FAIL\n");
      return false;

    }
    plc_free(L);

  }
  
  printf("pass\n");

  printf("\tdegenerate quadrilateral 2 - edges (1,1,1,1) - diagonal (0.0) - various dihedrals...");

  edge_lengths[0] = 1.0; edge_lengths[1] = 1.0; edge_lengths[2] = 1.0; edge_lengths[3] = 1.0;
  diagonal_lengths[0] = 0.0;
  
  for(i=0,theta=0;i<100;i++,theta+=2*PI/99.0) {
  
    dihedral_angles[0] = theta;
    L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedral_angles);
    if (!tsmcmc_polygon_embedding_ok(L,T,edge_lengths,diagonal_lengths,dihedral_angles)) { 

      printf("FAIL\n");
      return false;

    }
    plc_free(L);

  }

  printf("pass\n");

  printf("\tillegal quadrilateral 1 - edges (1,10,1,1) - diagonal (1.0) - various dihedrals...");
  
  edge_lengths[0] = 1.0; edge_lengths[1] = 10.0; edge_lengths[2] = 1.0; edge_lengths[3] = 1.0;
  diagonal_lengths[0] = 1.0;
  
  for(i=0,theta=0;i<100;i++,theta+=2*PI/99.0) {
  
    dihedral_angles[0] = theta;
    L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedral_angles);
    if (L != NULL) { 

      printf("FAIL\n");
      return false;

    }

  }
  
  printf("pass\n");

  printf("\tillegal quadrilateral 2 - edges (1,1,1,1) - diagonal (5.0) - various dihedrals...");
  
  edge_lengths[0] = 1.0; edge_lengths[1] = 1.0; edge_lengths[2] = 1.0; edge_lengths[3] = 1.0;
  diagonal_lengths[0] = 5.0;
  
  for(i=0,theta=0;i<100;i++,theta+=2*PI/99.0) {
  
    dihedral_angles[0] = theta;
    L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedral_angles);
    if (L != NULL) { 

      printf("FAIL\n");
      return false;

    }

  }
  
  printf("pass\n");

  printf("\tillegal quadrilateral 3 - edges (1,1,1,1) - diagonal (-1.5) - various dihedrals...");
  
  edge_lengths[0] = 1.0; edge_lengths[1] = 1.0; edge_lengths[2] = 1.0; edge_lengths[3] = 1.0;
  diagonal_lengths[0] = -1.5;
  
  for(i=0,theta=0;i<100;i++,theta+=2*PI/99.0) {
  
    dihedral_angles[0] = theta;
    L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedral_angles);
    if (L != NULL) { 

      printf("FAIL\n");
      return false;

    }

  }
  
  printf("pass\n");

  tsmcmc_triangulation_free(T);

  double Nedge_lengths[10] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  double Ndiagonal_lengths[7] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  double Ndihedral_angles[7] = {PI,PI,PI,PI,PI,PI,PI};
  tsmcmc_triangulation_t NT = tsmcmc_fan_triangulation(10);

  printf("\tnondegenerate 10-gon - edges all 1.0 - diagonals all 1.0 - dihedrals all PI...");

  L = tsmcmc_embed_polygon(NT,Nedge_lengths,Ndiagonal_lengths,Ndihedral_angles);
  if (!tsmcmc_polygon_embedding_ok(L,NT,Nedge_lengths,Ndiagonal_lengths,Ndihedral_angles)) { 

    printf("FAIL\n");
    return false;

  }
  plc_free(L);
  
  printf("pass\n");

  printf("\tnondegenerate 10-gon - edges all 1.0 - diagonals all 1.0 - dihedrals (PI,PI,PI,PI,0.1,PI/2,0.2)...");

  Ndihedral_angles[4] = 0.1; Ndihedral_angles[5] = PI/2.0; Ndihedral_angles[6] = 0.2;

  L = tsmcmc_embed_polygon(NT,Nedge_lengths,Ndiagonal_lengths,Ndihedral_angles);
  if (!tsmcmc_polygon_embedding_ok(L,NT,Nedge_lengths,Ndiagonal_lengths,Ndihedral_angles)) { 

    printf("FAIL\n");
    return false;

  }
  plc_free(L);
  tsmcmc_triangulation_free(NT);

  printf("pass\n");

  printf("-----------------------------\n");

  return true;

}

bool moment_polytope_step_tests(gsl_rng *rng) {

  printf("\nMoment Polytope Step Test Suite \n"
	 "---------------------------------\n");

  int i,j,n;
  double *edge_lengths, *diagonal_lengths, *dihedral_angles;
  tsmcmc_triangulation_t T;

  for(n=4;n<9;n++) {

    printf("\ttesting for legal diagonals over 10,000 steps on an equilateral fan-triangulated %d-gon...",n);

    T = tsmcmc_fan_triangulation(n);
    tsmcmc_equilateral_ngon(rng,T,&edge_lengths,&diagonal_lengths,&dihedral_angles);    
 
    for(i=0;i<10000;i++) { 

      tsmcmc_moment_polytope_step(rng,T,edge_lengths,diagonal_lengths);
    
      if (!tsmcmc_edgelengths_equilateral(T,edge_lengths) || !tsmcmc_diagonals_ok(T,edge_lengths,diagonal_lengths)) {
	
	printf("FAIL\n");
	return false;
	
      }
      
    }
  
    printf("pass\n");

    free(edge_lengths);
    free(diagonal_lengths);
    free(dihedral_angles);
    tsmcmc_triangulation_free(T);
  
  }

  for(n=4;n<9;n++) {

    double ftc = gsl_ran_flat(rng,0.01,n-0.01);

    printf("\ttesting for legal diagonals over 10,000 steps on a fan-triangulated %d-gon with failure-to-close %g...",n);

    T = tsmcmc_fan_triangulation(n);
    tsmcmc_failure_to_close_ngon(rng,T,ftc,&edge_lengths,&diagonal_lengths,&dihedral_angles);    
 
    for(i=0;i<10000;i++) { 

      tsmcmc_moment_polytope_step(rng,T,edge_lengths,diagonal_lengths);
    
      if (!tsmcmc_edgelengths_obey_ftc(T,ftc,edge_lengths) || !tsmcmc_diagonals_ok(T,edge_lengths,diagonal_lengths)) {
	
	printf("FAIL\n");
	return false;
	
      }
      
    }
  
    printf("pass\n");

    free(edge_lengths);
    free(diagonal_lengths);
    free(dihedral_angles);
    tsmcmc_triangulation_free(T);
  
  }

  printf("\tGenerating Mathematica file 5-gon-diagonals.dat of diagonal data over 10,000 steps...");
  
  FILE *outfile;
  outfile = fopen("5-gon-diagonals.dat","w");

  T = tsmcmc_fan_triangulation(5);
  tsmcmc_equilateral_ngon(rng,T,&edge_lengths,&diagonal_lengths,&dihedral_angles);    
  
  for(i=0;i<10000;i++) { 

    tsmcmc_moment_polytope_step(rng,T,edge_lengths,diagonal_lengths);
    
    if (!tsmcmc_edgelengths_equilateral(T,edge_lengths) || !tsmcmc_diagonals_ok(T,edge_lengths,diagonal_lengths)) {
      
      printf("FAIL\n");
      return false;
      
    }

    fprintf(outfile,"%g %g \n",diagonal_lengths[0],diagonal_lengths[1]);
    
  }

  fclose(outfile);
  
  printf("done\n");
  
  free(edge_lengths);
  free(diagonal_lengths);
  free(dihedral_angles);
  tsmcmc_triangulation_free(T);

  /* ftc = 0.5  */

  printf("\tGenerating Mathematica file 5-gon-ftc-half-diagonals.dat of diagonal data over 10,000 steps...");
  outfile = fopen("5-gon-ftc-half-diagonals.dat","w");

  T = tsmcmc_fan_triangulation(5);
  tsmcmc_failure_to_close_ngon(rng,T,0.5,&edge_lengths,&diagonal_lengths,&dihedral_angles);    
  
  for(i=0;i<10000;i++) { 

    tsmcmc_moment_polytope_step(rng,T,edge_lengths,diagonal_lengths);
    
    if (!tsmcmc_edgelengths_obey_ftc(T,0.5,edge_lengths) || !tsmcmc_diagonals_ok(T,edge_lengths,diagonal_lengths)) {
      
      printf("FAIL\n");
      return false;
      
    }

    fprintf(outfile,"%g %g \n",diagonal_lengths[0],diagonal_lengths[1]);
    
  }

  fclose(outfile);
  
  printf("done\n");
  
  free(edge_lengths);
  free(diagonal_lengths);
  free(dihedral_angles);
  tsmcmc_triangulation_free(T);

  /* ftc = 2.5  */

  printf("\tGenerating Mathematica file 5-gon-ftc-two-half-diagonals.dat of diagonal data over 10,000 steps...");
  
  outfile = fopen("5-gon-ftc-two-half-diagonals.dat","w");

  T = tsmcmc_fan_triangulation(5);
  tsmcmc_failure_to_close_ngon(rng,T,2.5,&edge_lengths,&diagonal_lengths,&dihedral_angles);    
  
  for(i=0;i<10000;i++) { 

    tsmcmc_moment_polytope_step(rng,T,edge_lengths,diagonal_lengths);
    
    if (!tsmcmc_edgelengths_obey_ftc(T,2.5,edge_lengths) || !tsmcmc_diagonals_ok(T,edge_lengths,diagonal_lengths)) {
      
      printf("FAIL\n");
      return false;
      
    }

    fprintf(outfile,"%g %g \n",diagonal_lengths[0],diagonal_lengths[1]);
    
  }

  fclose(outfile);
  
  printf("done\n");
  
  free(edge_lengths);
  free(diagonal_lengths);
  free(dihedral_angles);
  tsmcmc_triangulation_free(T);

  /* ftc = 3.5  */

  printf("\tGenerating Mathematica file 5-gon-ftc-half-diagonals.dat of diagonal data over 10,000 steps...");
  
  outfile = fopen("5-gon-ftc-three-half-diagonals.dat","w");

  T = tsmcmc_fan_triangulation(5);
  tsmcmc_failure_to_close_ngon(rng,T,3.5,&edge_lengths,&diagonal_lengths,&dihedral_angles);    
  
  for(i=0;i<10000;i++) { 

    tsmcmc_moment_polytope_step(rng,T,edge_lengths,diagonal_lengths);
    
    if (!tsmcmc_edgelengths_obey_ftc(T,ftc,edge_lengths) || !tsmcmc_diagonals_ok(T,edge_lengths,diagonal_lengths)) {
      
      printf("FAIL\n");
      return false;
      
    }

    fprintf(outfile,"%g %g \n",diagonal_lengths[0],diagonal_lengths[1]);
    
  }

  fclose(outfile);
  
  printf("done\n");
  
  free(edge_lengths);
  free(diagonal_lengths);
  free(dihedral_angles);
  tsmcmc_triangulation_free(T);


  printf("---------------------------------\n");
  return true;
}

bool confined_moment_polytope_step_tests(gsl_rng *rng) {

  printf("\nConfined Moment Polytope Step Test Suite \n"
	 "-------------------------------------------\n");

  int i,j,n;
  double *edge_lengths, *diagonal_lengths, *dihedral_angles;
  tsmcmc_triangulation_t T;
  double radius;


  for(n=4;n<9;n++) {

    printf("\ttesting for legal diagonals over 10,000 steps on a fan-triangulated %d-gon...",n);

    for(radius = (double)(n)/2.0; radius > 1; radius -= ((double)(n)/2.0 - 1.1)/9) { 
   
      printf("\t\tconfined in sphere of radius %g ... ",radius);
      
      T = tsmcmc_fan_triangulation(n);
      tsmcmc_confined_equilateral_ngon(rng,T,radius,&edge_lengths,&diagonal_lengths,&dihedral_angles);    
      for(i=0;i<10000;i++) { 
	
	tsmcmc_confined_moment_polytope_step(rng,T,radius,edge_lengths,diagonal_lengths);
	
	if (!tsmcmc_edgelengths_equilateral(T,edge_lengths) || 
	    !tsmcmc_confined_diagonals_ok(T,radius,edge_lengths,diagonal_lengths,i)) {
	  
	  printf("FAIL\n");
	  return false;
	  
	}
	
      }
      
      printf("pass\n");
      
      free(edge_lengths);
      free(diagonal_lengths);
      free(dihedral_angles);
      tsmcmc_triangulation_free(T);

    }
  
  }

  double radii[6] = {2.5,2.25,2.0,1.75,1.50,1.25};
  char   names[6][100] = {"5-gon-diagonals-250.dat",
			  "5-gon-diagonals-225.dat",
			  "5-gon-diagonals-200.dat",
			  "5-gon-diagonals-175.dat",
			  "5-gon-diagonals-150.dat",
			  "5-gon-diagonals-125.dat"};

  int ccount;

  for(ccount=0;ccount<6;ccount++) { 

    printf("\tGenerating Mathematica file %s of diagonal data over 10,000 steps...",names[ccount]);
    
    FILE *outfile;
    outfile = fopen(names[ccount],"w");

    T = tsmcmc_fan_triangulation(5);
    tsmcmc_confined_equilateral_ngon(rng,T,radii[ccount],
				     &edge_lengths,&diagonal_lengths,&dihedral_angles);    
  
    for(i=0;i<10000;i++) { 

      tsmcmc_confined_moment_polytope_step(rng,T,radii[ccount],edge_lengths,diagonal_lengths);
    
      if (!tsmcmc_edgelengths_equilateral(T,edge_lengths) || 
	  !tsmcmc_confined_diagonals_ok(T,radii[ccount],edge_lengths,diagonal_lengths,i)) {
      
	printf("FAIL\n");
	return false;
	
      }
      
      fprintf(outfile,"%g %g \n",diagonal_lengths[0],diagonal_lengths[1]);
    
    }

    fclose(outfile);
  
    printf("done\n");
  
    free(edge_lengths);
    free(diagonal_lengths);
    free(dihedral_angles);
    tsmcmc_triangulation_free(T);

  }

  printf("---------------------------------\n");
  return true;
}

bool edgepermute_step_tests(gsl_rng *rng) {

  printf("\nEdgepermute Step Test Suite \n"
	 "---------------------------------\n");
  int i,j,n;
  double *edge_lengths, *diagonal_lengths, *dihedral_angles;
  tsmcmc_triangulation_t T;

  for(n=4;n<9;n++) {

    printf("\ttesting for legal diagonals over 10,000 steps on a fan-triangulated %d-gon...",n);

    T = tsmcmc_fan_triangulation(n);
    tsmcmc_equilateral_ngon(rng,T,&edge_lengths,&diagonal_lengths,&dihedral_angles);    

    for(i=0;i<10000;i++) { 

      tsmcmc_edgepermute_step(rng,T,edge_lengths,diagonal_lengths,dihedral_angles);
    
      if (!tsmcmc_edgelengths_equilateral(T,edge_lengths) || !tsmcmc_diagonals_ok(T,edge_lengths,diagonal_lengths)) {
	
	printf("FAIL\n");
	return false;
	
      }
      
    }
  
    printf("pass\n");

    free(edge_lengths);
    free(diagonal_lengths);
    free(dihedral_angles);
    tsmcmc_triangulation_free(T);
  
  }

  printf("---------------------------------\n");
  return true;
}

bool fan_tri_from_chords(gsl_rng *rng,int n) { 

  int i;
  printf("fan triangulation of %d-gon from chords tests...\n",n);
  
  printf("\tgenerating triangulation from chord system...");
  tsmcmc_triangulation_t T = tsmcmc_fan_triangulation(n);
  printf("done\n");
  
  printf("\ttesting triangulation_ok...");
  if (tsmcmc_triangulation_ok(T)) { 
    printf("pass\n");
  } else {
    printf("FAIL\n");
    return false;
  }

  printf("\tassigning random dihedrals and building...");
  double *dihedrals;
  dihedrals = calloc(n-3,sizeof(double));
  for(i=0;i<n-3;i++) { 
    dihedrals[i] = gsl_ran_flat(rng,0,6.28);
  }

  double *edge_lengths;
  edge_lengths = calloc(n,sizeof(double));
  for(i=0;i<n;i++) { edge_lengths[i] = 1; }
  
  double *diagonal_lengths;
  diagonal_lengths = calloc(n-3,sizeof(double));
  for(i=0;i<n-3;i++) { 
    diagonal_lengths[i] = 1.0;
  }

  plCurve *L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedrals);

  if (tsmcmc_polygon_embedding_ok(L,T,edge_lengths,diagonal_lengths,dihedrals)) { 

    printf("pass\n");

  } else { 

    printf("FAIL\n");
    return false;

  }

  FILE *fant;
  char fan_name[256];
  sprintf(fan_name,"fan_triangulation_%d.csv",n);
  fant = fopen(fan_name,"w");
  tsmcmc_triangulation_print(fant,T);
  fclose(fant);
  printf("\tprinted triangulation to %s\n",fan_name);


  free(edge_lengths);
  free(diagonal_lengths);
  free(dihedrals);
  tsmcmc_triangulation_free(T);

  printf("fan triangulation of %d-gon from chords tests... pass\n\n",n);
  return true;
}

bool spiral_tri_from_chords(gsl_rng *rng,int n) { 

  int i;
  printf("spiral triangulation of %d-gon from chords tests...\n",n);
  
  printf("\tgenerating triangulation from chord system...");
  tsmcmc_triangulation_t T = tsmcmc_spiral_triangulation(n);
  printf("done\n");
  
  printf("\ttesting triangulation_ok...");
  if (tsmcmc_triangulation_ok(T)) { 
    printf("pass\n");
  } else {
    printf("FAIL\n");
    return false;
  }

  printf("\tassigning random dihedrals and building...");
  double *dihedrals;
  dihedrals = calloc(n-3,sizeof(double));
  for(i=0;i<n-3;i++) { 
    dihedrals[i] = gsl_ran_flat(rng,0,6.28);
  }

  double *edge_lengths;
  edge_lengths = calloc(n,sizeof(double));
  for(i=0;i<n;i++) { edge_lengths[i] = 1; }
  
  double *diagonal_lengths;
  diagonal_lengths = calloc(n-3,sizeof(double));
  for(i=0;i<n-3;i++) { 
    diagonal_lengths[i] = 1.0;
  }

  plCurve *L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedrals);

  FILE *looky;
  looky = fopen("looky.vect","w");
  plc_write(looky,L);
  fclose(looky);

  if (tsmcmc_polygon_embedding_ok(L,T,edge_lengths,diagonal_lengths,dihedrals)) { 

    printf("pass\n");

  } else { 

    printf("FAIL\n");
    return false;

  }

  FILE *spiralt;
  char spiral_name[256];
  sprintf(spiral_name,"spiral_triangulation_%d.csv",n);
  spiralt = fopen(spiral_name,"w");
  tsmcmc_triangulation_print(spiralt,T);
  fclose(spiralt);
  printf("\tprinted triangulation to %s\n",spiral_name);


  free(edge_lengths);
  free(diagonal_lengths);
  free(dihedrals);
  tsmcmc_triangulation_free(T);

  printf("spiral triangulation of %d-gon from chords tests... pass\n\n",n);
  return true;
}

bool random_tri_from_chords(gsl_rng *rng,int n) { 

  int i;
  printf("random triangulation of %d-gon from chords tests...\n",n);
  
  printf("\tgenerating triangulation from chord system...");
  tsmcmc_triangulation_t T = tsmcmc_random_triangulation(rng,n);
  printf("done\n");
  
  printf("\ttesting triangulation_ok...");
  if (tsmcmc_triangulation_ok(T)) { 
    printf("pass\n");
  } else {
    printf("FAIL\n");
    return false;
  }

  printf("\tassigning random dihedrals and building...");
  double *dihedrals;
  dihedrals = calloc(n-3,sizeof(double));
  for(i=0;i<n-3;i++) { 
    dihedrals[i] = gsl_ran_flat(rng,0,6.28);
  }

  double *edge_lengths;
  edge_lengths = calloc(n,sizeof(double));
  for(i=0;i<n;i++) { edge_lengths[i] = 1; }
  
  double *diagonal_lengths;
  diagonal_lengths = calloc(n-3,sizeof(double));
  for(i=0;i<n-3;i++) { 
    diagonal_lengths[i] = 1.0;
  }

  plCurve *L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedrals);

  FILE *looky;
  looky = fopen("looky.vect","w");
  plc_write(looky,L);
  fclose(looky);

  if (tsmcmc_polygon_embedding_ok(L,T,edge_lengths,diagonal_lengths,dihedrals)) { 

    printf("pass\n");

  } else { 

    printf("FAIL\n");
    return false;

  }

  FILE *randomt;
  char random_name[256];
  sprintf(random_name,"random_triangulation_%d.csv",n);
  randomt = fopen(random_name,"w");
  tsmcmc_triangulation_print(randomt,T);
  fclose(randomt);
  printf("\tprinted triangulation to %s\n",random_name);


  free(edge_lengths);
  free(diagonal_lengths);
  free(dihedrals);
  tsmcmc_triangulation_free(T);

  printf("random triangulation of %d-gon from chords tests... pass\n\n",n);
  return true;
}

bool teeth_tri_from_chords(gsl_rng *rng,int n) { 

  int i;
  printf("teeth triangulation of %d-gon from chords tests...\n",n);
  
  printf("\tgenerating triangulation from chord system...");
  tsmcmc_triangulation_t T = tsmcmc_teeth_triangulation(n);
  printf("done\n");
  
  printf("\ttesting triangulation_ok...");
  if (tsmcmc_triangulation_ok(T)) { 
    printf("pass\n");
  } else {
    printf("FAIL\n");
    return false;
  }

  printf("\tassigning random dihedrals and building...");
  double *dihedrals;
  dihedrals = calloc(n-3,sizeof(double));
  for(i=0;i<n-3;i++) { 
    dihedrals[i] = gsl_ran_flat(rng,0,6.28);
  }

  double *edge_lengths;
  edge_lengths = calloc(n,sizeof(double));
  for(i=0;i<n;i++) { edge_lengths[i] = 1; }
  
  double *diagonal_lengths;
  diagonal_lengths = calloc(n-3,sizeof(double));
  for(i=0;i<n-3;i++) { 
    diagonal_lengths[i] = 1.0;
  }

  plCurve *L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedrals);

  if (tsmcmc_polygon_embedding_ok(L,T,edge_lengths,diagonal_lengths,dihedrals)) { 

    printf("pass\n");

  } else { 

    printf("FAIL\n");
    return false;

  }

  FILE *teetht;
  char teeth_name[256];
  sprintf(teeth_name,"teeth_triangulation_%d.csv",n);
  teetht = fopen(teeth_name,"w");
  tsmcmc_triangulation_print(teetht,T);
  fclose(teetht);
  printf("\tprinted triangulation to %s\n",teeth_name);


  free(edge_lengths);
  free(diagonal_lengths);
  free(dihedrals);
  tsmcmc_triangulation_free(T);

  printf("teeth triangulation of %d-gon from chords tests... pass\n\n",n);
  return true;
}

bool tsmcmc_generate_triangulation_tests(gsl_rng *rng) { 

  printf("\ntsmcmc_generate_triangulation Test Suite \n"
	 "------------------------------------------------\n");
  
  int i;

  for(i=4;i<17;i++) {

    if (!teeth_tri_from_chords(rng,i)) { 

       printf("------------------------------------------------\n"
	      "FAIL\n\n");
       return false;

    }

  }
  

  for(i=4;i<17;i++) {

    if (!fan_tri_from_chords(rng,i)) { 

       printf("------------------------------------------------\n"
	      "FAIL\n\n");
       return false;

    }

  }

  for(i=4;i<17;i++) {

    if (!spiral_tri_from_chords(rng,i)) { 

       printf("------------------------------------------------\n"
	      "FAIL\n\n");
       return false;

    }

  }

  for(i=4;i<17;i++) {

    if (!random_tri_from_chords(rng,i)) { 

       printf("------------------------------------------------\n"
	      "FAIL\n\n");
       return false;

    }

  }
  

  printf("-----------------------------------------------\n"
	 "tsmcmc_generate_triangulation Test Suite: PASS\n");

  return true;
}

bool tsmcmc_triangulation_string_tests(gsl_rng *rng) {

 printf("\ntsmcmc_triangulation_string Test Suite \n"
	"------------------------------------------------\n");

 printf("strings for fan triangulations 4 to 17 ...\n");

 tsmcmc_triangulation_t T;
 int nedges;

 for(nedges = 4;nedges < 17;nedges++) { 

   T = tsmcmc_fan_triangulation(nedges);
   char *tstring = tsmcmc_triangulation_MathematicaForm(T);
   printf("\t%s\n",tstring);
   free(tstring);
   tsmcmc_triangulation_free(T);

 }

 printf("strings for spiral triangulations 4 to 17 ...\n");

 for(nedges = 4;nedges < 17;nedges++) { 

   T = tsmcmc_spiral_triangulation(nedges);
   char *tstring = tsmcmc_triangulation_MathematicaForm(T);
   printf("\t%s\n",tstring);
   free(tstring);
   tsmcmc_triangulation_free(T);

 }

 printf("strings for random triangulations 4 to 17 ...\n");

 for(nedges = 4;nedges < 17;nedges++) { 

   T = tsmcmc_random_triangulation(rng,nedges);
   char *tstring = tsmcmc_triangulation_MathematicaForm(T);
   printf("\t%s\n",tstring);
   free(tstring);
   tsmcmc_triangulation_free(T);

 }
 
 printf("------------------------------------------------\n"
	"tsmcmc_triangulation_string Test Suite: PASS \n");

 return true;

}

bool tsmcmc_nedge_performance_test(int n,gsl_rng *rng) {

  int TRIALS = 1000;

  tsmcmc_triangulation_t T;
  T = tsmcmc_spiral_triangulation(n);
  double *edge_lengths, *diagonal_lengths, *dihedral_angles;
  tsmcmc_equilateral_ngon(rng,T,&edge_lengths,&diagonal_lengths,&dihedral_angles);

  clock_t start,end;
  double  cpu_seconds_used; 

  printf("\ttesting %d edges...\n",n);
  int i;

  start = clock();
  for(i=0;i<TRIALS;i++) { 
    tsmcmc_dihedrals_step(rng,T,dihedral_angles);
  }
  end = clock();
  cpu_seconds_used = ((double)(end - start))/CLOCKS_PER_SEC;
  printf("\t\tdihedral step = %g millisec\n",1000*cpu_seconds_used/(double)(TRIALS));

  start = clock();
  for(i=0;i<TRIALS;i++) { 
    tsmcmc_moment_polytope_step(rng,T,edge_lengths,diagonal_lengths);
  }
  end = clock();
  cpu_seconds_used = ((double)(end - start))/CLOCKS_PER_SEC;
  printf("\t\tmoment polytope step = %g millisec\n",1000*cpu_seconds_used/(double)(TRIALS));
 
  
  start = clock();
  for(i=0;i<TRIALS;i++) { 
    tsmcmc_edgepermute_step(rng,T,edge_lengths,diagonal_lengths,dihedral_angles);
  }
  end = clock();
  cpu_seconds_used = ((double)(end - start))/CLOCKS_PER_SEC;
  printf("\t\tpermute step = %g millisec\n",1000*cpu_seconds_used/(double)(TRIALS));
  printf("\n");

  tsmcmc_triangulation_free(T);
  free(edge_lengths); free(diagonal_lengths); free(dihedral_angles);

  return true;
}

bool tsmcmc_embed_performance_test(int n,gsl_rng *rng) {

  int TRIALS = 1000;

  tsmcmc_triangulation_t T;
  T = tsmcmc_spiral_triangulation(n);
  double *edge_lengths, *diagonal_lengths, *dihedral_angles;
  tsmcmc_equilateral_ngon(rng,T,&edge_lengths,&diagonal_lengths,&dihedral_angles);

  clock_t start,end;
  double  cpu_seconds_used; 
  plCurve *L;

  printf("\ttesting %d edges...\n",n);
  int i;

  start = clock();
  for(i=0;i<TRIALS;i++) { 
    L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedral_angles);
    plc_free(L);
  }
  end = clock();
  cpu_seconds_used = ((double)(end - start))/CLOCKS_PER_SEC;
  printf("\t\tembed polygon = %g millisec\n",1000*cpu_seconds_used/(double)(TRIALS));

  L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedral_angles);

  double tc;
  start = clock();
  for(i=0;i<TRIALS;i++) { 
    tc = plc_totalcurvature(L,NULL);
  }
  end = clock();
  cpu_seconds_used = ((double)(end - start))/CLOCKS_PER_SEC;
  printf("\t\ttotal curvature = %g millisec\n",1000*cpu_seconds_used/(double)(TRIALS));
  plc_free(L);

  tsmcmc_triangulation_free(T);
  free(edge_lengths); free(diagonal_lengths); free(dihedral_angles);

  printf("\n");

  return true;
}

bool tsmcmc_performance_tests(gsl_rng *rng) {

 printf("\ntsmcmc performance Test Suite \n"
	"------------------------------------------------\n");

 printf("testing relative speed of dihedral, moment-polytope, and permute steps ...\n\n");

 tsmcmc_nedge_performance_test(23,rng);
 tsmcmc_nedge_performance_test(64,rng);

 printf("testing speed of polygon embedding ... \n\n");

 tsmcmc_embed_performance_test(23,rng);
 tsmcmc_embed_performance_test(64,rng);


}

bool tsmcmc_almost_singular_polygons()
{
 printf("\ntsmcmc almost singular polygons Test Suite \n"
	"------------------------------------------------\n");

 printf("testing embedding of polygon with spiral triangulation and...\n\n");
 printf("\t edge lengths: all 1.0 \n"
	"\n"
	"\t diagonal lengths: \n"
	"\t\t 0.8203 1.39919 1.36268 1.34788 0.258789 1.84865 1.2851 \n"
	"\t\t 1.81077 1.14015 0.533886 1.38649 8.17225e-06 1.9399 1.88337 \n"
	"\t\t 1.1497 1.26663 0.815301 1.79551 1.4216 1.70461 1.39776 1.57686 \n"
	"\t\t 0.857435 1.14579 1.39937 0.777705 1.39289 1.18298 1.1659 \n"
	"\t\t 1.39259 0.855614 1.39824 1.7176 1.67816 1.70466 2.45296 1.46251 \n"
	"\t\t 1.38648 3.35318 1.00364 2.01539 1.96247 2.31103 1.67508 1.5026 \n"
	"\t\t 0.778535 1.31492 0.867065 1.90322 0.906939 1.83155 3.54076 2.87585 \n"
	"\t\t 3.57509 1.25227 1.33213 2.73829 3.53427 3.73631 2.22089 5.73045 \n"
	"\n"
	"\t dihedral angles: \n"
	"\t\t3.94435 4.27367 5.71782 1.62429 4.18749 3.02999 2.77182 5.39566\n"
	"\t\t3.15617 5.62749 1.05454 5.12623 1.53307 4.67421 1.86586 3.20549 \n"
	"\t\t4.83249 5.61363 3.75876 5.552 3.78259 5.19444 1.68043 0.935621 \n"
	"\t\t4.06108 0.691943 0.806483 3.79551 4.27225 3.54175 5.33094 1.85065 \n"
	"\t\t3.45624 0.0421008 3.63652 4.54139 3.61467 3.4112 6.11496 2.83757 \n"
	"\t\t1.57716 1.91662 4.24315 4.26672 0.370204 4.98458 5.54396 3.20882 \n"
	"\t\t1.10438 5.97226 4.35798 5.14567 3.05829 2.45555 2.88507 0.891443 \n"
	"\t\t1.01377 1.44409 1.88013 5.77135 5.81233 \n"
	"\n");

 printf("\ttrying to embed this polygon...");

 double edge_lengths[64] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 
			    ,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 
			    ,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

 double diagonal_lengths[61] = {0.8203,1.39919,1.36268,1.34788,0.258789,1.84865,
				1.2851,1.81077,1.14015,0.533886,1.38649,8.17225e-06,
				1.9399,1.88337,1.1497,1.26663,0.815301,1.79551,
				1.4216,1.70461,1.39776,1.57686,0.857435,1.14579,
				1.39937,0.777705,1.39289,1.18298,1.1659,1.39259,
				0.855614,1.39824,1.7176,1.67816,1.70466,2.45296,
				1.46251,1.38648,3.35318,1.00364,2.01539,1.96247,
				2.31103,1.67508,1.5026,0.778535,1.31492,0.867065,
				1.90322,0.906939,1.83155,3.54076,2.87585,3.57509,
				1.25227,1.33213,2.73829,3.53427,3.73631,2.22089,5.73045};
 
 double dihedral_angles[61] = {3.94435,4.27367,5.71782,1.62429,4.18749,3.02999,2.77182,
			       5.39566,3.15617,5.62749,1.05454,5.12623,1.53307,4.67421,
			       1.86586,3.20549,4.83249,5.61363,3.75876,5.552,3.78259,
			       5.19444,1.68043,0.935621,4.06108,0.691943,0.806483,3.79551,
			       4.27225,3.54175,5.33094,1.85065,3.45624,0.0421008,3.63652,
			       4.54139,3.61467,3.4112,6.11496,2.83757,1.57716,1.91662,4.24315,
			       4.26672,0.370204,4.98458,5.54396,3.20882,1.10438,5.97226,4.35798,
			       5.14567,3.05829,2.45555,2.88507,0.891443,1.01377,1.44409,1.88013,
			       5.77135,5.81233};
  
 tsmcmc_triangulation_t T = tsmcmc_spiral_triangulation(64);
 plCurve *L;
 L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedral_angles);

 if (L == NULL) { 

   printf("FAIL\n\n");
   return false;

 }

 printf("pass\n");
 printf("\tcomputing total curvature of this polygon...");

 double kappa = plc_totalcurvature(L,NULL);

 if (isfinite(kappa)) { printf("pass (%g)\n\n",kappa); }
 else { printf("FAIL (%g)\n",kappa); return false; }

 return true;

}

int main(int argc,char *argv[]) {

  int nerrors;

  void *argtable[] = 
    {
      seed = arg_int0(NULL,"seed","<int>","seed for random number generator"),
      verbose = arg_lit0(NULL,"verbose","print debugging information"),
      quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"), 
      help = arg_lit0(NULL,"help","display help message"),
      end = arg_end(20)
    };
  
  void *helptable[] = {help,helpend = arg_end(20)};

  /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("tsmcmc_triangulation_tests: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      fprintf(stderr,"tsmcmc_triangulation_tests compiled " __DATE__ " " __TIME__ "\n");
      arg_print_errors(stdout,end,"tsmcmc_triangulation_tests");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      fprintf(stderr,"tsmcmc_triangulation_tests compiled " __DATE__ " " __TIME__ "\n");
      printf("tsmcmc_triangulation_tests runs unit test code for the triangulation code in the toric symplectic markov chain monte carlo (tsmcmc) package for polygon spaces\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }
  
  clock_t start,end,bigstart,bigend;
  double cpu_time_used;

  const gsl_rng_type * T;
     
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  
  int seedi;
  
  if (seed->count > 0) { seedi = seed->ival[0]; }
  else { seedi = time(0); }
  gsl_rng_set(rng,seedi);
 
  printf("tsmcmc_triangulation test suite\n");   
  printf("with %s random number gen, seeded with %d.\n",gsl_rng_name(rng),seedi);

  if (tsmcmc_almost_singular_polygons() == false) {

    printf("Embedding Almost Singular Polygons Test Suite ... FAIL\n");
    exit(1);

  } else {

    printf("Embedding Almost Singular Polygons Test Suite ... pass\n");
    
  }


  if (tsmcmc_generate_triangulation_tests(rng) == false) {

    printf("Generate Triangulation Test Suite ... FAIL\n");
    exit(1);

  } else {

    printf("Generate Triangulation Test Suite ... pass\n");
    
  }

  if (tsmcmc_triangulation_string_tests(rng) == false) {

    printf("Triagnulation string Test Suite ... FAIL\n");
    exit(1);

  } else {

    printf("Triangulation string Test Suite ... pass\n");
    
  }


  if (fan_triangulation_tests() == false) { 

    printf("Fan Triangulation Test Suite ... FAIL\n");
    exit(1);

  } else {

    printf("Fan Triangulation Test Suite ... pass\n");

  }

  if (triangle_from_edgelengths_tests() == false) { 

    printf("Triangle From Edgelengths Test Suite ... FAIL\n");
    exit(1);
    
  } else {

    printf("Triangle From Edgelengths Test Suite ... pass\n");

  }

  if (embed_polygon_tests() == false) { 

    printf("Embed Polygon Test Suite ... FAIL\n");
    exit(1);
    
  } else {

    printf("Embed Polygon Test Suite ... pass\n");

  }

  if (moment_polytope_step_tests(rng) == false) { 

    printf("Moment Polytope Step Test Suite ... FAIL\n");
    exit(1);
    
  } else {

    printf("Moment Polytope Step Test Suite ... pass\n");

  }
  
  if (edgepermute_step_tests(rng) == false) { 

    printf("Edge Permute Step Test Suite ... FAIL\n");
    exit(1);
    
  } else {

    printf("Edge Permute Step Test Suite ... pass\n");

  }

  if (confined_moment_polytope_step_tests(rng) == false) {
    
    printf("Confined Moment Polytope Step Test Suite ... FAIL\n");

  } else {

    printf("Confined Moment Polytope Step Test Suite ... pass\n");

  }

  if (tsmcmc_performance_tests(rng) == false) {

    printf("Performance Test Suite ... FAIL\n");
    exit(1);
    
  } else {

    printf("Performance Test Suite ... pass\n");
    
  }

  exit(0);

}
  
