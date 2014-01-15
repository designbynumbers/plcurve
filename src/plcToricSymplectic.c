#include<plCurve.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <assert.h>
/*  #include <Accelerate/Accelerate.h> */

#include "tsmcmc.h"
#include "geyer_estimators.h"

struct dihedral_struct { 

  int triA[3];
  int triB[3];
  
};  /* Used to represent vertex numbers for two triangles inscribed in 
       a polygon. We will compute dihedrals by computing the (oriented)
       angle between the normals to these triangles. */


double  *tsmcmc_random_dir(gsl_rng *rng,unsigned int n);   
/* Generates a random direction in R^n */

/* Actual code. */

bool     tsmcmc_equilateral_check(plCurve *L)
/* Check edgelengths of polygon to make sure the polygon is still length-1 equilateral. */
{
  double max,min,mean,moment2;
  plc_edgelength_stats(L, &max, &min, &mean, &moment2);
  return (max < 1.0 + 1e-8 && min > 1.0 - 1e-8);
}


bool     tsmcmc_failure_to_close_check(plCurve *L,double ftc)
/* Check edgelengths of polygon to make sure that the first edge has length ftc, the rest are length-1 equilateral. */
{
  int i;
  if (fabs(plc_distance(L->cp[0].vt[0],L->cp[0].vt[1]) - ftc) > 1e-8) { return false; }
  for(i=1;i<L->cp[0].nv;i++) { 
    if (fabs(plc_distance(L->cp[0].vt[i],L->cp[0].vt[i+1]) - 1.0) > 1e-8) { return false; }
  }
  return true;
}

void tsmcmc_regular_ngon(tsmcmc_triangulation_t T,
			 double **edge_lengths,double **diagonal_lengths,double **dihedral_angles)
/* Generate regular planar ngon. */
{
  
  int n = T.nedges;
  double ngon_radius;
  double TWO_PI = 6.2831853071795864769;
  double phi,theta;
  int i;

  plCurve *L;
  bool open = false;
  int  nv = (int)(n);
  int  cc = 0;

  L = plc_new(1,&nv,&open,&cc);

  /* For the regular n-gon with side length 1, we have that the radius */
  /* of the circumscribed circle obeys:

     (1/2) / r = sin (theta/2), 

     or r = (1/2) / sin (theta/2), where theta = TWO_PI/n_edges. 
  */

  theta = TWO_PI/(double)(nv);
  ngon_radius = (1.0/2.0) / sin(theta/2.0);

  for(i=0,phi=0;i<L->cp[0].nv;i++,phi += theta) { 

    L->cp[0].vt[i].c[0] = ngon_radius*cos(phi);
    L->cp[0].vt[i].c[1] = ngon_radius*sin(phi);
    L->cp[0].vt[i].c[2] = 0;

  }

  plc_fix_wrap(L);
  assert(tsmcmc_equilateral_check(L));

  /* Now we convert to edgelengths, diagonals, and dihedrals. */

  *edge_lengths     = calloc(T.nedges,sizeof(double)); assert(edge_lengths != NULL);
  *diagonal_lengths = calloc(T.ndiags,sizeof(double)); assert(diagonal_lengths != NULL);
  *dihedral_angles  = calloc(T.ndiags,sizeof(double)); assert(dihedral_angles != NULL);
  bool *dihedrals_defined = calloc(T.ndiags,sizeof(bool)); assert(dihedrals_defined != NULL);

  tsmcmc_compute_edgelengths(L,T,*edge_lengths);
  tsmcmc_compute_diagonals(L,T,*diagonal_lengths);
  tsmcmc_compute_dihedral_angles(L,T,*dihedral_angles,dihedrals_defined);

  free(dihedrals_defined);
  plc_free(L);
  
}

void tsmcmc_regular_failure_to_close_ngon(tsmcmc_triangulation_t T,
					  double ftc,
					  double **edge_lengths,
					  double **diagonal_lengths,
					  double **dihedral_angles)
/* Generate planar ngon inscribed in circle with first edgelength ftc and all other edgelengths 1. */
{
  
  int n = T.nedges;
  double ngon_radius;
  double TWO_PI = 6.2831853071795864769;
  double phi,theta;
  int i;

  plCurve *L;
  bool open = false;
  int  nv = (int)(n);
  int  cc = 0;

  L = plc_new(1,&nv,&open,&cc);

  /* To solve for the n-gon with (n-1) sides of length 1 and one side of 
     length ftc, we must solve a system of equations for the angle theta
     subtended by the side of length 1, the angle theta' subtended by 
     the side of length ftc, and the radius r of the circle:

     (n-1) theta + theta' = 2 pi

     2 r sin(theta/2) = 1

     2 r sin(theta'/2) = ftc

     with the condition that 0 < theta < 2 pi/3 (assuming that n>=3).
     Solving the first equation for theta', we get 

     theta' = 2 pi - (n-1) theta

     and so we can write the last equation

     2 r sin(pi - (n-1) theta/2) = ftc

     and expanding with sin(a+b) = sin a cos b + sin b cos a,
     we get 

     2 r sin (-(n-1) theta/2) cos pi = ftc

     or 

     2 r sin( (n-1) theta/2) = ftc

     Since ftc = ftc * 1 = ftc 2 r sin(theta/2), we can write
     the desired solution as 

     2 r sin( (n-1) theta/2) = ftc 2 r sin(theta/2)

     and cancelling the 2 r terms, we get a single equation for theta:

     sin( (n-1) theta/2 ) = ftc sin(theta/2).

     Now this isn't algebraically solvable, but we can solve it 
     by Newton's method (or another rootfinding method). If we 
     make the approximation sin x ~ x - x^3/3!, then Mathematica
     gives us the solution

     (2*Sqrt(6)*Sqrt(1 + ftc - n))/Sqrt(1 + ftc - 3*n + 3*Power(n,2) - Power(n,3))

     which we can use as a starting point for the Newton iteration.   
     But be careful! This will generally involve negative square
     roots. So we flip the signs in the square roots.

  */

  double numSqrt,denSqrt;
  numSqrt = (double)(n) - ftc - 1.0;
  denSqrt = -1.0 - ftc + 3.0*n - 3*pow(n,2.0) + pow(n,3.0);

  if (numSqrt > 0 && denSqrt > 0) {

    theta = (2.0*sqrt(6.0)*sqrt(numSqrt))/sqrt(denSqrt);

  } else {

    theta = 0.3; // Just a guess, we've got to start somewhere!

  }
  
  for(i=0;i<10;i++) { 

    theta = theta - (ftc*sin(theta/2.) - sin(((-1 + n)*theta)/2.))/
      ((ftc*cos(theta/2.))/2. - ((-1 + n)*cos(((-1 + n)*theta)/2.))/2.);

  }

  if (!isfinite(theta)) { 

    fprintf(stderr,"tsmcmc_regular_failure_to_close_ngon: Newton iteration failed, resulting in theta = %g.\n",theta);
    exit(1);

  }

  /* Now that we know theta, we can solve for r and check the edgelengths using 
     the equation:

     2 r sin(theta/2) = 1.0

  */

  ngon_radius = 1.0/(2.0 * sin(theta/2));
  
  double theta_prime;
  theta_prime = TWO_PI - (n-1)*theta;

  /* Now we check that the edgelengths work out. */

  if (fabs(2*ngon_radius*sin(theta/2.0) - 1.0) > 1e-8 || fabs(2*ngon_radius*sin(theta_prime/2.0) - ftc) > 1e-8) { 

    fprintf(stderr,"tsmcmc: Newton iteration failed to produce initial n-gon with (n-1) edges of length 1 and one of length %g.\n",
	    ftc);
    exit(1);

  }

  /* We'll build the first two vertices manually. */

  L->cp[0].vt[0] = plc_build_vect(ngon_radius*cos(0),ngon_radius*sin(0),0);
  L->cp[0].vt[1] = plc_build_vect(ngon_radius*cos(theta_prime),ngon_radius*sin(theta_prime),0);
  
  for(i=2,phi=theta_prime+theta;i<L->cp[0].nv;i++,phi += theta) { 

    L->cp[0].vt[i].c[0] = ngon_radius*cos(phi);
    L->cp[0].vt[i].c[1] = ngon_radius*sin(phi);
    L->cp[0].vt[i].c[2] = 0;

  }

  plc_fix_wrap(L);
  assert(tsmcmc_failure_to_close_check(L));

  /* Now we convert to edgelengths, diagonals, and dihedrals. */

  *edge_lengths     = calloc(T.nedges,sizeof(double)); assert(edge_lengths != NULL);
  *diagonal_lengths = calloc(T.ndiags,sizeof(double)); assert(diagonal_lengths != NULL);
  *dihedral_angles  = calloc(T.ndiags,sizeof(double)); assert(dihedral_angles != NULL);
  bool *dihedrals_defined = calloc(T.ndiags,sizeof(bool)); assert(dihedrals_defined != NULL);

  tsmcmc_compute_edgelengths(L,T,*edge_lengths);
  tsmcmc_compute_diagonals(L,T,*diagonal_lengths);
  tsmcmc_compute_dihedral_angles(L,T,*dihedral_angles,dihedrals_defined);

  free(dihedrals_defined);
  plc_free(L);
  
}

void tsmcmc_equilateral_ngon(gsl_rng *rng,tsmcmc_triangulation_t T,
			     double **edge_lengths,double **diagonal_lengths,double **dihedral_angles) {
/* This equilateral unit edge length polygon is used as a default starting point. It is a regular
   n-gon after one moment polytope step, one dihedral step, and one permutation step. */

  tsmcmc_regular_ngon(T,edge_lengths,diagonal_lengths,dihedral_angles);
  tsmcmc_moment_polytope_step(rng,T,*edge_lengths,*diagonal_lengths);
  tsmcmc_dihedrals_step(rng,T,*dihedral_angles);
  tsmcmc_edgepermute_step(rng,T,*edge_lengths,*diagonal_lengths,*dihedral_angles);

}

void tsmcmc_failure_to_close_ngon(gsl_rng *rng,tsmcmc_triangulation_t T,double ftc,double **edge_lengths,double **diagonal_lengths,double **dihedral_angles) {
  /* This is an n-gon with edge 0 having edgelength ftc and the remaining edges with length 1. 
     It has one moment polytope, one dihedral step, and one permutation step applied. */
  
  tsmcmc_regular_failure_to_close_ngon(T,ftc,edge_lengths,diagonal_lengths,dihedral_angles);
  tsmcmc_moment_polytope_step(rng,T,*edge_lengths,*diagonal_lengths);
  tsmcmc_dihedrals_step(rng,T,*dihedral_angles);
  tsmcmc_edgepermute_step(rng,T,*edge_lengths,*diagonal_lengths,*dihedral_angles);

}   

void tsmcmc_folded_triangle(tsmcmc_triangulation_t T,
			    double **edge_lengths,double **diagonal_lengths,double **dihedral_angles)
/* Generate regular planar ngon. */
{
  
  int n = T.nedges;
  int i;

  plCurve *L;
  bool open = false;
  int  nv = (int)(n);
  int  cc = 0;

  L = plc_new(1,&nv,&open,&cc);

  /* 
     The basic idea is that we're going to consider three points:

     (0,0,0) - where the polygon starts and ends
     (1,0,0) - the odd vertices
     (1/2,sqrt(3)/2,0) - the even vertices

     Notice that these three points form the vertices of an
     equilateral triangle. The polygon will just oscillate between
     these two vertices.
  */

  plc_vector verts[2] = {{{1,0,0}},{{0.5,sqrt(3.0)/2.0,0}}};
  
  L->cp[0].vt[0] = plc_build_vect(0,0,0);
 
  for(i=1;i<L->cp[0].nv;i++) { 

    L->cp[0].vt[i] = verts[i%2];
   
  }

  plc_fix_wrap(L);
  assert(tsmcmc_equilateral_check(L));

  /* Now we convert to edgelengths, diagonals, and dihedrals. */

  *edge_lengths     = calloc(T.nedges,sizeof(double)); assert(edge_lengths != NULL);
  *diagonal_lengths = calloc(T.ndiags,sizeof(double)); assert(diagonal_lengths != NULL);
  *dihedral_angles  = calloc(T.ndiags,sizeof(double)); assert(dihedral_angles != NULL);
  bool *dihedrals_defined = calloc(T.ndiags,sizeof(bool)); assert(dihedrals_defined != NULL);

  tsmcmc_compute_edgelengths(L,T,*edge_lengths);
  tsmcmc_compute_diagonals(L,T,*diagonal_lengths);
  tsmcmc_compute_dihedral_angles(L,T,*dihedral_angles,dihedrals_defined);

  free(dihedrals_defined);
  plc_free(L);
  
}

void tsmcmc_confined_equilateral_ngon(gsl_rng *rng,tsmcmc_triangulation_t T,double confinement_radius,
				      double **edge_lengths,double **diagonal_lengths,double **dihedral_angles) {
  
  /* This equilateral unit edge length polygon is used as a default
     starting point. It is a "folded triangle" after one moment
     polytope step and one dihedral step. */

  assert(confinement_radius > 1.0); 

  tsmcmc_folded_triangle(T,edge_lengths,diagonal_lengths,dihedral_angles);
  tsmcmc_confined_moment_polytope_step(rng,T,confinement_radius,*edge_lengths,*diagonal_lengths);
  tsmcmc_dihedrals_step(rng,T,*dihedral_angles);

}

double *tsmcmc_random_dir(gsl_rng *rng,unsigned int n) 
/* Constructs random direction in R^n, using the fact that the Gaussian
   is spherically symmetric. There may be a faster way to do this for small n. */
{
  double *d;
  int i;

  d = malloc(n*sizeof(double));
  assert(d != NULL);

  for(i=0;i<n;i++) { 

    d[i] = gsl_ran_ugaussian(rng);

  }

  return d;
}

double   tsmcmc_chord_length(int chord,chordtype_t type,
	     		     double *diagonal_lengths,double *edge_lengths);

double   tsmcmc_variation(int chord,chordtype_t type,double *dv) 
{
  if (type == edge) { return 0; }
  else { return dv[chord]; }
}

void     tsmcmc_apply_triangle_inequality(int chordA,chordtype_t typeA,
					  int chordB,chordtype_t typeB,
					  int chordC,chordtype_t typeC,
					  double *edge_lengths,double *diagonal_lengths,double *dv,
					  double *tupper,double *tlower) 
  /* 
    Apply the triangle inequality:

    (a + t da) + (b + t db) - (c + t dc) >= 0
    
    or 
    
    (a + b - c) + t (da + db - dc) >= 0
    
    or 
    
    t (da + db - dc) >= -(a + b - c).
    
    to tupper, tlower. */

{
  double num, denom;

  num = -tsmcmc_chord_length(chordA,typeA,diagonal_lengths,edge_lengths) 
        -tsmcmc_chord_length(chordB,typeB,diagonal_lengths,edge_lengths) 
        +tsmcmc_chord_length(chordC,typeC,diagonal_lengths,edge_lengths);

  denom = tsmcmc_variation(chordA,typeA,dv) + tsmcmc_variation(chordB,typeB,dv) - tsmcmc_variation(chordC,typeC,dv);

  if (!isfinite(num) || !isfinite(denom)) { 

    printf("tsmcmc_apply_triangle_inequality: num = %g or denom = %g is not finite.\n"
	   "-----------------------------------------------------------------------\n"
	   "chordA = %d, chordB = %d, chordC = %d \n"
	   "typeA = %s, typeB = %s, typeC = %s\n"
	   "variation = %g, variation = %g, variation = %g\n"
	   "-----------------------------------------------------------------------\n",
	   num,denom,
	   chordA,chordB,chordC,
	   (typeA == edge) ? "edge" : "diagonal",
	   (typeB == edge) ? "edge" : "diagonal",
	   (typeC == edge) ? "edge" : "diagonal",
	   tsmcmc_variation(chordA,typeA,dv),
	   tsmcmc_variation(chordB,typeB,dv),
	   tsmcmc_variation(chordC,typeC,dv));
    exit(1);

  }	   

  if (fabs(denom) < 1e-8) { /* This is a degenerate case, just return. */ 

    return;

  } else if (denom > 0) { /* We have t >= num/denom, and since num <= 0, we know num/denom < 0. */
                          /* This means the solution only affects tlower. */

    if (num/denom > *tlower) { *tlower = num/denom; }
    return;

  } else { /* We know denom < 0, so the inequality is t <= num/denom, and num/denom > 0. */
           /* This case only affects tupper. */

    if (num/denom < *tupper) { *tupper = num/denom; }
    return;

  }

}    

void     tsmcmc_moment_polytope_worker(tsmcmc_triangulation_t T,double *edge_lengths,
				       double *diagonal_lengths, double *dv,
				       int tri,double *tupper,double *tlower)

/* Apply the triangle inequalities for this triangle to the current bounds for t,
   tightening them if needed, using the current edge_lengths and diagonal_lengths,
   and step direction dv. This triangle is composed of three chords we'll call 
   X, Y, and Z. */

{
  int X,Y,Z;
  chordtype_t Xtype, Ytype, Ztype;

  X = T.triangles[tri].parent_chord; Xtype = T.triangles[tri].parent_type;
  Y = T.triangles[tri].daughter_chord[0]; Ytype = T.triangles[tri].daughter_type[0];
  Z = T.triangles[tri].daughter_chord[1]; Ztype = T.triangles[tri].daughter_type[1];

  tsmcmc_apply_triangle_inequality(X,Xtype,Y,Ytype,Z,Ztype,edge_lengths,diagonal_lengths,dv,tupper,tlower);
  tsmcmc_apply_triangle_inequality(Y,Ytype,Z,Ztype,X,Xtype,edge_lengths,diagonal_lengths,dv,tupper,tlower);
  tsmcmc_apply_triangle_inequality(Z,Ztype,X,Xtype,Y,Ytype,edge_lengths,diagonal_lengths,dv,tupper,tlower);

  /* Now recurse */
  
  int i;

  for(i=0;i<2;i++) { 

    if (T.triangles[tri].daughter_tri[i] != -1) { 

      tsmcmc_moment_polytope_worker(T,edge_lengths,diagonal_lengths,dv,
				    T.triangles[tri].daughter_tri[i],tupper,tlower);
      
    }
    
  }

}

void     tsmcmc_moment_polytope_step(gsl_rng *rng,tsmcmc_triangulation_t T,
				     double *edge_lengths,double *diagonal_lengths)
/*  Make a hit-and-run step in the moment polytope using the triangulation T, which alters diagonal_lengths. */
{
  
  double *dv;
  dv = tsmcmc_random_dir(rng,T.ndiags);
  double tupper = 1e12, tlower = -1e12;
  
  tsmcmc_moment_polytope_worker(T,edge_lengths,diagonal_lengths,dv,0,&tupper,&tlower);
  double t; 
  t = gsl_ran_flat(rng,tlower,tupper);
  
  int i;
  //cblas_daxpy(T.ndiags,t,dv,1,diagonal_lengths,1);
  for(i=0;i<T.ndiags;i++) {
    diagonal_lengths[i] += t*dv[i];
  }
  
  bool ok_flag=true;
  for(i=0;i<T.ndiags;i++) { if (!isfinite(diagonal_lengths[i])) { ok_flag = false; } }

  if (!ok_flag) { 

    printf("tsmcmc_moment_polytope_step: A diagonal is not a finite value. \n"
	   "-------------------------------------------------------------- \n"
	   "tlower = %g   t = %g   tupper = %g \n"
	   "dv = (",tlower,t,tupper);

    for(i=0;i<T.ndiags-1;i++) { printf("%g,",dv[i]); }; printf("%g)\n",dv[i]);

    printf("diagonal_lengths = (");
    for(i=0;i<T.ndiags-1;i++) { printf("%g,",diagonal_lengths[i]); }; printf("%g)\n",diagonal_lengths[i]);
    printf("------------------------------------------------------------- \n");
    exit(1);

  }

  free(dv);

}

void     tsmcmc_confined_moment_polytope_step(gsl_rng *rng,tsmcmc_triangulation_t T,
					      double confinement_radius,
					      double *edge_lengths,double *diagonal_lengths)
/*  Make a hit-and-run step in the moment polytope using the triangulation T, which alters diagonal_lengths. */
{

  double *dv;
  dv = tsmcmc_random_dir(rng,T.ndiags);
  double tupper = 1e12, tlower = -1e12;

  tsmcmc_moment_polytope_worker(T,edge_lengths,diagonal_lengths,dv,0,&tupper,&tlower);
  double t; 

  /* We've now set the diagonal_lengths array so that we stay in the
     moment polytope. We need to further confine tlower and tupper so 
     that we stay in the confinement sphere. To do that, we need 

     diagonal_lengths[i] + t*dv[i] <= confinement_radius 

     or 

     t*dv[i] <= confinement_radius - diagonal_lengths[i]

     In fact, we need a little more, since the diagonals are bounded
     below by zero. But these inequalities were already taken care of
     in the moment polytope worker. Whether this inequality
     (potentially) affects tupper or tlower at each stage is
     determined by the sign of dv[i].

  */

  int i;

  for(i=0;i<T.ndiags;i++) {

    if (fabs(dv[i]) > 1e-8) { /* It's entirely possible to get a dv equal to zero, but it shouldn't affect tupper or tlower. */

      if (dv[i] > 0) { 
	
	tupper = ((confinement_radius-diagonal_lengths[i])/dv[i] < tupper) ? 
	  (confinement_radius-diagonal_lengths[i])/dv[i] : tupper;
	
      } else { 
	
	tlower = ((confinement_radius-diagonal_lengths[i])/dv[i] > tlower) ? 
	  (confinement_radius-diagonal_lengths[i])/dv[i] : tlower;
	
      }

    }

  }

  t = gsl_ran_flat(rng,tlower,tupper);

  for(i=0;i<T.ndiags;i++) {
    diagonal_lengths[i] += t*dv[i];
  }
  //cblas_daxpy(T.ndiags,t,dv,1,diagonal_lengths,1);

  bool ok_flag=true;
  for(i=0;i<T.ndiags;i++) { if (!isfinite(diagonal_lengths[i])) { ok_flag = false; } }

  if (!ok_flag) { 

    printf("tsmcmc_confined_moment_polytope_step: A diagonal is not a finite value. \n"
	   "-------------------------------------------------------------- \n"
	   "tlower = %g   t = %g   tupper = %g \n"
	   "dv = (",tlower,t,tupper);

    for(i=0;i<T.ndiags-1;i++) { printf("%g,",dv[i]); }; printf("%g)\n",dv[i]);

    printf("diagonal_lengths = (");
    for(i=0;i<T.ndiags-1;i++) { printf("%g,",diagonal_lengths[i]); }; printf("%g)\n",diagonal_lengths[i]);
    printf("------------------------------------------------------------- \n");
    exit(1);
    
  }

  free(dv);

}

bool   tsmcmc_diagonals_ok_worker(tsmcmc_triangulation_t T,double *edge_lengths,double *diagonal_lengths,int tri)

/* Check triangle inequalities for this triangle; recurse. */

{
  double a,b,c;

  a = tsmcmc_chord_length(T.triangles[tri].parent_chord,
			  T.triangles[tri].parent_type,
			  diagonal_lengths,edge_lengths);
 
  b = tsmcmc_chord_length(T.triangles[tri].daughter_chord[0],
			  T.triangles[tri].daughter_type[0],
			  diagonal_lengths,edge_lengths);

  c = tsmcmc_chord_length(T.triangles[tri].daughter_chord[1],
			  T.triangles[tri].daughter_type[1],
			  diagonal_lengths,edge_lengths);

  /* Now check triangle inequalities */

  if (!isfinite(a) || !isfinite(b) || !isfinite(c) || a + b < c || b + c < a || c + a < b) { 

    printf("tsmcmc_diagonal_lengths_ok: Lengths are bad (or not numbers) for tri %d \n"
	   "------------------------------------------------------------------\n"
	   "side a: side length %g, chord_type: %s, chord_number: %d \n"
	   "side b: side length %g, chord_type: %s, chord_number: %d \n"
	   "side c: side length %g, chord_type: %s, chord_number: %d \n"
	   "------------------------------------------------------------------\n",
	   tri,
	   a, (T.triangles[tri].parent_type == edge) ? "edge    " : "diagonal",
	   T.triangles[tri].parent_chord,
	   
	   b, (T.triangles[tri].daughter_type[0] == edge) ? "edge    " : "diagonal",
	   T.triangles[tri].daughter_chord[0],

	   c, (T.triangles[tri].daughter_type[1] == edge) ? "edge    " : "diagonal",
	   T.triangles[tri].daughter_chord[1]);

    return false; }

  int i;

  for(i=0;i<2;i++) { 

    if (T.triangles[tri].daughter_tri[i] != -1) { 

      if (!tsmcmc_diagonals_ok_worker(T,edge_lengths,diagonal_lengths,T.triangles[tri].daughter_tri[i])) { 

	return false; 

      }

    }

  }

  return true;

}

bool   tsmcmc_diagonals_ok(tsmcmc_triangulation_t T,double *edge_lengths,double *diagonal_lengths)
/* Check to make sure that the diagonals obey the triangle inequalities from the given triangulation */

{
  bool ok_flag = true;

  /* First, do some consistency checks on the triangulation values */

  assert(T.nedges > 0 && T.nedges < 1e10);
  assert(T.ndiags >= 0 && T.ndiags < 1e10);
  assert(T.ndiags == T.nedges - 3);

  int i;
  /* Next, we scan to check the edge lengths and diagonal lengths for garbage values */

  for(i=0;i<T.nedges;i++) { if (!isfinite(edge_lengths[i])) { ok_flag = false; } }

  if (ok_flag) { 

    for(i=0;i<T.ndiags;i++) { if (!isfinite(diagonal_lengths[i])) { ok_flag = false; } }

    if (ok_flag) { 

      if (!tsmcmc_diagonals_ok_worker(T,edge_lengths,diagonal_lengths,0)) { ok_flag = false; }

      if (ok_flag) { return true; }

    }

  }

  /* Something is wrong; drop a (hopefully informative) error message. */

  printf("tsmcmc_diagonals_ok: Triangulation has %d edges and %d == (%d - 3) diagonals \n",
	 T.nedges,T.ndiags,T.nedges);

  printf("                     Edge lengths: ");
  for(i=0;i<T.nedges;i++) { printf("%g ",edge_lengths[i]); } 
  printf("\n");

  printf("                     Diagonal lengths: ");
  for(i=0;i<T.ndiags;i++) { printf("%g ",diagonal_lengths[i]); } 
  printf("\n");

  printf("\n");

  return false;
}

bool   tsmcmc_confined_diagonals_ok(tsmcmc_triangulation_t T,double confinement_radius,
				    double *edge_lengths,double *diagonal_lengths,int step)
/* Check to make sure that the diagonals obey the triangle inequalities from the given triangulation */

{
  if (tsmcmc_diagonals_ok_worker(T,edge_lengths,diagonal_lengths,0) == false) { 

    return false;

  } 

  int i;
  for(i=0;i<T.ndiags;i++) { 

    if (diagonal_lengths[i] > confinement_radius) { 

      printf("tsmcmc_triangulation_ok: Failed confinement check at step %d \n"
	     "-----------------------------------------------------\n"
	     "diagonal_length %d == %g > confinement_radius == %g  \n"
	     "-----------------------------------------------------\n",
	     step,i,diagonal_lengths[i],confinement_radius);

      return false;

    }

  }

  return true;
}

bool     tsmcmc_edgelengths_equilateral(tsmcmc_triangulation_t T,double *edge_lengths)
/* Check to make sure all edges are length 1.0 */
{
  int i;
  for(i=0;i<T.nedges;i++) {

    if (fabs(edge_lengths[i] - 1.0) > 1e-8) { 

      printf("tsmcmc_edgelengths_equilateral: edge %d has length %g != 1.0\n",
	     i,edge_lengths[i]);

      return false;

    }

  }

  return true;

}

bool     tsmcmc_edgelengths_obey_ftc(tsmcmc_triangulation_t T,double ftc, double *edge_lengths)
/* Check to see that edge 0 has length ftc and all other edges have length 1. */
{
  int i;
  if (fabs(ftc - edge_lengths[0]) > 1e-8) { 
    
    printf("tsmcmc_edgelengths_obey_ftc: edge 0 has length %g != ftc (%g)\n",edge_lengths[0],ftc);
    return false; 

  }

  for(i=1;i<T.nedges;i++) { 

    if (fabs(edge_lengths[i] - 1.0) > 1e-8) { 

      printf("tsmcmc_edgelengths_obey_ftc: edge %i has length %g != 1.0\n",i,edge_lengths[i]);
      return false;

    }

  }

  return true;
}


void     tsmcmc_dihedrals_step(gsl_rng *rng,tsmcmc_triangulation_t T,double *dihedral_angles)
/* Reset dihedral angles uniformly. */
{
  double PI = 3.141592653589793;
  int i;

  for(i=0;i<T.ndiags;i++) { 

    dihedral_angles[i] = gsl_ran_flat(rng,0,2.0*PI);

  }
}
    

void tsmcmc_edgepermute_step(gsl_rng *rng,tsmcmc_triangulation_t T,
			     double *edge_lengths, double *diagonal_lengths, 
			     double *dihedral_angles)
/* Generate a polygon, permute edges, and then recompute diagonal_lengths, dihedral_angles, edge_lengths. */

{
  plCurve *L;
  plc_vector *edgeset;

  edgeset = calloc(T.nedges,sizeof(plc_vector));
  assert(edgeset != NULL);

  L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedral_angles);

  int i;
  for(i=0;i<T.nedges;i++) { edgeset[i] = plc_vect_diff(L->cp[0].vt[i+1],L->cp[0].vt[i]); }

  if (fabs(edge_lengths[0] - edge_lengths[1]) > 1e-8) { // Edgelength 0 is different, use pointer arithmetic to shuffle only part of the buffer of edges

    gsl_ran_shuffle(rng,edgeset+1,T.nedges-1,sizeof(plc_vector));

  } else { // We assume that the polygon is equilateral. 

    gsl_ran_shuffle(rng,edgeset,T.nedges,sizeof(plc_vector));

  }

  plCurve *newL;

  /* Recreate a new vertex set from these edges */

  int  nv = L->cp[0].nv;
  bool open = L->cp[0].open;
  int  cc = 0;

  newL = plc_new(1,&nv,&open,&cc);
  assert(newL != NULL);

  newL->cp[0].vt[0] = plc_build_vect(0,0,0);
  for(i=1;i<T.nedges;i++) { newL->cp[0].vt[i] = plc_vect_sum(edgeset[i],newL->cp[0].vt[i-1]); }
  plc_fix_wrap(newL);

  free(edgeset);
  plc_free(L);

  /* Now recompute diagonals, dihedrals, and edgelengths. */

  bool *dihedral_defined = calloc(T.ndiags,sizeof(bool));
  
  tsmcmc_compute_diagonals(newL,T,diagonal_lengths);
  tsmcmc_compute_dihedral_angles(newL,T,dihedral_angles,dihedral_defined);
  tsmcmc_compute_edgelengths(newL,T,edge_lengths);

  plc_free(newL);
  free(dihedral_defined);
}

tsmcmc_step_t  tsmcmc_dihedral_diagonal_step(gsl_rng *rng, tsmcmc_triangulation_t T,
					     double *edge_lengths, double *diagonal_lengths,double *dihedral_angles, 
					     double beta)
/* Randomly chooses dihedral or diagonal step. */
{
  assert(beta <= 1.0 && beta >= 0.0);
  double b = gsl_ran_flat(rng,0.0,1.0);
  
  if (b < beta) { 

    tsmcmc_moment_polytope_step(rng,T,edge_lengths,diagonal_lengths);
    return moment_polytope;
    
  } else {

    tsmcmc_dihedrals_step(rng,T,dihedral_angles);
    return dihedral;

  } 

}

tsmcmc_step_t tsmcmc_dihedral_diagonal_permute_step(gsl_rng *rng, tsmcmc_triangulation_t T,
						    double *edge_lengths, double *diagonal_lengths,
						    double *dihedral_angles, 
						    double beta, double delta, int moment_polytope_repeat)
/* Randomly chooses dihedral, diagonal, or permute step. */
{
  assert(beta <= 1.0 && beta >= 0.0); assert(delta <= 1.0 && delta >= 0.0); 
  double b = gsl_ran_flat(rng,0.0,1.0);
  double d = gsl_ran_flat(rng,0.0,1.0);

  if (d < delta) { 

    tsmcmc_edgepermute_step(rng,T,edge_lengths,diagonal_lengths,dihedral_angles);
    return permute;

  } else {

    if (b < beta) { 

      int i;
      for(i=0;i<moment_polytope_repeat;i++) { 
	tsmcmc_moment_polytope_step(rng,T,edge_lengths,diagonal_lengths);
      }
      return moment_polytope;

    } else {

      tsmcmc_dihedrals_step(rng,T,dihedral_angles);
      return dihedral;
      
    } 

  }

}

tsmcmc_step_t tsmcmc_confined_dihedral_diagonal_step(gsl_rng *rng, tsmcmc_triangulation_t T,
						     double confinement_radius,
						     double *edge_lengths, double *diagonal_lengths,
						     double *dihedral_angles, 
						     double beta, int moment_polytope_repeat)
/* Randomly chooses dihedral, diagonal, or permute step. */
{
  assert(beta <= 1.0 && beta >= 0.0);
  double b = gsl_ran_flat(rng,0.0,1.0);
  
  if (b < beta) { 
    
    int i;
    for(i=0;i<moment_polytope_repeat;i++) {
      tsmcmc_confined_moment_polytope_step(rng,T,confinement_radius,edge_lengths,diagonal_lengths);
    }
    return moment_polytope;

  } else {

    tsmcmc_dihedrals_step(rng,T,dihedral_angles);
    return dihedral;

  } 

}


char    *tsmcmc_run_stats_MathematicaForm(tsmcmc_run_stats run_stats)
{

  char *ret = calloc(2048,sizeof(char));
  assert(ret != NULL);

  sprintf(ret,
	  "{\"max lagged covariance used\" -> %d, \n"
	  " \"lagged covariances available\" -> %d, \n"
	  " \"dihedral steps\" -> %d, \n"
	  " \"moment polytope steps\" -> %d, "
	  " \"permute steps\" -> %d, "
	  " \"markov chain time\" -> %g, "
	  " \"error estimate time\" -> %g }",
	  run_stats.max_lagged_covariance_used,
	  run_stats.lagged_covariances_available,
	  run_stats.dihedral_steps,
	  run_stats.mp_steps,
	  run_stats.permute_steps,
	  run_stats.total_seconds,
	  run_stats.geyer_ips_seconds);

  return ret;
}

char    *tsmcmc_run_params_MathematicaForm(tsmcmc_run_parameters run_params)

{
  char *ret = calloc(1024,sizeof(char));
  assert(ret != NULL);

  sprintf(ret,
	  "{ \"burn in steps\" -> %d, \n"
	  "  \"delta (permute fraction)\" -> %g, \n"
	  "  \"beta (mp fraction)\" -> %g, \n"
	  "  \"moment polytope repeat\" -> %d, \n"
	  "  \"log file\"->\"%s\" }",
	  run_params.burn_in,
	  run_params.delta,
	  run_params.beta,
	  run_params.moment_polytope_repeat,
	  run_params.logfile_name);

  return ret;
}

char    *tsmcmc_configuration_MathematicaForm(tsmcmc_triangulation_t T,double *edge_lengths,double *diagonal_lengths,double *dihedral_angles)
{

  char *ret = calloc(2048,sizeof(char));
  assert(ret != NULL);

  sprintf(ret,
	  "{\"nedges\"->%d, \"ntri\"->%d, \"ndiags\"->%d,\"edge_lengths\"->{",
	  T.nedges,T.ntri,T.ndiags);

  int i;
  for(i=0;i<T.nedges-1;i++) { sprintf(ret + strlen(ret),"%g,",edge_lengths[i]); } 
  sprintf(ret+strlen(ret),"%g},",edge_lengths[i]);

  sprintf(ret+strlen(ret),"\"diagonal_lengths\"->{");
  for(i=0;i<T.ndiags-1;i++) { sprintf(ret + strlen(ret),"%g,",diagonal_lengths[i]); }
  sprintf(ret+strlen(ret),"%g},",diagonal_lengths[i]);

  sprintf(ret+strlen(ret),"\"dihedral_angles\"->{");
  for(i=0;i<T.ndiags-1;i++) { sprintf(ret + strlen(ret),"%g,",dihedral_angles[i]); }
  sprintf(ret+strlen(ret),"%g}}",dihedral_angles[i]);

  return ret;
}

double   tsmcmc_equilateral_expectation(gsl_rng *rng,double integrand(plCurve *L),
					int max_steps,int max_seconds,
					tsmcmc_triangulation_t T,
					tsmcmc_run_parameters run_params,
					tsmcmc_run_stats *run_stats,
					double *error)
/* 
   This is the "master" driver function for computing an expectation
   over equilateral (unconfined) polygons. It uses the Geyer ips
   estimator to compute error bars.

   We set parameters for the algorithm with run_params (intended to be
   one of the predefined profiles for a run), and return a lot of detailed
   information about the run in run_stats (optional, set to NULL if you 
   don't care). 

   Note that the number of edges is set (implicitly) by the triangulation.
*/


{
  clock_t start,end;
  double  cpu_seconds_used = 0;

  start = clock();
  assert(max_steps > 0); assert(max_seconds > 0);
  assert(run_params.burn_in > 0); assert(run_params.burn_in < max_steps); 

  geyer_ips_t *gips;
  int bufsize = 1000000 < max_steps+10 ? 1000000 : max_steps+10; 
  gips = geyer_ips_new(bufsize);

  if (run_stats != NULL) {

    run_stats->dihedral_steps = 0; 
    run_stats->mp_steps = 0;
    run_stats->permute_steps = 0;

  }
    
  double *edge_lengths, *diagonal_lengths, *dihedral_angles, ival;
  int i;
  plCurve *L;

  tsmcmc_equilateral_ngon(rng,T,&edge_lengths,&diagonal_lengths,&dihedral_angles);
  //tsmcmc_folded_triangle(T,&edge_lengths,&diagonal_lengths,&dihedral_angles);

  for(i=0;i<max_steps && cpu_seconds_used < max_seconds;i++) { 

    tsmcmc_step_t step;
    step = tsmcmc_dihedral_diagonal_permute_step(rng,T,edge_lengths,diagonal_lengths,dihedral_angles,
						 run_params.beta,run_params.delta,run_params.moment_polytope_repeat);

    if (run_stats != NULL) {

      switch (step) {

      case moment_polytope: 
	(run_stats->mp_steps)++; break;
      case dihedral : 
	(run_stats->dihedral_steps)++; break;
      case permute : 
	(run_stats->permute_steps)++; break;
	
      };

    }
      
    if (i > run_params.burn_in) { 

      L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedral_angles);
      if (L != NULL) { 

	ival =integrand(L);
	geyer_ips_add(gips,ival);
	plc_free(L);

      } else {

	printf("Warning! Could not embed a polygon.\n");
	
	printf("Edge lengths were:\n");
	int j;
	for(j=0;j<T.nedges;j++) { printf("%g ",edge_lengths[j]); }
	printf("\n");

	printf("Diagonal lengths were:\n");
	for(j=0;j<T.ndiags;j++) { printf("%17.17g \n",diagonal_lengths[j]); }
	printf("\n");

	printf("Dihedral angles were:\n");
       	for(j=0;j<T.ndiags;j++) { printf("%17.17g \n",dihedral_angles[j]); }
	printf("\n");

	exit(1);
      }

    }

    if (i%1000 == 0) {  /* Check the time every 1,000 steps */

      end = clock();
      cpu_seconds_used = ((double)(end - start))/CLOCKS_PER_SEC;

    }

    if (run_params.logfile != NULL && i%run_params.log_interval == 0 ) {

	char *configuration = tsmcmc_configuration_MathematicaForm(T,edge_lengths,diagonal_lengths,dihedral_angles);
	fprintf(run_params.logfile,"{ \"step number\"->%d, %s, \"integrand\"->%g, \"sample mean\"->%g }\n",
		i,configuration,ival,(double)(gips->sample_sum/(long double)(i)));
	fflush(run_params.logfile);
    }
    
  }
  
  end = clock();
  if (run_stats != NULL) { run_stats->total_seconds = ((double)(end - start))/CLOCKS_PER_SEC; }

  free(edge_lengths); free(diagonal_lengths); free(dihedral_angles);

  double ret,err;
  bool ok;
  
  start = clock();
  ret = geyer_ips_value(gips,&err,&ok);
  end = clock();

  if (error != NULL) { if (ok) { *error = err; } else { *error = 1.0/0.0; } }
  
  if (run_stats != NULL) { 

    run_stats->geyer_ips_seconds = ((double)(end - start))/CLOCKS_PER_SEC; 
    run_stats->max_lagged_covariance_used   = 2*gips->N;
    run_stats->lagged_covariances_available = gips->nldots > 0 ? gips->nldots : gips->buffer_size; 

  }

  geyer_ips_free(&gips);
  
  return ret;

}
  
double   tsmcmc_confined_equilateral_expectation(gsl_rng *rng,double integrand(plCurve *L),
						 double confinement_radius, int nedges,
						 int max_steps,int max_seconds,
						 tsmcmc_run_parameters run_params,
						 tsmcmc_run_stats *run_stats,
						 double *error)

/* This is the "master" driver function for computing an expectation
   over equilateral polygons in "rooted" spherical confinement (the
   confinement is "rooted" when the first vertex of the polygon is at
   the center of the sphere). Since only one triangulation is
   possible, we don't pass in a triangulation.  It uses the Geyer ips
   estimator to compute error bars.

   We set the run parameters with the usual run_params struct, but
   notice that permutation steps aren't possible, so we ignore delta
   (if set). Again, this is usually intended to be one of the
   predefined defaults.

   Detailed information about the run is returned run_stats (optional,
   set to NULL if you don't care).
*/

{
  clock_t start,end;
  double  cpu_seconds_used = 0;

  start = clock();
  assert(max_steps > 0); assert(max_seconds > 0);
  assert(run_params.burn_in > 0); assert(run_params.burn_in < max_steps);

  geyer_ips_t *gips;
  int buffer_size;
  buffer_size = max_steps < 1000000 ? max_steps + 10 : 1000000;

  gips = geyer_ips_new(buffer_size);

  if (run_stats != NULL) {

    run_stats->mp_steps = 0;
    run_stats->dihedral_steps = 0;
    run_stats->permute_steps = 0;

  }

  double *edge_lengths, *diagonal_lengths, *dihedral_angles, ival = 0;
  int i;
  plCurve *L;
  tsmcmc_triangulation_t T;
 
  T = tsmcmc_fan_triangulation(nedges);
  tsmcmc_confined_equilateral_ngon(rng,T,confinement_radius,&edge_lengths,&diagonal_lengths,&dihedral_angles);
  
  for(i=0;i<max_steps && cpu_seconds_used < max_seconds;i++) { 

    tsmcmc_step_t step;

    step = tsmcmc_confined_dihedral_diagonal_step(rng,T,confinement_radius,
						  edge_lengths,diagonal_lengths,dihedral_angles,
						  run_params.beta,run_params.moment_polytope_repeat);

    if (run_stats != NULL) { 
 
      switch (step) {

      case moment_polytope: 
	(run_stats->mp_steps)++; break;
      case dihedral : 
	(run_stats->dihedral_steps)++; break;
      case permute :
	break; // This should never happen in the confined setting.

      };

    }
      
    assert(tsmcmc_edgelengths_equilateral(T,edge_lengths) == true);
    assert(tsmcmc_confined_diagonals_ok(T,confinement_radius,edge_lengths,diagonal_lengths,i) == true);

    if (i > run_params.burn_in) { 

      L = tsmcmc_embed_polygon(T,edge_lengths,diagonal_lengths,dihedral_angles);
      
      if (L != NULL) {

	ival = integrand(L);
	geyer_ips_add(gips,ival);
	plc_free(L);

      } else {

	printf("Warning! Could not embed a polygon.\n");
	
	printf("Edge lengths were:\n");
	int j;
	for(j=0;j<T.nedges;j++) { printf("%g ",edge_lengths[j]); }
	printf("\n");

	printf("Diagonal lengths were:\n");
	for(j=0;j<T.ndiags;j++) { printf("%17.17g ",diagonal_lengths[j]); }
	printf("\n");

	printf("Dihedral angles were:\n");
       	for(j=0;j<T.ndiags;j++) { printf("%17.17g ",dihedral_angles[j]); }
	printf("\n");

	exit(1);

      }

    }

    if (i%1000 == 0) {  /* Check the time every 1,000 steps */

      end = clock();
      cpu_seconds_used = ((double)(end - start))/CLOCKS_PER_SEC;

    }

    if (run_params.logfile != NULL && i%run_params.log_interval == 0 ) {

      char *configuration = tsmcmc_configuration_MathematicaForm(T,edge_lengths,diagonal_lengths,dihedral_angles);
      fprintf(run_params.logfile,"{ \"step number\"->%d, %s, \"integrand\"->%g, \"sample mean\"->%g }\n",
	      i,configuration,ival,(double)(gips->sample_sum/(long double)(i)));
      fflush(run_params.logfile);
      free(configuration);
      
    }
    
  }

  end = clock();

  if (run_stats != NULL) { 
    
    run_stats->total_seconds = ((double)(end - start))/CLOCKS_PER_SEC;

  }

  free(edge_lengths); free(diagonal_lengths); free(dihedral_angles);

  double ret,err;
  bool ok;
  
  start = clock();
  ret = geyer_ips_value(gips,&err,&ok);
  end = clock();
  
  if (run_stats != NULL) { 

    run_stats->geyer_ips_seconds = ((double)(end-start))/CLOCKS_PER_SEC;
    run_stats->max_lagged_covariance_used = 2*gips->N;
    run_stats->lagged_covariances_available = (gips->nldots > 0) ? gips->nldots : gips->buffer_size;

  }

  if (error != NULL) { *error = err; } 
  geyer_ips_free(&gips);
  tsmcmc_triangulation_free(T);

  return ret;

}
