/*
 *  Routines to generate random self-avoiding polygons as part of plCurve.
 *
 */

/* Copyright 2023 The University of Georgia. */

/* This file is part of plCurve.

plCurve is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

plCurve is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with plCurve; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#include"plcRandomSap.h"

uint64_t *plc_xoshiro_init(uint64_t init) {

  uint64_t *seed = malloc(4*sizeof(uint64_t));
  
  uint64_t z;
  int i;
  
  for(i=0;i<4;i++ ) {

    z = (init += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
  
    seed[i] = z ^ (z >> 31);

  }

  return seed;
  
}

void plc_xoshiro_free(uint64_t *rng) {

  free(rng);

}

inline uint64_t rotl(const uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}

uint64_t plc_xoshiro_next(uint64_t *s) {

  const uint64_t result = s[0] + s[3];
  
  const uint64_t t = s[1] << 17;
  
  s[2] ^= s[0];
  s[3] ^= s[1];
  s[1] ^= s[2];
  s[0] ^= s[3];
  
  s[2] ^= t;
  
  s[3] = rotl(s[3], 45);
  
  return result;
}

inline double plc_xoshiro_uniform(uint64_t *state) {
  const uint64_t r = plc_xoshiro_next(state);
  const double _nrm=1.0/(1ull<<52);
  return (r>>12)*_nrm;
}
 
/* This is the jump function for the generator. It is equivalent 
    to 2^128 calls to next(); it can be used to generate 2^128 
    non-overlapping subsequences for parallel computations. */ 

void plc_xoshiro_jump(uint64_t *s) {
  static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };
  
  uint64_t s0 = 0;
  uint64_t s1 = 0;
  uint64_t s2 = 0;
  uint64_t s3 = 0;

  for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
    for(int b = 0; b < 64; b++) {
      if (JUMP[i] & UINT64_C(1) << b) {
	s0 ^= s[0];
	s1 ^= s[1];
	s2 ^= s[2];
	s3 ^= s[3];
      }
      plc_xoshiro_next(s);
    }
  
  s[0] = s0;
  s[1] = s1;
  s[2] = s2;
  s[3] = s3;
}


/* /\* This is the long-jump function for the generator. It is equivalent to */
/*    2^192 calls to next(); it can be used to generate 2^64 starting points, */
/*    from each of which jump() will generate 2^64 non-overlapping */
/*    subsequences for parallel distributed computations. *\/ */

/* void long_jump(void) { */
/*   static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635 }; */
  
/*   uint64_t s0 = 0; */
/*   uint64_t s1 = 0; */
/*   uint64_t s2 = 0; */
/*   uint64_t s3 = 0; */
/*   for(int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++) */
/*     for(int b = 0; b < 64; b++) { */
/*       if (LONG_JUMP[i] & UINT64_C(1) << b) { */
/* 	s0 ^= s[0]; */
/* 	s1 ^= s[1]; */
/* 	s2 ^= s[2]; */
/* 	s3 ^= s[3]; */
/*       } */
/*       next();	 */
/*     } */
  
/*   s[0] = s0; */
/*   s[1] = s1; */
/*   s[2] = s2; */
/*   s[3] = s3; */
/* } */


bool plc_is_sap_internal(plCurve *L, bool verbose) {

  /* This is an internal debugging function which double-checks that all
     vertices of L are at least (almost) distance 1.0 away from each other.
     It's written for the general case where L has many components, even
     though we don't have code to generate such SAPs yet. */

  int cp1,cp2;
  int vt1,vt2;

  for(cp1=0;cp1<L->nc;cp1++) {

    for(cp2=0;cp2<=cp1;cp2++) {

      if (cp1 == cp2) {
	
	for(vt1=1;vt1<L->cp[cp1].nv;vt1++) {
	  
	  for(vt2=0;vt2<vt1;vt2++) {
	    
	    if (plc_sq_dist(L->cp[cp1].vt[vt1],L->cp[cp2].vt[vt2]) < 1.0-1e-10) {

	      if (verbose) {

		printf("plc_is_sap: Squared distance between L->cp[%d].vt[%d] and L->cp[%d].vt[%d] is %g < %g\n",
		       cp1,vt1,cp2,vt2,plc_sq_dist(L->cp[cp1].vt[vt1],L->cp[cp2].vt[vt2]),1.0-1e-10);
		       
		return false;

	      }
	      
	    }
	    
	  }
	  
	}
	
      } else {
	
	for(vt1=0;vt1<L->cp[cp1].nv;vt1++) {
	  
	  for(vt2=0;vt2<L->cp[cp2].nv;vt2++) {
	    
	     if (verbose) {

		printf("plc_is_sap: Squared distance between L->cp[%d].vt[%d] and L->cp[%d].vt[%d] is %g < %g\n",
		       cp1,vt1,cp2,vt2,plc_sq_dist(L->cp[cp1].vt[vt1],L->cp[cp2].vt[vt2]),1.0-1e-10);
		       
		return false;

	     }

	  }

	}

      }

    }

  }

  return true;

}

bool plc_is_sap(plCurve *L) {

  return plc_is_sap_internal(L,false);

}


plCurve *plc_random_equilateral_closed_self_avoiding_polygon(uint64_t *xos,int n)
/* 
   Generates random closed polygons where vertices are surrounded by a 
   disjoint spheres of radius 1/2 (the "string of pearls" model) by 
   rejection sampling. This is exponentially slow, so it's limited both 
   in number of attempts and by the max and min n which will be attempted. 

   We also do a fair amount of one-off optimization here in order to increase 
   performance, including the use of an internal random number generator 
   instead of our usual calls to gsl. 

   To call this function, we need to use 

   uint64_t *xos =  plc_xoshiro_init((uint64_t)(time(0)));

   and then pass the resulting xos pointer to this function with each call, 
   freeing the xoshiro state pointer xos with 

   plc_xoshiro_free(xos);

   after runs are complete.

*/
{
  int nv = n,cc = {0};
  bool open = {false};
  plCurve *L = plc_new(1,&nv,&open,&cc);

  /* We'll copy over our sap to the buffer in L once we generate it. */

  double *X[3];  /* This will be our working buffer of vertex positions. */
 
  assert(n >= 5);
  if (n < 5) {

    printf("Error. Can't call plc_random_sap\n"
	   "with n < 5\n");
    exit(1);

  }

  assert(n <= 25);
  if (n > 25) {

    printf("Error. Can't call plc_random_sap\n"
	   "with n > 25\n");
    exit(1);
  }

  /* The design of this code is basically straight from the Mathematica notebook spaam-sampling.nb */
  /* These are buffers of vertex coordinates. In order to address the individual vectors, we use */
  /* the **seemingly backwards** notation X[k][i] to get the k-th entry of the i-th vector. */

  X[0] = malloc(n * sizeof(double));
  X[1] = malloc(n * sizeof(double));
  X[2] = malloc(n * sizeof(double));

  assert(X[0] != NULL); assert(X[1] != NULL); assert(X[2] != NULL);

  X[0][0] = 0.; X[1][0] = 0.; X[2][0] = 0.;
  X[0][1] = 1.; X[1][1] = 0.; X[2][1] = 0.;

  double frame[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};

  /* To match the memory usage in the X buffers, we'll address these "backwards" too: frame[k][i] is */
  /* the k-th coordinate of the i-th vector. */
  
  double last_diag, this_diag = 0.;
  double normal[3];
  int i,j;
  int trials=0,max_trials = 50000000;

  double TWOPI = 6.2831853071795864769;
  double cos_phi, sin_phi;
  double theta, cos_theta, sin_theta;
  double nm;

  restart_trial:

  trials++;
  assert(trials < max_trials);
  
  normal[0] = 0.; normal[1] = 0.; normal[2] = 1.;
  last_diag = 1.0;

  for(i=2;i<(n-1); i++, last_diag = this_diag) {

    this_diag = last_diag + 2.0*(plc_xoshiro_uniform(xos)-0.5);
    if ((this_diag + last_diag < 1.0) || (this_diag <= 0.0)) { goto restart_trial; };
    
    /* frame[-][0] = X[-][i-1]/||X[-][i-1]|| */
    
    nm = sqrt(X[0][i-1]*X[0][i-1] + X[1][i-1]*X[1][i-1] + X[2][i-1]*X[2][i-1]);
    frame[0][0] = X[0][i-1]/nm; frame[1][0] = X[1][i-1]/nm; frame[2][0] = X[2][i-1]/nm;
    
    /* frame[-][1] = normal x frame[-][0] */
      
    frame[0][1] = normal[1]*frame[2][0] - normal[2]*frame[1][0];  
    frame[1][1] = normal[2]*frame[0][0] - normal[0]*frame[2][0];
    frame[2][1] = normal[0]*frame[1][0] - normal[1]*frame[0][0];
    
    /* frame[-][1] *= 1/||f[-][1]|| */
    
    nm = sqrt(frame[0][1]*frame[0][1] + frame[1][1]*frame[1][1] + frame[2][1]*frame[2][1]);
    frame[0][1] /= nm;
    frame[1][1] /= nm;
    frame[2][1] /= nm;
    
    cos_phi = (this_diag * this_diag + last_diag * last_diag - 1.0)/(2.0 * this_diag * last_diag);
    sin_phi = sqrt(1 - cos_phi * cos_phi);
      
    /* X[-][i] = this_diag (cos_phi * frame[-][0] + sin_phi * frame[-][1]); */
    
    X[0][i] = this_diag * (cos_phi * frame[0][0] + sin_phi * frame[0][1]);
    X[1][i] = this_diag * (cos_phi * frame[1][0] + sin_phi * frame[1][1]);
    X[2][i] = this_diag * (cos_phi * frame[2][0] + sin_phi * frame[2][1]);
    
    /* frame[-][2] = -sin_phi frame[-][0] + cos_phi frame[-][1] */
    
    frame[0][2] = -sin_phi * frame[0][0] + cos_phi * frame[0][1];
    frame[1][2] = -sin_phi * frame[1][0] + cos_phi * frame[1][1];
    frame[2][2] = -sin_phi * frame[2][0] + cos_phi * frame[2][1];
    
    /* Now pick a dihedral angle theta to try. */
    
    theta = TWOPI*plc_xoshiro_uniform(xos);
    cos_theta = cos(theta);
    sin_theta = sin(theta);
      
    /* normal = cos(theta) normal + sin(theta) frame[-][2]; */
      
    normal[0] = cos_theta * normal[0] + sin_theta * frame[0][2];
    normal[1] = cos_theta * normal[1] + sin_theta * frame[1][2];
    normal[2] = cos_theta * normal[2] + sin_theta * frame[2][2];
    
    /* Now check distances. We go backwards because the new sphere is most likely */
    /* to collide with nearby spheres. */
    
    for(j=i-1;j>=0;j--) {

      if (((X[0][i]-X[0][j])*(X[0][i]-X[0][j]) +
	   (X[1][i]-X[1][j])*(X[1][i]-X[1][j]) +
	   (X[2][i]-X[2][j])*(X[2][i]-X[2][j])) < 1.0-(1e-9)) { goto restart_trial; }
	
    }

  }

  /* If we got to this point, we're at i = n - 1; the last vertex! */

  assert(i==n-1); 
  last_diag = this_diag;
  
  if (last_diag <= 0.0 || last_diag >= 2.0) { goto restart_trial; }
    
  /* frame[-][0] = X[-][i-1]/||X[-][i-1]|| */
  
  nm = sqrt(X[0][i-1]*X[0][i-1] + X[1][i-1]*X[1][i-1] + X[2][i-1]*X[2][i-1]);
  frame[0][0] = X[0][i-1]/nm; frame[1][0] = X[1][i-1]/nm; frame[2][0] = X[2][i-1]/nm;
  
  /* frame[-][1] = normal x frame[-][0] */
  
  frame[0][1] = normal[1]*frame[2][0] - normal[2]*frame[1][0];  
  frame[1][1] = normal[2]*frame[0][0] - normal[0]*frame[2][0];
  frame[2][1] = normal[0]*frame[1][0] - normal[1]*frame[0][0];
  
  /* frame[-][1] *= 1/||f[-][1]|| */
  
  nm = sqrt(frame[0][1]*frame[0][1] + frame[1][1]*frame[1][1] + frame[2][1]*frame[2][1]);
  frame[0][1] /= nm;
  frame[1][1] /= nm;
  frame[2][1] /= nm;
  
  /* Since we know that this_diag is 1.0, this formula simplifies... */
	
  cos_phi = last_diag/2.0;
  sin_phi = sqrt(1.0 - cos_phi * cos_phi);
  
  /* X[-][i] = this_diag (cos_phi * frame[-][0] + sin_phi * frame[-][1]); */
	
  X[0][i] = cos_phi * frame[0][0] + sin_phi * frame[0][1];
  X[1][i] = cos_phi * frame[1][0] + sin_phi * frame[1][1];
  X[2][i] = cos_phi * frame[2][0] + sin_phi * frame[2][1];
  
  /* This is the last vertex, so there's no need to update the frame, but we still */
  /* check distances. We go backwards because the new sphere is most likely */
  /* to collide with nearby spheres. */
	
  for(j=i-1;j>=0;j--) {
	  
    if (((X[0][i]-X[0][j])*(X[0][i]-X[0][j]) +
	 (X[1][i]-X[1][j])*(X[1][i]-X[1][j]) +
	 (X[2][i]-X[2][j])*(X[2][i]-X[2][j])) < 1.0-(1e-9)) { goto restart_trial; };
	
  }
	
  /* If we make it to here, it's time to copy the contents of the buffers into the plCurve data structure */

  for(i=0;i<n;i++) {
    
    L->cp[0].vt[i].c[0] = X[0][i];
    L->cp[0].vt[i].c[1] = X[1][i];
    L->cp[0].vt[i].c[2] = X[2][i];
    
  }

  plc_fix_wrap(L);

  /* Now we need to clean up the temporary buffers */

  free(X[0]); free(X[1]); free(X[2]);

  return L;

} 
