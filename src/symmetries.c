/* 

   symmetries.c : This file is part of the ridgerunner distribution. 
   
   The purpose is to allow us to symmetrize plCurves and variations over
   subgroups of SO(3) as a way to enforce symmetries during RR runs.

*/

#include<plCurve.h>
#include<matrix.h>

#include<config.h>

#ifdef HAVE_MATH_H
  #include<math.h>
#endif

#ifdef HAVE_ASSERT_H
  #include<assert.h>
#endif

#ifdef HAVE_STDLIB_H
  #include<stdlib.h>
#endif

void plc_identity_matrix(plc_matrix *A)
{
  int i,j;
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      (*A)[i][j] = (i == j) ? 1 : 0;
    }
  }
}

void plc_rotation_matrix(plc_vector axis, double theta, plc_matrix *A)

/* Generates the matrix in SO(3) corresponding to rotation around axis with angle theta. */

{
  plc_vector u; /* This will be the normalized axis. */
  u = plc_normalize_vect(axis,NULL);
  
  (*A)[0][0] = cos(theta) + u.c[0]*u.c[0]*(1 - cos(theta));
  (*A)[0][1] = u.c[0]*u.c[1]*(1 - cos(theta)) - u.c[2]*sin(theta);
  (*A)[0][2] = u.c[0]*u.c[2]*(1 - cos(theta)) + u.c[1]*sin(theta);

  (*A)[1][0] = u.c[1]*u.c[0]*(1 - cos(theta)) + u.c[2]*sin(theta);
  (*A)[1][1] = cos(theta) + u.c[1]*u.c[1]*(1 - cos(theta));
  (*A)[1][2] = u.c[1]*u.c[2]*(1 - cos(theta)) - u.c[0]*sin(theta);

  (*A)[2][0] = u.c[2]*u.c[0]*(1 - cos(theta)) - u.c[1]*sin(theta);
  (*A)[2][1] = u.c[2]*u.c[1]*(1 - cos(theta)) + u.c[0]*sin(theta);
  (*A)[2][2] = cos(theta) + u.c[2]*u.c[2]*(1 - cos(theta));
}
  

void plc_reflection_matrix(plc_vector axis, plc_matrix *A)

/* Generates the matrix in O(3) corresponding to reflection across the plane normal to axis. */
/* This turns out to be 

   A[i][j] = delta(i,j) - 2 axis[i]*axis[j]/||axis||^2

*/

{
  plc_vector u;
  int i,j;

  u = plc_normalize_vect(axis,NULL);
  
  for(i=0;i<3;i++) {
   
    for(j=0;j<3;j++) {

      (*A)[i][j] = (i == j ? 1:0) - 2*u.c[i]*u.c[j];

    }
    
  }

}

plc_symmetry *plc_symmetry_new(plCurve *model)

/* Generate an (empty) symmetry which will be applied to the curve in model. 
   We are basically using the model to pick off number-of-components and 
   vertices-per-component data. */

{
  plc_symmetry *sym;
  int cp;

  sym = calloc(1,sizeof(plc_symmetry));
  assert(sym != NULL);
  sym->curve = model;

  sym->target = calloc(model->nc,sizeof(struct plc_vertex_loc *));
  sym->transform = NULL; /* We leave this unallocated, since we expect to get the memory elsewhere. */

  for(cp=0;cp<model->nc;cp++) {

    sym->target[cp] = calloc(model->cp[cp].nv,sizeof(struct plc_vertex_loc));
    assert(sym->target[cp] != NULL);

  }

  return sym;

}
  
void plc_symmetry_free(plc_symmetry **Sym) 

/* Kills a plc_symmetry. Safe for multiple calls on 
   the same pointer, which is set to NULL after the call. */

{
  int cp;
  plc_symmetry *sym;
  sym = *Sym;

  if (sym != NULL) {

    if (sym->transform != NULL) {

      free(sym->transform);
      sym->transform = NULL;

    }

    if (sym->target != NULL) {
      
      for(cp=0;cp<sym->curve->nc;cp++) {
	
	if (sym->target[cp] != NULL) { 

	  free(sym->target[cp]);
	  sym->target[cp] = NULL;

	}
	
      }
    
      free(sym->target);
      sym->target = NULL;
    
    }

    free(sym);
    *Sym = NULL;
  }

}   

plc_symmetry *plc_build_symmetry(plc_matrix *A, plCurve *L)

/* Uses plc_nearest_vertex to try to figure out the "intended" target of each vertex under the symmetry.
   This can fail if the curve is not symmetric to begin with to within < (1/2) an edgelength. In this case, 
   you should build the symmetry map yourself and then symmetrize the results. */

{
  int cp,afterAcp;
  int vt,afterAvt;
  plc_vector afterA;
  bool *used;

  struct plc_nearest_vertex_pc_data *pc_data = NULL;

  /* We are going to mark off "used" vertex targets as they come in. */
  /* If a vertex target is hit twice, we fail. */

  used = calloc(plc_num_verts(L),sizeof(bool));

  /* Now we figure out the (putative) target of each vertex under A */

  plc_symmetry *build;
  build = plc_symmetry_new(L);
  build->transform = plc_matrix_copy(A);  
  
  for(cp=0;cp<L->nc;cp++) {

    for(vt=0;vt<L->cp[cp].nv;vt++) {

      afterA = plc_matrix_vector_multiply(A,L->cp[cp].vt[vt]);
      plc_nearest_vertex(afterA,L,&afterAcp,&afterAvt,&pc_data,NULL);
      
      build->target[cp][vt].cp = afterAcp;
      build->target[cp][vt].vt = afterAvt;

      if (used[plc_vertex_num(L,afterAcp,afterAvt)]) { /* We fail. Cleanup and quit. */
	
	plc_nearest_vertex_pc_data_free(&pc_data);
	plc_symmetry_free(&build);
	free(used);
	return NULL;

      } else { /* Mark this off for the future */

	used[plc_vertex_num(L,afterAcp,afterAvt)] = true;

      }
     
    }
      
  }

  /* We survived this far, so the build must have worked. */

  plc_nearest_vertex_pc_data_free(&pc_data);
  free(used);
  return build;

}

plc_symmetry *plc_symmetry_copy(plc_symmetry *A)
/* Make a new-memory copy of A */
{

  if (A == NULL) { return NULL; } /* Must be safe for no symmetry passed. */
  
  plc_symmetry *build;
  build = plc_symmetry_new(A->curve);
  build->transform = plc_matrix_copy(A->transform);
  
  int cp,vt;
  for(cp=0;cp<A->curve->nc;cp++) {
    for(vt=0;vt<A->curve->cp[cp].nv;vt++) {
      build->target[cp][vt] = A->target[cp][vt];
    }
  }

  return build;

} 

plc_symmetry *plc_compose_symmetries(plc_symmetry *A,plc_symmetry *B)

/* Combine two symmetries in the order A first, then B (or matrix multiplication BA). */

{
  plc_symmetry *build;
  plCurve *L;

  L = A->curve;
  build = plc_symmetry_new(L);

  int cp,vt;

  for(cp=0;cp<L->nc;cp++) {

    for(vt=0;vt<L->cp[cp].nv;vt++) {

      build->target[cp][vt] = B->target[A->target[cp][vt].cp][A->target[cp][vt].vt];

    }
    
  }

  build->transform = plc_matrix_matrix_multiply(B->transform,A->transform);

  return build;

}

/************************** Symmetry GROUP code **************************/

plc_symmetry_group *plc_symmetry_group_new(int n)
{

  if (n <= 0) { return NULL; } /* A little bit of sanity checking */

  plc_symmetry_group *build;
  build = calloc(1,sizeof(plc_symmetry_group));

  build->n = n;
  build->sym = calloc(n,sizeof(plc_symmetry *));
  build->inverse = calloc(n,sizeof(int));

  return build;

}

void plc_symmetry_group_free(plc_symmetry_group **G) 
{
  int i;

  if (G != NULL) {

    if ((*G) != NULL) {

      for(i=0;i<(*G)->n;i++) {   /* Free each actual symmetry */

	plc_symmetry_free(&((*G)->sym[i]));
	
      }

      free((*G)->sym);         
      (*G)->sym = NULL;      /* Free the buffers of group information */
      
      free((*G)->inverse);     
      (*G)->inverse = NULL;

      free(*G);
      *G = NULL;
    
    }

  }

}

plc_symmetry_group *plc_symmetry_group_copy(plc_symmetry_group *G) 
{

  if (G == NULL) {  /* The procedure must be safe for no group passed. */

    return NULL;

  }

  plc_symmetry_group *build;
  build = plc_symmetry_group_new(G->n);
  
  int i;
  for(i=0;i<G->n;i++) {
    build->sym[i] = plc_symmetry_copy(G->sym[i]);
  }

  for(i=0;i<G->n;i++) {
    build->inverse[i] = G->inverse[i];
  }

  return build;
}


plc_symmetry_group *plc_rotation_group(plCurve *L,plc_vector axis, int n)

/* Creates the symmetry group Z/nZ of rotations around axis. */

{
  double TWO_PI = 6.2831853071795864769;
  plc_symmetry_group *build;

  build = plc_symmetry_group_new(n);

  /* We now try to set up the group. */

  plc_matrix A;
  plc_identity_matrix(&A);
  build->sym[0] = plc_build_symmetry(&A,L);

  int i;
  for(i=1;i<n;i++) {
    plc_rotation_matrix(axis,i*(TWO_PI/(double)(n)),&A);
    build->sym[i] = plc_build_symmetry(&A,L);
  }

  build->inverse[0] = 0;
  for(i=1;i<n;i++) {
    build->inverse[i] = n - i;  /* The group law in Z/nZ: i + (n-i) = n = 0 (mod n) */
  }

  /* We need to scan and make sure the builds succeeded. */

  for(i=0;i<n;i++) {

    if (build->sym[i] == NULL) {

      plc_symmetry_group_free(&build);
      return NULL;

    }

  }

  return build;

}

plc_symmetry_group *plc_reflection_group(plCurve *L,plc_vector axis)

/* Creates the Z/2Z group of reflections over the plane normal to axis. */

{
  plc_symmetry_group *build;
  build = plc_symmetry_group_new(2);

  /* We now try to set up the group. */

  plc_matrix A;
  plc_identity_matrix(&A);
  build->sym[0] = plc_build_symmetry(&A,L);

  plc_reflection_matrix(axis,&A);
  build->sym[1] = plc_build_symmetry(&A,L);

  build->inverse[0] = 0; build->inverse[1] = 1; /* Both elements are their own inverses */

  /* We need to scan and make sure the builds succeeded. */

  int i;
  for(i=0;i<2;i++) {

    if (build->sym[i] == NULL) {

      plc_symmetry_group_free(&build);
      return NULL;

    }

  }

  return build;

}

plc_symmetry_group *plc_coordplanes_reflection_group(plCurve *L)

/* Creates reflections over the x-y, y-z, and z-x planes. */

{
  plc_symmetry_group *build;
  build = plc_symmetry_group_new(4);

  /* We now try to set up the group. */

  plc_matrix A;
  plc_identity_matrix(&A);
  build->sym[0] = plc_build_symmetry(&A,L);

  plc_reflection_matrix(plc_build_vect(1,0,0),&A);
  build->sym[1] = plc_build_symmetry(&A,L);

  plc_reflection_matrix(plc_build_vect(0,1,0),&A);
  build->sym[2] = plc_build_symmetry(&A,L);
  
  plc_reflection_matrix(plc_build_vect(0,0,1),&A);
  build->sym[3] = plc_build_symmetry(&A,L);

  /* All elements are their own inverses */

  build->inverse[0] = 0; build->inverse[1] = 1; 
  build->inverse[2] = 2; build->inverse[3] = 3;

  /* We need to scan and make sure the builds succeeded. */

  int i;
  for(i=0;i<2;i++) {

    if (build->sym[i] == NULL) {

      plc_symmetry_group_free(&build);
      return NULL;

    }

  }

  return build;

}

/* We are at last ready to symmetrize a curve over a group! The algorithm */
/* is simple: for each point, we collect its G-orbit, average those points, */
/* and update the point accordingly. Then we push the average around the orbit. */

/* It takes a little work to make sure that we don't hit points multiple times. */
/* Since the symmetry group is stored in L->G, we don't pass it separately */

void plc_symmetrize(plCurve *L)

{
  bool *point_covered;
  plc_vector symmetrized_pos, orbit_vec;
  int cp,vt;
  struct plc_vertex_loc target;
  int i;

  if (L->G == NULL) {  /* If the plCurve has no associated symmetry group, we do nothing */

    return;

  }

  point_covered = calloc(plc_num_verts(L),sizeof(bool));
  
  for(cp=0;cp<L->nc;cp++) {
    
    for(vt=0;vt<L->cp[cp].nv;vt++) {
      
      if (!point_covered[plc_vertex_num(L,cp,vt)]) {  /* We have not covered this orbit yet. */
	
	for(i=0,symmetrized_pos = plc_build_vect(0,0,0);i<L->G->n;i++) {  /* Collect the average vector over orbit at cur pt */
	  
	  target = L->G->sym[i]->target[cp][vt];
	  orbit_vec = plc_matrix_vector_multiply(L->G->sym[L->G->inverse[i]]->transform,
						 L->cp[target.cp].vt[target.vt]);

	  plc_M_vmadd(symmetrized_pos,1/(double)(L->G->n),orbit_vec);

	}

	/* Now push that vector forward around the orbit, marking the vertices as covered */

	for(i=0;i<L->G->n;i++) {

	  target = L->G->sym[i]->target[cp][vt];
	  L->cp[target.cp].vt[target.vt] = plc_matrix_vector_multiply(L->G->sym[i]->transform,symmetrized_pos);
	  point_covered[plc_vertex_num(L,target.cp,target.vt)] = true;

	}

      }

    }

  }
	  
  free(point_covered);

}

void plc_symmetrize_variation(plCurve *L,plc_vector *buffer)

{
  bool *point_covered;
  plc_vector symmetrized_pos, orbit_vec;
  int cp,vt;
  struct plc_vertex_loc target;
  int i;

  if (L->G == NULL) { return; } /* If L has no associated symmetry group, there is nothing to do. */

  point_covered = calloc(plc_num_verts(L),sizeof(bool));

  for(cp=0;cp<L->nc;cp++) {

    for(vt=0;vt<L->cp[cp].nv;vt++) {

      if (!point_covered[plc_vertex_num(L,cp,vt)]) {  /* We have not covered this orbit yet. */
	
	for(i=0,symmetrized_pos = plc_build_vect(0,0,0);i<L->G->n;i++) {  /* Collect the average vector over orbit at cur pt */

	  target = L->G->sym[i]->target[cp][vt];
	  orbit_vec = plc_matrix_vector_multiply(L->G->sym[L->G->inverse[i]]->transform,
						 buffer[plc_vertex_num(L,target.cp,target.vt)]);

	  plc_M_vmadd(symmetrized_pos,1/(double)(L->G->n),orbit_vec);

	}

	/* Now push that vector forward around the orbit, marking the vertices as covered */

	for(i=0;i<L->G->n;i++) {

	  target = L->G->sym[i]->target[cp][vt];
	  buffer[plc_vertex_num(L,target.cp,target.vt)] = plc_matrix_vector_multiply(L->G->sym[i]->transform,symmetrized_pos);
	  point_covered[plc_vertex_num(L,target.cp,target.vt)] = true;

	}

      }

    }

  }
	  
  free(point_covered);

}


double plc_symmetry_check(plCurve *L,plc_symmetry *A)

{
  int cp,vt;
  double maxdist = 0.0,thisdist;
  plc_vector afterA;
  
  if (L == NULL || A == NULL) { return 0; }  /* No curve or no symmetry => no error */

  for(cp=0;cp<L->nc;cp++) {

    for(vt=0;vt<L->cp[cp].nv;vt++) {

      afterA = plc_matrix_vector_multiply(A->transform,L->cp[cp].vt[vt]);
      thisdist =  plc_distance(afterA,L->cp[A->target[cp][vt].cp].vt[A->target[cp][vt].vt]);
      maxdist = (maxdist > thisdist) ? maxdist : thisdist;

    }

  }

  return maxdist;
}

double plc_symmetry_variation_check(plCurve *L,plc_vector *buffer,plc_symmetry *A)

{
  int cp,vt;
  double maxdist = 0.0,thisdist;
  plc_vector afterA;
  
  if (L == NULL || A == NULL) { return 0; }  /* No curve or no symmetry => no error */

  for(cp=0;cp<L->nc;cp++) {

    for(vt=0;vt<L->cp[cp].nv;vt++) {

      afterA = plc_matrix_vector_multiply(A->transform,buffer[plc_vertex_num(L,cp,vt)]);
      thisdist =  plc_distance(afterA,buffer[plc_vertex_num(L,A->target[cp][vt].cp,A->target[cp][vt].vt)]);
      maxdist = (maxdist > thisdist) ? maxdist : thisdist;

    }

  }

  return maxdist;
}

  /* To check an entire group, use plc_symmetry_group_check, which checks the entire 
     group L->G and returns the maximum error. */

double plc_symmetry_group_check(plCurve *L)
{

  int i;
  double maxdist = 0,thisdist;

  if (L->G == NULL) { return 0; } /* If there is no symmetry group, there is no error. */

  for(i=0;i<L->G->n;i++) {

    thisdist = plc_symmetry_check(L,L->G->sym[i]);
    maxdist = (maxdist > thisdist) ? maxdist : thisdist;

  } 

  return maxdist;

}

double plc_symmetry_group_variation_check(plCurve *L,plc_vector *buffer)
{
  int i;
  double maxdist = 0,thisdist;

  if (L->G == NULL) { return 0; } /* If there is no symmetry group, there is no error. */

  for(i=0;i<L->G->n;i++) {

    thisdist = plc_symmetry_variation_check(L,buffer,L->G->sym[i]);
    maxdist = (maxdist > thisdist) ? maxdist : thisdist;

  } 

  return maxdist;

}
