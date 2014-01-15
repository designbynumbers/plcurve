/* 

   Code to manage the triangulation data structure. 

*/

#include "tsmcmc.h"
#include <assert.h>
#include <math.h>
#include <string.h>



tsmcmc_triangulation_t tsmcmc_fan_triangulation(int nedges)
/* Constructs the fan triangulation for an n-edge polygon. */
{

  tsmcmc_triangulation_t T;
  tsmcmc_chord_t *chordsys = calloc(nedges-3,sizeof(tsmcmc_chord_t));
  assert(chordsys != NULL);

  int i;
  for(i=0;i<nedges-3;i++) { 

    chordsys[i].vt[0] = 0; chordsys[i].vt[1] = i+2;

  }
  assert(tsmcmc_chord_system_ok(nedges-3,chordsys));

  T = tsmcmc_generate_triangulation(nedges,chordsys);
  free(chordsys);

  return T;
}

tsmcmc_triangulation_t tsmcmc_spiral_triangulation(int nedges)
/* Constructs a triangulation by recursively joining adjacent edges with a diagonal. */
{
  bool *vert_active = calloc(nedges,sizeof(bool));
  assert(vert_active != NULL);

  tsmcmc_chord_t *chordsys = calloc(nedges-3,sizeof(tsmcmc_chord_t));
  assert(chordsys != NULL);

  /* Start with all vertices active. */

  int i,cd;
  for(i=0;i<nedges;i++) { vert_active[i] = true; }

  for(cd=0,chordsys[cd].vt[0] = 0,i=0;cd < nedges-3;cd++) { 

    /* Start assuming that chordsys[cd].vt[0] is set, and equal to i. */
    /* Scan forward until the next active vertex */
    for(i = (i+1)%nedges;!vert_active[i];i = (i+1)%nedges);

    /* Now mark that vertex off; we're skipping it. */
    vert_active[i] = false;

    /* Scan forward again, looking for the next active vertex. */
    for(i = (i+1)%nedges;!vert_active[i];i = (i+1)%nedges);

    /* This is the other end of this chord, so book it as the end
       of the chord (and if there's another chord, the start of 
       the next one). */

    chordsys[cd].vt[1] = i;
    if(cd < nedges-4) { chordsys[cd+1].vt[0] = i; }

    /* We'll keep the loop up until we triangulate the polygon. */
  
  }
  
  assert(tsmcmc_chord_system_ok(nedges-3,chordsys));

  tsmcmc_triangulation_t T = tsmcmc_generate_triangulation(nedges,chordsys);
  assert(tsmcmc_triangulation_ok(T));

  free(chordsys);
  free(vert_active);

  return T;
}

void tsmcmc_random_triangulation_worker(gsl_rng *rng,int *chords_used,tsmcmc_chord_t *chordsys,int *polygon,int polsize)

/* Split polygon along a random diagonal, recurse on the two sub-polygons. */

{
  if (polsize == 3) { 

    return; 

  } /* A single triangle has no diagonals */

  /* Choose a diagonal randomly on a polygon of size polsize.  We will
     generate the diagonal in sorted order. Note that if we pick the
     0th vertex, we have to select the next from 2, ..., polsize - 2
     (inclusive), but if we pick any other vert for chord[0], we get
     to go from that + 2 to polsize - 1 (inclusive). */
  
  int chord[2];
  chord[0] = gsl_rng_uniform_int(rng,polsize-2); 
  chord[1] = chord[0] + 2 + gsl_rng_uniform_int(rng,polsize - chord[0] - (chord[0] == 0 ? 3 : 2));

  /* Now translate this into a new chord for the chord system */

  chordsys[*chords_used].vt[0] = polygon[chord[0]];
  chordsys[*chords_used].vt[1] = polygon[chord[1]];
  (*chords_used)++;

  /* Now come up with new daughter polygons. */

  int daughterlenA,daughterlenB; 
  daughterlenB = polsize - (chord[1] - chord[0] - 1); /* Polsize - #verts between chord[0] and chord[1] */
  daughterlenA = polsize + 2 - daughterlenB;          /* Everybody else, remembering that chord[0] and chord[1] 
							 appear in both daughters */
  int *daughterA, *daughterB; 
  daughterA = calloc(daughterlenA,sizeof(int)); assert(daughterA != NULL);
  daughterB = calloc(daughterlenB,sizeof(int)); assert(daughterB != NULL);

  int i;
  for(i=0;i<daughterlenA;i++) { daughterA[i] = polygon[chord[0]+i]; }
  for(i=0;i<daughterlenB;i++) { daughterB[i] = polygon[(chord[1]+i)%polsize]; }

  tsmcmc_random_triangulation_worker(rng,chords_used,chordsys,daughterA,daughterlenA);
  tsmcmc_random_triangulation_worker(rng,chords_used,chordsys,daughterB,daughterlenB);

  free(daughterA); free(daughterB); 

}

tsmcmc_triangulation_t tsmcmc_random_triangulation(gsl_rng *rng,int nedges) 

/* Constructs a triangulation by randomly filling in chords. The
   algorithm is fairly simple: we keep splitting the polygon in two by
   choosing a chord randomly, and recursing on each of the sub-polygons. */

{
  int *start_polygon = calloc(nedges,sizeof(int));
  assert(start_polygon != NULL);
  int i;
  for(i=0;i<nedges;i++) { start_polygon[i] = i; }

  tsmcmc_chord_t *chordsys = calloc(nedges-3,sizeof(tsmcmc_chord_t *));
  assert(chordsys != NULL);
  int chords_used = 0;
  
  tsmcmc_random_triangulation_worker(rng,&chords_used,chordsys,start_polygon,nedges);
  assert(tsmcmc_chord_system_ok(nedges-3,chordsys));
  
  tsmcmc_triangulation_t T = tsmcmc_generate_triangulation(nedges,chordsys);
  assert(tsmcmc_triangulation_ok(T));

  free(chordsys);
  free(start_polygon);

  return T;
}

tsmcmc_triangulation_t tsmcmc_teeth_triangulation(int nedges)

/* Constructs a triangulation by "walking" across the polygon. */

{
  tsmcmc_chord_t *chordsys = calloc(nedges-3,sizeof(tsmcmc_chord_t *));
  assert(chordsys != NULL);
  
  int chords_used;
  chordsys[0].vt[0] = nedges-1; chordsys[0].vt[1] = 1; /* We start by skipping vertex 0. */
  chords_used = 1;

  for(;chords_used < nedges-3;chords_used++) { 

    chordsys[chords_used] = chordsys[chords_used-1];

    if (chords_used % 2 != 0) {  /* On ODD chords, move vertex 1 forward */
      
      chordsys[chords_used].vt[1]++; 

    } else { /* On EVEN chords, move vertex 0 backwards */

      chordsys[chords_used].vt[0]--; 

    }

  }
  
  tsmcmc_triangulation_t T = tsmcmc_generate_triangulation(nedges,chordsys);
  assert(tsmcmc_triangulation_ok(T));

  free(chordsys);
  
  return T;
}
 

int tsmcmc_vert_chord(tsmcmc_triangulation_t T,int vt,int chord,chordtype_t type) 
{
  assert(vt == 0 || vt == 1);
  assert(type == diagonal || type == edge);

  if (type == diagonal) {

    assert(chord >= 0 && chord <= T.ndiags);
    return T.diags[chord].vt[vt];

  } else { 
    
    assert(chord >= 0 && chord <= T.nedges);
    return T.edges[chord].vt[vt];

  }
}

double tsmcmc_chord_length(int chord,chordtype_t type,
			   double *diagonal_lengths,double *edge_lengths) 
{

  if (type == diagonal) { 

    return diagonal_lengths[chord];

  } else if (type == edge) {

    return edge_lengths[chord];

  } else {

    printf("tsmcmc_chord_length: illegal chord type %d is not diagonal = %d or edge = %d\n",
	   type,diagonal,edge);
    exit(1);

  }

}

  
bool tsmcmc_vert_on_chord(tsmcmc_triangulation_t T,int vt,int chord,chordtype_t type)
/* Checks to see whether the vertex vt is on the given chord. */
{
  if (type == diagonal) { 

    return (T.diags[chord].vt[0] == vt || T.diags[chord].vt[1] == vt);

  } else if (type == edge) {

    return (T.edges[chord].vt[0] == vt || T.edges[chord].vt[1] == vt);

  } else {

    printf("illegal chord type %d is not diagonal = %d or edge = %d \n",
	   type,diagonal,edge);
    exit(1);

  }

}

bool tsmcmc_chord_ok(tsmcmc_triangulation_t T,int chord,chordtype_t type) 
{
/* Self-checks, mostly for indices in bounds for a chord in a triangulation. */
  if (type != diagonal && type != edge) {

    printf("bad chord: chord type %d is not diagonal = %d or edge = %d\n",
	   type,diagonal,edge);
    return false;

  } 

  if (type == edge) {

    if (chord < 0 || chord > T.nedges-1) { 

      printf("bad chord: type edge, index %d is < 0 or > T.nedges-1 = %d.\n",
	     chord,T.nedges-1);
      return false;
      
    }

    int i;

    for(i=0;i<2;i++) {

      if (T.edges[chord].vt[i] < 0 || T.edges[chord].vt[i] > T.nedges-1) {
	
	printf("bad chord: type edge, vertex %d = %d is < 0 or > T.nedges - 1 = %d\n",
	       i,T.edges[chord].vt[i],T.nedges-1);
	return false;

      }

    }

  } else {

    if (chord < 0 || chord > T.ndiags-1) { 

      printf("bad chord: type diagonal, index %d is < 0 or > T.ndiags-1 = %d.\n",
	     chord,T.ndiags-1);
      return false;
      
    }

    int i;

    for(i=0;i<2;i++) {

      if (T.diags[chord].vt[i] < 0 || T.diags[chord].vt[i] > T.nedges-1) {
	
	printf("bad chord: type diagonal, vertex %d = %d is < 0 or > T.nedges - 1 = %d\n",
	       i,T.diags[chord].vt[i],T.nedges-1);
	return false;

      }
      
    }

  }

  return true;

}
	   

void tsmcmc_print_chord(tsmcmc_triangulation_t T,int chord,chordtype_t type)
{
  if (type == edge) {

    printf("edge %d (%d <-> %d)",chord,T.edges[chord].vt[0],T.edges[chord].vt[1]);

  } else {

    printf("diagonal %d (%d <-> %d)",chord,T.diags[chord].vt[0],T.diags[chord].vt[1]);

  }
}

bool tsmcmc_triangulation_ok(tsmcmc_triangulation_t T)
/* Self-tests on a triangulation. */
{
  /* First, we check that the numbers of triangles, diagonals, and verts
     are compatible with the claimed number of edges. */

  if (T.ntri != T.nedges-2) {
    
    printf("T.ntri = %d != T.nedges-2 = %d.\n",T.ntri,T.nedges-2);
    return false;

  }

  if (T.ndiags != T.nedges-3) { 

    printf("T.diags = %d != T.nedges-3 = %d.\n",T.ndiags,T.nedges-3);
    return false;

  }

  /* Now check for unallocated buffers */

  if (T.triangles == NULL) {

    printf("T.triangles == NULL\n");
    return false;

  } 

  if (T.diags == NULL) {

    printf("T.diags == NULL\n");
    return false;

  }

  if (T.edges == NULL) {

    printf("T.edges == NULL\n");
    return false;

  }

  /* Now check the vertex indices in the edge chords */

  int i,j;

  for(i=0;i<T.nedges;i++) {

    if (!tsmcmc_chord_ok(T,i,edge)) { return false; }

  }
    
  /* Now check vertex indices in the "diagonal" chords. */

  for(i=0;i<T.ndiags;i++) {

    if (!tsmcmc_chord_ok(T,i,diagonal)) { return false; }

  }

  /* Now check the triangles, mostly for legal indices. We will also check that all but */
  /* one of the triangles appears as a daughter exactly once in the list of daughters. */

  int *daughterof;
  daughterof = calloc(T.ntri,sizeof(int));
  for(i=0;i<T.ntri;i++) { daughterof[i] = -1; }

  for(i=0;i<T.ntri;i++) {

    if (!tsmcmc_chord_ok(T,T.triangles[i].parent_chord,T.triangles[i].parent_type)) { 

      printf("T.triangles[%d].parent chord not ok.\n",i);
      return false;

    }
      
    for(j=0;j<2;j++) { /* Check each of the daughter chords and daughter triangles for legal indices. */

      if (!tsmcmc_chord_ok(T,T.triangles[i].daughter_chord[j],T.triangles[i].daughter_type[j])) {
	
	printf("T.triangles[%d].daughter_chord[%d] not ok.\n",i,j);
	return false;

      }

      if (T.triangles[i].daughter_tri[j] < -1 || T.triangles[i].daughter_tri[j] > T.ntri-1) {

	printf("T.triangles[%d].daughter_tri[%d] = %d illegal (<-1 or >T.ntri-1 = %d)\n",
	       i,j,T.triangles[i].daughter_tri[j],T.ntri-1);

	return false;

      }

      if (T.triangles[i].daughter_tri[j] != -1) { /* This is a real daughter, enter in the daughterof list. */

	if (daughterof[T.triangles[i].daughter_tri[j]] != -1) { /* Already used! Trouble. */

	  printf("Triangle %d is a daughter of BOTH triangle %d and triangle %d.\n",
		 T.triangles[i].daughter_tri[j],i,daughterof[T.triangles[i].daughter_tri[j]]);
	  return false;
	  
	}
	
	daughterof[T.triangles[i].daughter_tri[j]] = i;

      }	

    }

    /* Check orientation. The parent chord is stored REVERSED, so we should check it's 0th vt against 
       the 0th vert of daughter[0]. */
    
    if (i > 0) { /* The root triangle stores the parent chord (which is really an edge) in positive orientation */

      if (tsmcmc_vert_chord(T,0,T.triangles[i].parent_chord,T.triangles[i].parent_type) != 
	  tsmcmc_vert_chord(T,0,T.triangles[i].daughter_chord[0],T.triangles[i].daughter_type[0])) {
	
	printf("orientation mismatch at daughter_chord 0/parent_chord transition of triangle %d\n",i);
	printf("daughter 0:"); 
	tsmcmc_print_chord(T,T.triangles[i].daughter_chord[0],T.triangles[i].daughter_type[0]); printf("\n");
	printf("parent_chord:"); 
	tsmcmc_print_chord(T,T.triangles[i].parent_chord,T.triangles[i].parent_type); printf("\n");
	return false;
	
      } 

    } else {

       if (tsmcmc_vert_chord(T,1,T.triangles[i].parent_chord,T.triangles[i].parent_type) != 
	  tsmcmc_vert_chord(T,0,T.triangles[i].daughter_chord[0],T.triangles[i].daughter_type[0])) {
	
	printf("orientation mismatch at daughter_chord 0/parent_chord transition of triangle %d\n",i);
	printf("daughter 0:"); 
	tsmcmc_print_chord(T,T.triangles[i].daughter_chord[0],T.triangles[i].daughter_type[0]); printf("\n");
	printf("parent_chord:"); 
	tsmcmc_print_chord(T,T.triangles[i].parent_chord,T.triangles[i].parent_type); printf("\n");
	return false;
	
      } 

    }
      
    /* Orientations of daughter chords AGREE, so we check tail == head, as expected. */

    if (tsmcmc_vert_chord(T,1,T.triangles[i].daughter_chord[0],T.triangles[i].daughter_type[0]) != 
	tsmcmc_vert_chord(T,0,T.triangles[i].daughter_chord[1],T.triangles[i].daughter_type[1])) {

      printf("orientation mismatch at daughter_chord 0/daughter_chord 1 transition of triangle %d\n",i);
      printf("daughter 0:"); 
      tsmcmc_print_chord(T,T.triangles[i].daughter_chord[0],T.triangles[i].daughter_type[0]); printf("\n");
      printf("daughter 1:"); tsmcmc_print_chord(T,T.triangles[i].daughter_chord[1],T.triangles[i].daughter_type[1]); printf("\n");

      return false;

    } 

    /* Again, orientation of parent chord is REVERSED, so we check tail of daughter 1 against tail of parent */

    if (i > 0) {

      if (tsmcmc_vert_chord(T,1,T.triangles[i].parent_chord,T.triangles[i].parent_type) != 
	  tsmcmc_vert_chord(T,1,T.triangles[i].daughter_chord[1],T.triangles[i].daughter_type[1])) {
	
	printf("orientation mismatch at daughter_chord 1/parent_chord transition of triangle %d\n",i);
	printf("daughter 1:"); 
	tsmcmc_print_chord(T,T.triangles[i].daughter_chord[1],T.triangles[i].daughter_type[1]); printf("\n");
	printf("parent_chord:"); 
	tsmcmc_print_chord(T,T.triangles[i].parent_chord,T.triangles[i].parent_type); printf("\n");
	
	return false;
	
      } 

    } else {

       if (tsmcmc_vert_chord(T,0,T.triangles[i].parent_chord,T.triangles[i].parent_type) != 
	  tsmcmc_vert_chord(T,1,T.triangles[i].daughter_chord[1],T.triangles[i].daughter_type[1])) {
	
	printf("orientation mismatch at daughter_chord 1/parent_chord transition of triangle %d\n",i);
	printf("daughter 1:"); 
	tsmcmc_print_chord(T,T.triangles[i].daughter_chord[1],T.triangles[i].daughter_type[1]); printf("\n");
	printf("parent_chord:"); 
	tsmcmc_print_chord(T,T.triangles[i].parent_chord,T.triangles[i].parent_type); printf("\n");
	
	return false;
	
      } 

    }

  } /* The end of the loop through the triangles. */
    
  /* We now check that all but one of the triangles appear as daughter triangles somewhere. */
  
  int roottri = -1;
  
  for(i=0;i<T.ntri;i++) { 
    
    if (daughterof[i] == -1) { /* This should be the root vertex */
      
      if (roottri != -1) { /* Two roots! */
	
	printf("BOTH triangle %d and triangle %d don't appear as a daughter.\n",
	       i,roottri);
	return false;
	
      } else {
	
	roottri = i;
	
      }
      
    }
    
  }
  
  free(daughterof);
  
  /* These are all of the checks that I can think of right now, but there's probably an
     orientation check that could also be done. */
  
  return true;
  
}


void tsmcmc_triangulation_free(tsmcmc_triangulation_t T)
{
  if (T.triangles != NULL) { free(T.triangles); }
  T.triangles = NULL;
  T.ntri = 0;

  if (T.diags != NULL) { free(T.diags); }
  T.diags = NULL;
  T.ndiags = 0;

  if (T.edges != NULL) { free(T.edges); }
  T.edges = NULL;
  T.nedges = 0;

}

bool     tsmcmc_triangle_from_edgelengths(plc_vector *a, plc_vector *b, plc_vector *c,
					  double A, double B, double C)

/* Construct a triangle from edgelengths, or return false if you can't. */

/* 
         c
        / \
       /   \
      B     A
     /       \
    /         \
   a ----C-----b
 
*/

{
  /* We need first to check that the triangle inequalities hold in order to actually embed 
     the triangle. */

  if (A < 0 || B < 0 || C < 0) { return false; }
  if (A + B < C || B + C < A || C + A < B) { return false; }

  /* Now we need to worry about various special cases. */

  if (A < 1e-8 && B < 1e-8 && C < 1e-8) { /* Very tiny triangle: place at zero */

    *a = plc_build_vect(0,0,0);
    *b = plc_build_vect(0,0,0);
    *c = plc_build_vect(0,0,0);

    return true;

  }

  if (B < 1e-8) { /* Very short leg -- flat triangle. */

    *a = plc_build_vect(0,0,0);
    *b = plc_build_vect(C,0,0);
    *c = plc_build_vect(0,0,0);

    return true;

  } 

  if (C < 1e-8) { /* Very short leg -- flat triangle */

    *a = plc_build_vect(0,0,0);
    *b = plc_build_vect(0,0,0);
    *c = plc_build_vect(B,0,0);

    return true;

  }

  if (A < 1e-8) {

    *a = plc_build_vect(0,0,0);
    *b = plc_build_vect(C,0,0);
    *c = plc_build_vect(C,0,0);

    return true;

  }
  
  double alpha;

  alpha = acos((B*B + C*C - A*A)/(2.0*B*C));

  *a = plc_build_vect(0,0,0);
  *b = plc_build_vect(C,0,0);
  *c = plc_build_vect(B*cos(alpha),B*sin(alpha),0);

  return true;

}

bool tsmcmc_embed_polygon_worker(plCurve *L,tsmcmc_triangulation_t T,
				 double  *edge_lengths,
				 double  *diagonal_lengths,
				 double  *dihedral_angles,
				 int tri,
				 plc_vector parent_normal)

/* Assume that the parent triangle is embedded (that is, those verts in L are filled in) 
   with parent_normal. Embed triangle i using the edge, diagonal, and dihedral data. */

{ 
  if (fabs(plc_norm(parent_normal) - 1.0) >= 1e-8) { 
  
    printf("tsmcmc_embed_polygon_worker: parent_normal has norm %g != 1.0 (tri = %d of %d)\n",plc_norm(parent_normal),tri,T.ntri);
    return false;

  }
  
  //assert(fabs(plc_norm(parent_normal) - 1.0) < 1e-8); /* The parent normal is supposed to be unit. */

  plc_vector current_normal;

  /* The first job is to come up with the normal for the new triangle. If the previous edge is nonzero, 
     this is well-defined; if not, we can just pick one. */

  plc_vector frame[3];

  frame[0] = parent_normal;

  frame[1] = plc_vect_diff(L->cp[0].vt[tsmcmc_vert_chord(T,0,
							 T.triangles[tri].parent_chord,
							 T.triangles[tri].parent_type)],
			   L->cp[0].vt[tsmcmc_vert_chord(T,1,
							 T.triangles[tri].parent_chord,
							 T.triangles[tri].parent_type)]);
  /* Keep in mind that the parent chord is stored in reverse orientation. */

  bool ok;
  frame[1] = plc_normalize_vect(frame[1],&ok);
  double PI = 3.141592653589793;

  if (!ok) {

    current_normal = parent_normal;

  } else { /* Use the dihedral data. */

    frame[2] = plc_cross_prod(frame[0],frame[1]);

    if (fabs(plc_norm(frame[2]) - 1.0) >= 1e-8) {
      
      printf("tsmcmc_embed_polygon_worker: frame[2] has norm %g != 1.0 (tri = %d of %d)\n",plc_norm(frame[2]),tri,T.ntri);
      printf("                             dot of frame[0], frame[1] = %g\n",plc_dot_prod(frame[0],frame[1]));
      printf("                             norm(frame[0]) = %g, norm(frame[1]) = %g\n",plc_norm(frame[0]),plc_norm(frame[1]));
      return false;

    }

    //assert(fabs(plc_norm(frame[2]) - 1.0) < 1e-8); /* frame[0] and frame[1] should have been orthogonal */

    assert(T.triangles[tri].parent_type == diagonal);

    /* We need to be careful here about the angle between the normals, and the dihedral angle. */
    /* The angle between from the previous normal TO the next normal, denoted theta, is given by: 

       theta = dihedral - PI, 

       since when the dihedral is PI, theta is zero, when the dihedral is 0, theta is -PI,
       and when the dihedral is PI/2, theta is -PI/2. */

    double theta = dihedral_angles[T.triangles[tri].parent_chord] - PI;

    current_normal = plc_vlincomb(cos(theta),frame[0],
				  sin(theta),frame[2]);
  }
    
  /* We now generate a new frame for the current triangle using the new normal. */

  plc_vector newframe[3];

  newframe[0] = frame[1];         /* That's going to be "x" axis, along the parent edge. */
  newframe[2] = current_normal;   /* That's going to be the "z" axis, normal to the triangle */
  newframe[1] = plc_cross_prod(newframe[2],newframe[0]); /* Y = Z x X */

  /* Now we construct a triangle in the x-y plane using the triangle builder. */ 

  plc_vector embedded_tri[3];

  if (!tsmcmc_triangle_from_edgelengths(&embedded_tri[0],&embedded_tri[1],&embedded_tri[2],
					tsmcmc_chord_length(T.triangles[tri].daughter_chord[0],
							    T.triangles[tri].daughter_type[0],
							    diagonal_lengths,edge_lengths),
					tsmcmc_chord_length(T.triangles[tri].daughter_chord[1],
							    T.triangles[tri].daughter_type[1],
							    diagonal_lengths,edge_lengths),
					tsmcmc_chord_length(T.triangles[tri].parent_chord,
							    T.triangles[tri].parent_type,
							    diagonal_lengths,edge_lengths))) {

    return false; /* These diagonals and edgelengths are not compatible with the triangulation */

  }
      
  /* Now we actually embed the free vertex of the triangle. The second ("1") vertex of the parent */
  /* chord is actually the "origin" of the new coordinate system (again, the parent chord is REVERSED) */
  /* We then use the "x" and "y" coordinates of the third vertex of the embedded triangle, together */
  /* with the new frame to place the last vertex in space. */

  plc_vector new_vert;

  new_vert = plc_vect_sum(L->cp[0].vt[tsmcmc_vert_chord(T,1,T.triangles[tri].parent_chord,T.triangles[tri].parent_type)],
			  plc_vlincomb(embedded_tri[2].c[0],newframe[0],
				       embedded_tri[2].c[1],newframe[1]));
  
  L->cp[0].vt[tsmcmc_vert_chord(T,1,T.triangles[tri].daughter_chord[0],T.triangles[tri].daughter_type[0])] = new_vert;

  /* We're done this vertex, so the only thing left to do is recurse to the daughters */

  int i;
  
  for(i=0;i<2;i++) { 

    if (T.triangles[tri].daughter_tri[i] != -1) { 

      if (!tsmcmc_embed_polygon_worker(L,T,edge_lengths,diagonal_lengths,dihedral_angles,
				       T.triangles[tri].daughter_tri[i],current_normal)) {

	return false; /* Pass up a failure from any lower level in the tree. */

      }

    }

  }

  return true;
  
} 
  
  
plCurve *tsmcmc_embed_polygon(tsmcmc_triangulation_t T,double *edge_lengths,
			      double *diagonal_lengths,double *dihedral_angles)

/* Construct a plCurve from moment polytope data. */

{
  plCurve *L;
  int nv = T.nedges,cc = 0;
  bool open = false;

  double min_diag = 1e100, min_edgelength = 1e100;
  int i;
  for(i=0;i<T.ndiags;i++) { min_diag = (diagonal_lengths[i] < min_diag) ? diagonal_lengths[i] : min_diag; }
  for(i=0;i<T.nedges;i++) { min_edgelength = (edge_lengths[i] < min_edgelength) ? edge_lengths[i] : min_edgelength; }

  //if (min_diag < 1e-5 || min_edgelength < 1e-5) { 

    /* We can't embed data where there's a zero length diagonal. */
    //return NULL;

  //}
  
  L = plc_new(1,&nv,&open,&cc);

  /* Start by embedding the root triangle. */
  
  if (!tsmcmc_triangle_from_edgelengths(				   
					&L->cp[0].vt[tsmcmc_vert_chord(T,0,T.triangles[0].parent_chord,
								       T.triangles[0].parent_type)], 
					
					&L->cp[0].vt[tsmcmc_vert_chord(T,1,T.triangles[0].parent_chord,
								       T.triangles[0].parent_type)],
					
					&L->cp[0].vt[tsmcmc_vert_chord(T,1,T.triangles[0].daughter_chord[0],
								       T.triangles[0].daughter_type[0])],
					
					tsmcmc_chord_length(T.triangles[0].daughter_chord[0],
							    T.triangles[0].daughter_type[0],
							    diagonal_lengths,edge_lengths),
					
					tsmcmc_chord_length(T.triangles[0].daughter_chord[1],
							    T.triangles[0].daughter_type[1],
							    diagonal_lengths,edge_lengths),
					
					tsmcmc_chord_length(T.triangles[0].parent_chord,T.triangles[0].parent_type,
							    diagonal_lengths,edge_lengths))) {

    plc_free(L);
    return NULL; /* The given diagonal and edge lengths are incompatible with the triangulation. */

  }
      
  plc_vector root_normal = {{0,0,1}}; /* Normal to first triangle */

  for(i=0;i<2;i++) { 

    if (T.triangles[0].daughter_tri[i] != -1) { /* If there is a daughter in this direction */

      if (!tsmcmc_embed_polygon_worker(L,T,edge_lengths,diagonal_lengths,dihedral_angles,
				       T.triangles[0].daughter_tri[i],root_normal)) { 

	plc_free(L);
	return NULL;

      }

    }

  }

  plc_fix_wrap(L);
  return L;

}

void tsmcmc_compute_diagonals(plCurve *L,tsmcmc_triangulation_t T,double *diagonal_lengths)

{
  int i;
  
  for(i=0;i<T.ndiags;i++) {

    int vtA,vtB;
    
    vtA = T.diags[i].vt[0];
    vtB = T.diags[i].vt[1];

    diagonal_lengths[i] = plc_distance(L->cp[0].vt[vtA],L->cp[0].vt[vtB]);

  }

}

void tsmcmc_compute_edgelengths(plCurve *L,tsmcmc_triangulation_t T,double *edge_lengths) 

{
  int i;

  for(i=0;i<T.nedges;i++) {
    
    int vtA,vtB;
    
    vtA = T.edges[i].vt[0];
    vtB = T.edges[i].vt[1];

    edge_lengths[i] = plc_distance(L->cp[0].vt[vtA],L->cp[0].vt[vtB]);

  }

}

struct dihedral { 
    
    int triA[3];
    int triB[3];

};

void tsmcmc_compute_dihedrals_worker(tsmcmc_triangulation_t T,int tri,struct dihedral *dihedrals) {
  
  /* Record the vertices corresponding to the dihedrals of one or both daughter triangles. */
  
  int i;

  for(i=0;i<2;i++) {

    if (T.triangles[tri].daughter_tri[i] != -1) { 

      assert(T.triangles[tri].daughter_type[i] == diagonal);
      
      struct dihedral *di;
      di = &(dihedrals[T.triangles[tri].daughter_chord[i]]); /* Convenience pointer */

      /* triA is the current triangle. The assignments on the right are in parent -> daughter[0] -> daughter[1] order
	 but we assign them in triA so that triA always starts with the shared chord. */

      int a,b,c; /* These are the verts of the current triangle in order. */

      a = tsmcmc_vert_chord(T,1,T.triangles[tri].daughter_chord[1],T.triangles[tri].daughter_type[1]);
      b = tsmcmc_vert_chord(T,0,T.triangles[tri].daughter_chord[0],T.triangles[tri].daughter_type[0]);
      c = tsmcmc_vert_chord(T,1,T.triangles[tri].daughter_chord[0],T.triangles[tri].daughter_type[0]);

      /* Now if we're talking about the 0th daughter, these should be stored b, c, a */

      if (i==0) { di->triA[0] = b; di->triA[1] = c; di->triA[2] = a; } 
      
      /* if we're talking about the 1st daughter, these should be stored c, a, b */

      else { di->triA[0] = c; di->triA[1] = a; di->triA[2] = b; }

      /* triB is a little harder to get at, since it's the daughter
	 triangle. It's always stored in the same order. */

      int dtri;
      dtri = T.triangles[tri].daughter_tri[i];

      di->triB[0] = tsmcmc_vert_chord(T,1,T.triangles[dtri].daughter_chord[1],T.triangles[dtri].daughter_type[1]);
      di->triB[1] = tsmcmc_vert_chord(T,0,T.triangles[dtri].daughter_chord[0],T.triangles[dtri].daughter_type[0]);
      di->triB[2] = tsmcmc_vert_chord(T,1,T.triangles[dtri].daughter_chord[0],T.triangles[dtri].daughter_type[0]);
      
      /* Finally, recurse on the daughter */
      
      tsmcmc_compute_dihedrals_worker(T,dtri,dihedrals);

    }

  }

}

void tsmcmc_compute_dihedral_angles(plCurve *L,tsmcmc_triangulation_t T,
				    double *dihedral_angles, bool *dihedral_defined)

{
  struct dihedral *dihedrals;

  dihedrals = calloc(T.ndiags,sizeof(struct dihedral));
  assert(dihedrals != NULL); 

  /* Traverse the tree, filling in data in the dihedrals array. */

  tsmcmc_compute_dihedrals_worker(T,0,dihedrals);

  /* Now we can actually compute the dihedrals. */

  int i;
  for(i=0;i<T.ndiags;i++) { 

    /* Now the dihedrals worker has been rigged so that we get triA in the order "shared 1, shared 2, unshared"
       and triB */

    bool ok;
    double diangle = plc_dihedral_angle(L->cp[0].vt[dihedrals[i].triA[1]],
					L->cp[0].vt[dihedrals[i].triA[2]],
					L->cp[0].vt[dihedrals[i].triA[0]],
					L->cp[0].vt[dihedrals[i].triB[2]],&ok);						  
    if (ok) { 

      dihedral_angles[i] = diangle;
      dihedral_defined[i] = true;

    } else { /* This is a degenerate triangle, so simply enter 0 for the dihedral. */

      dihedral_angles[i] = 0;
      dihedral_defined[i] = false;

    }

  }

  free(dihedrals);

}
           
bool tsmcmc_polygon_embedding_ok(plCurve *L, 
				 tsmcmc_triangulation_t T,double *edge_lengths,
				 double *diagonal_lengths, double *dihedral_angles)
/* Checks to see whether edgelengths, diagonallengths, and dihedrals are correct. */

{

  int i;

  /* We start with a bunch of basically elementary checks on the polygon. */

  if (L == NULL) { 
    printf("tsmcmc_polygon_embedding_ok: polygon embedding failed.\n");
    return false; 
  }
  if (L->nc != 1) { 
    printf("tsmcmc_polygon_embedding_ok: wrong number of components.\n");
    return false; 
  }
  if (L->cp[0].nv != T.nedges) {
    printf("tsmcmc_polygon_embedding_ok: number of verts (%d) doesn't match triangulation (%d)\n",
	   L->cp[0].nv,T.nedges);
    return false; 
  }
  if (L->cp[0].open != false) {
    printf("tsmcmc_polygon_embedding_ok: polygon is not closed\n");
    return false; 
  }
    
  /* Now we check the edgelengths. */

  double *new_edge_lengths = calloc(T.nedges,sizeof(double));
  assert(new_edge_lengths != NULL);

  tsmcmc_compute_edgelengths(L,T,new_edge_lengths);

  for(i=0;i<T.nedges;i++) {
    
    if (fabs(new_edge_lengths[i] - edge_lengths[i]) > 1e-8) { 

      printf("tsmcmc_polygon_embedding_ok: edge %d has actual length %g != edge_lengths[%d] == %g\n",
	     i,new_edge_lengths[i],i,edge_lengths[i]);
      return false; 

    }

  }

  free(new_edge_lengths);
  
  /* And the diagonal lengths. */

  double *new_diagonal_lengths = calloc(T.ndiags,sizeof(double));
  assert(new_diagonal_lengths != NULL);

  tsmcmc_compute_diagonals(L,T,new_diagonal_lengths);

  for(i=0;i<T.ndiags;i++) {

    if (fabs(new_diagonal_lengths[i] - diagonal_lengths[i]) > 1e-8) { 
      
      printf("tsmcmc_polygon_embedding_ok: diagonal %d has actual length %g != diagonal_lengths[%d] == %g\n",
	     i,new_diagonal_lengths[i],i,diagonal_lengths[i]);
      return false; 
      
    }

  }

  free(new_diagonal_lengths);

  /* And the dihedral angles */

  double *new_dihedral_angles = calloc(T.ndiags,sizeof(double));
  bool   *dihedral_defined = calloc(T.ndiags,sizeof(bool));

  tsmcmc_compute_dihedral_angles(L,T,new_dihedral_angles,dihedral_defined);

  for(i=0;i<T.ndiags;i++) {

    if (dihedral_defined[i]) { /* We can only check if the dihedral makes sense in the first place. */

      if (plc_angle_dist(new_dihedral_angles[i],dihedral_angles[i]) > 1e-8) { 
	
	printf("tsmcmc_polygon_embedding_ok: dihedral %d has angle %g != dihedral_angles[%d] == %g \n",
	       i,new_dihedral_angles[i],i,dihedral_angles[i]);
	return false; 
	
      }

    }

  }

  free(new_dihedral_angles);
  free(dihedral_defined);
   
  /* We have checked everything */

  return true;

}

struct endpoint { 

  int vt;
  char chord;

};

int endpoint_cmp(const void *A,const void *B) {
  struct endpoint *seA = (struct endpoint *)(A);
  struct endpoint *seB = (struct endpoint *)(B);
  return seA->vt - seB->vt;
}

bool tsmcmc_chords_cross(tsmcmc_chord_t chordA,tsmcmc_chord_t chordB) 
/* Decide whether chordA and chordB cross on the n-gon. The idea is that 
   if we sort both sets of endpoints, then a crossing has to look like A->B->A->B. */
{
  /* If the chords are identical, there's a problem. */

  if ((chordA.vt[0] == chordB.vt[0] && chordA.vt[1] == chordB.vt[1]) ||
      (chordA.vt[0] == chordB.vt[1] && chordA.vt[1] == chordB.vt[0])) { 

    return true;

  }

  /* If they have a single common vertex, they are ok. */

  if (chordA.vt[0] == chordB.vt[0] || chordA.vt[0] == chordB.vt[1] ||
      chordA.vt[1] == chordB.vt[0] || chordA.vt[1] == chordB.vt[1]) { 

    return false; /* Chords with a common vertex don't cross */

  }

  struct endpoint epbuf[4];
  
  epbuf[0].vt = chordA.vt[0]; epbuf[0].chord = 'A';
  epbuf[1].vt = chordA.vt[1]; epbuf[1].chord = 'A';
  epbuf[2].vt = chordB.vt[0]; epbuf[2].chord = 'B';
  epbuf[3].vt = chordB.vt[1]; epbuf[3].chord = 'B';

  qsort(epbuf,4,sizeof(struct endpoint),endpoint_cmp);

  if ((epbuf[0].chord == 'A' && epbuf[1].chord == 'B' && epbuf[2].chord == 'A' && epbuf[3].chord == 'B') ||
      (epbuf[0].chord == 'B' && epbuf[1].chord == 'A' && epbuf[2].chord == 'B' && epbuf[3].chord == 'A')) {

    return true;

  } else {

    return false;

  }

}


bool tsmcmc_chord_system_ok(int nchords,tsmcmc_chord_t *chord)
/* Checks to make sure that no pair of chords cross. */
{
  int i,j;

  for(i=1;i<nchords;i++) { 
    
    for(j=0;j<i;j++) { 

      if (tsmcmc_chords_cross(chord[i],chord[j])) { return false; }

    }

  }

  return true;
}

/* We now introduce a data structure for dealing with triangulations
   called the vertex-incidence table or vitable. The vitable stores,
   by vertex, a sorted list of all the chords and edges incident to
   the vertex in counterclockwise order. We use this to build a
   triangulation from a chord system. */

struct typed_chord {
  
  int                num;
  chordtype_t type;

};

struct vert_incidence {

  int                 nchords;
  struct typed_chord *tchords;

};

struct vitable {

  int                    nedges;
  struct vert_incidence *vi;

};


int tsmcmc_vit_position(tsmcmc_triangulation_t T,int chord,chordtype_t type,struct vitable *vit,int vt);

/* To sort the chords (and edges) incident to a vertex in counterclockwise order, 
   we're going to sort the array of typed chords. This requires a (hidden) pointer
   to the global triangulation and (for convenience) a global vertex. */

tsmcmc_triangulation_t *glob_tri_for_sorting;
int glob_vert_num_for_sorting;

int tchord_cmp(const void *A,const void *B) {

  struct typed_chord *tcA = (struct typed_chord *)(A);
  struct typed_chord *tcB = (struct typed_chord *)(B);

  /* The first thing to take care of is edges, which either come first or last 
     no matter what. */

  if (tcA->type == edge || tcB->type == edge) { 

    /* We can figure out the edge chords which come first and last in
       the sort by knowing the vertex we're sorting at. We have to be
       careful about "before" and after at vertex 0, so there's a
       weird modular arithmetic thing here. */

    int eBefore = (glob_vert_num_for_sorting + glob_tri_for_sorting->nedges - 1) % glob_tri_for_sorting->nedges;
    int eAfter = glob_vert_num_for_sorting;  
    eAfter += 0;

    /* Ok, now we know which edge should come first (eAfter) and last (eBefore) */

    if (tcA->type == edge) { 

      assert(tcA->num == eBefore || tcA->num == eAfter); /* It's either eB or eA */
      return (tcA->num == eBefore) ? 1 /* tcB wins */ : -1  /* tcA wins */;

    } else { /* Must have been that tcB->type == edge */

      assert(tcB->num == eBefore || tcB->num == eAfter); /* Should be either eB or eA */
      return (tcB->num == eBefore) ? -1 /* tcA wins */ : 1 /* tcB wins */;

    }

  }

  /* Ok, now we now that we're comparing diagonal chords to each other. */
  /* We're no longer assuming any particular ordering for the chords, so 
     we'll first have to identify the second vertices in each chord */
  
  tsmcmc_chord_t cA, cB;
  cA = glob_tri_for_sorting->diags[tcA->num];
  cB = glob_tri_for_sorting->diags[tcB->num];

  int svA = (cA.vt[0] == glob_vert_num_for_sorting) ? cA.vt[1] : cA.vt[0];
  int svB = (cB.vt[0] == glob_vert_num_for_sorting) ? cB.vt[1] : cB.vt[0];

  /* Now we know the second vertices for each chord. We want to compare them 
     to the given vertex, so we renumber as if the given vertex was vertex 0. */

  svA -= glob_vert_num_for_sorting; if (svA < 0) { svA += glob_tri_for_sorting->nedges; }
  svB -= glob_vert_num_for_sorting; if (svB < 0) { svB += glob_tri_for_sorting->nedges; }

  return svA - svB;

}

struct vitable tsmcmc_build_vitable(tsmcmc_triangulation_t T)
/* 
   We expect the triangulation to be partially filled in with diags,
   edges, nedges, ntri, and ndiags in place. We use this data to
   create a table of incident chords for each vertex of the polygon,
   which is ordered counterclockwise and includes the edges incident
   to the vertex as first and last entries.

*/
{
  int vt,cd;
  struct vitable vit;

  vit.nedges = T.nedges;
  vit.vi = calloc(T.nedges,sizeof(struct vert_incidence));
  assert(vit.vi != NULL);

  /* First, we scan to count the diagonals incident to each vertex */

  for(vt=0;vt<T.nedges;vt++) { vit.vi[vt].nchords = 2; } /* Every vert is incident to two edges */

  for(cd=0;cd<T.ndiags;cd++) { 

    vit.vi[T.diags[cd].vt[0]].nchords++;
    vit.vi[T.diags[cd].vt[1]].nchords++;

  }

  /* Now scan to allocate buffers and put in edges */

  int *chords_used = calloc(T.nedges,sizeof(int));
  assert(chords_used != NULL);

  for(vt=0;vt<T.nedges;vt++) { 

    vit.vi[vt].tchords = calloc(vit.vi[vt].nchords,sizeof(struct typed_chord));
    assert(vit.vi[vt].tchords != NULL); 

    vit.vi[vt].tchords[0].num = (vt + T.nedges - 1) % T.nedges; /* The edge before */
    vit.vi[vt].tchords[0].type = edge;

    vit.vi[vt].tchords[1].num =  vt; /* The edge after */
    vit.vi[vt].tchords[1].type = edge;

    chords_used[vt] = 2;
    
  }
 
  /* Now scan to fill in diagonal typed chords for each vertex */

  int i;

  for(cd=0;cd<T.ndiags;cd++) { 

    for (i=0;i<2;i++) { 
      
      int vt = T.diags[cd].vt[i];
      vit.vi[vt].tchords[chords_used[vt]].num = cd;
      vit.vi[vt].tchords[chords_used[vt]].type = diagonal;
      chords_used[vt]++;

      assert(chords_used[vt] <= vit.vi[vt].nchords); /* Just make sure we don't overrun the buffer */

    }

  }

  /* Last, sort each of the buffers of typed chords */

  glob_tri_for_sorting = &T;
  
  for(vt=0;vt<T.nedges;vt++) {

    glob_vert_num_for_sorting = vt;
    qsort(vit.vi[vt].tchords,vit.vi[vt].nchords,sizeof(struct typed_chord),tchord_cmp);
    
  }

  /* Now free memory and return the vitable */

  free(chords_used);

  return vit;

}

void tsmcmc_vitable_free(struct vitable vit) 

{
  /* Free memory associated with vitable */
  
  int vt;

  if (vit.vi != NULL) {

    for(vt=0;vt<vit.nedges;vt++) {
      
      if (vit.vi[vt].tchords != NULL) { free(vit.vi[vt].tchords); }
      
    }

    free(vit.vi);
    
  }

}

void tsmcmc_orient_daughters(tsmcmc_triangulation_t T,int tri);


void tsmcmc_build_triangulation_worker(tsmcmc_triangulation_t T,
				       int parent_chord,chordtype_t parent_type,
				       int *triangles_used,struct vitable *vit)

/* Given the parent chord, build the next triangle in the triangulation using the vertex incidence table. */

{
  int tri;
  tri = *triangles_used; (*triangles_used)++;

  T.triangles[tri].parent_chord = parent_chord;
  T.triangles[tri].parent_type =  parent_type;

  /* The vertices vA (end of parent chord in this orientation) and vB
     (start of parent chord in this orienation) depend on whether the
     parent is an edge (and hence stored IN order) or a diagonal (and
     hence stored REVERSED). */
 
  int vA = (parent_type == diagonal) ? T.diags[parent_chord].vt[0] : (parent_chord + 1) % T.nedges;
  int vB = (parent_type == diagonal) ? T.diags[parent_chord].vt[1] : parent_chord;

  /* The first daughter chord is the previous chord to the parent at vertex 0 in the vitable. */

  int vp = tsmcmc_vit_position(T,parent_chord,parent_type,vit,vA);
  assert(vp != -1 && vp != 0); 

  T.triangles[tri].daughter_chord[0] = vit->vi[vA].tchords[vp-1].num;
  T.triangles[tri].daughter_type[0] = vit->vi[vA].tchords[vp-1].type;
  T.triangles[tri].daughter_tri[0] = (T.triangles[tri].daughter_type[0] == edge) ? -1 : 0;

  /* The second daughter chord is the next chord to the parent at vertex 1 in the vitable */

  vp = tsmcmc_vit_position(T,parent_chord,parent_type,vit,vB);
  assert(vp != -1 && vp != vit->vi[vB].nchords-1); 

  T.triangles[tri].daughter_chord[1] = vit->vi[vB].tchords[vp+1].num;
  T.triangles[tri].daughter_type[1] = vit->vi[vB].tchords[vp+1].type;
  T.triangles[tri].daughter_tri[1] = (T.triangles[tri].daughter_type[1] == edge) ? -1 : 0;

  int i;

  /* Fix the orientations of the daughter chords */

  tsmcmc_orient_daughters(T,tri);

  for(i=0;i<2;i++) {

    if (T.triangles[tri].daughter_tri[i] != -1) { 

      T.triangles[tri].daughter_tri[i] = *triangles_used;
      tsmcmc_build_triangulation_worker(T,T.triangles[tri].daughter_chord[i],T.triangles[tri].daughter_type[i],triangles_used,vit);

    }

  }

}

int tsmcmc_vit_position(tsmcmc_triangulation_t T,int chord,chordtype_t type,struct vitable *vit,int vt)

/* Finds the given chord in the vitable or returns -1 if not found */

{
  int i;
  
  for(i=0;i<vit->vi[vt].nchords;i++) { 

    if (vit->vi[vt].tchords[i].type == type && vit->vi[vt].tchords[i].num == chord) { return i; }

  }

  return -1;

}

int tsmcmc_triangulation_chordnum(tsmcmc_triangulation_t T,int a,int b)
/* Finds the index of the chord with endpoints a and b, return -1 if not found */
{
  int i;

  for(i=0;i<T.ndiags;i++) { 

    if ((T.diags[i].vt[0] == a && T.diags[i].vt[1] == b) || (T.diags[i].vt[0] == b && T.diags[i].vt[1] == a)) { 

      return i;

    }

  }

  return -1;

}

void tsmcmc_orient_daughters(tsmcmc_triangulation_t T,int tri) 

/* Orients daughter chords of the triangle correctly,
   keeping in mind that

   a) The orientation of the parent chord is BACKWARDS unless the parent chord 
      is an edge (in which case it should be triangle 0).

   b) Edges are always correctly oriented as chords.

   c) If this chord is a diagonal, it either shares the first or the second vertex, 
      with the parent chord, and this tells you which way it should be oriented. 
*/
{
  tsmcmc_chord_t parent_pos;

  if (T.triangles[tri].parent_type == edge) {  /* Stored FORWARDS, since edges always are */
    parent_pos = T.edges[T.triangles[tri].parent_chord];

  } else { /* The parent is BACKWARDS with respect to this triangle, so reverse it. */

    parent_pos.vt[0] = T.diags[T.triangles[tri].parent_chord].vt[1];
    parent_pos.vt[1] = T.diags[T.triangles[tri].parent_chord].vt[0];

  }

  int i;

  for(i=0;i<2;i++) { 

    if (T.triangles[tri].daughter_type[i] == diagonal) { /* diagonals _CAN_ be reversed, edges can't */

      /* There are two cases where we reverse the diagonal.  

	 1) We include parent_pos.vt[0], meaning we're the SECOND edge, 
	 but not in position 1, meaning we have an incorrect orientation. 
	 
	 2) We include parent_pos.vt[1], meaning we're the FIRST edge, 
	 but not in position 0, meaning we have an incorrect orientation. 
	 
      */

      int cd = T.triangles[tri].daughter_chord[i]; /* Just for readability */
      
      if (parent_pos.vt[0] == T.diags[cd].vt[0] || parent_pos.vt[1] == T.diags[cd].vt[1]) { 
	
	int swap; swap = T.diags[cd].vt[0]; T.diags[cd].vt[0] = T.diags[cd].vt[1]; T.diags[cd].vt[1] = swap;
	
      }
      
    }

  }

}

tsmcmc_triangulation_t tsmcmc_generate_triangulation(int nedges,tsmcmc_chord_t *chord)
/* 
   Generate triangulations from a system of chords. 

   The first thing we do is construct a list of chords incident to every vertex in order.
   Once this is done, we can start to construct the triangulation recursively. 

 */
{
  tsmcmc_triangulation_t T;

  /* First, allocate the various arrays. */

  T.ntri   = nedges-2;
  T.ndiags = nedges-3;
  T.nedges = nedges;

  T.triangles = calloc(T.ntri,sizeof(tsmcmc_triangle_t));
  assert(T.triangles != NULL);

  T.diags = calloc(T.ndiags,sizeof(tsmcmc_chord_t));
  assert(T.diags != NULL);

  T.edges = calloc(T.nedges,sizeof(tsmcmc_chord_t));
  assert(T.edges != NULL);

  /* Now fill in the easy part: edges. */

  int i;

  for(i=0;i<T.nedges;i++) {

    T.edges[i].vt[0] = i; T.edges[i].vt[1] = (i+1)%nedges;

  }

  /* The diagonals come from the master list of chords. */

  for(i=0;i<T.ndiags;i++) { 

    T.diags[i] = chord[i];

  }

  if (nedges == 3) { /* This is really a special case. */

    T.triangles[0].parent_chord = 0;
    T.triangles[0].parent_type = edge;
    T.triangles[0].daughter_chord[0] = 1;
    T.triangles[0].daughter_type[0] = edge;
    T.triangles[0].daughter_chord[1] = 2;
    T.triangles[0].daughter_type[1] = edge;
    T.triangles[0].daughter_tri[0] = -1;
    T.triangles[0].daughter_tri[1] = -1;
	
    return T;

  } 
  
  /* The triangles are more difficult and will require real work. */
 
  struct vitable vit;
  int triangles_used = 0;

  vit = tsmcmc_build_vitable(T);
  tsmcmc_build_triangulation_worker(T,0,edge,&triangles_used,&vit);
  tsmcmc_vitable_free(vit);

  return T;

}

void tsmcmc_triangulation_print(FILE *outfile,tsmcmc_triangulation_t T)
/* Prints to a text file to be picked up by Mathematica */

{
  int i;

  for(i=0;i<T.ndiags;i++) { 

    fprintf(outfile,"%d , %d \n",T.diags[i].vt[0],T.diags[i].vt[1]);

  }

}

char *tsmcmc_triangulation_MathematicaForm(tsmcmc_triangulation_t T)

  /* This generates a string giving the Mathematica description of the triangulation */
  
{
  char **diag_strings;
  int i;
  
  diag_strings = calloc(T.ndiags,sizeof(char *));
  assert(diag_strings != NULL);

  for(i=0;i<T.ndiags;i++) { 
    
    diag_strings[i] = calloc(256,sizeof(char));
    assert(diag_strings[i] != NULL);

  }

  int strsize = 32;

  for(i=0;i<T.ndiags;i++) { 

    strsize += snprintf(diag_strings[i],256,"{%d,%d}",T.diags[i].vt[0],T.diags[i].vt[1]) + 5;

  }

  char *outstring = calloc(strsize,sizeof(char)); 
  assert(outstring != NULL);

  char start[16] = "{";
  char comma[16] = ",";
  char end[16] = "}";

  strcat(outstring,start);
  
  for(i=0;i<T.ndiags;i++) { 

    strcat(outstring,diag_strings[i]);
    if(i<T.ndiags-1) { strcat(outstring,comma); }

  }

  strcat(outstring,end);

  for(i=0;i<T.ndiags;i++) { free(diag_strings[i]); }
  free(diag_strings);

  return outstring;
}
  
  
