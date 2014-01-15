/*

  This is a combined header file for the toric symplectic mcmc experiments. 

*/

#ifndef __tsmcmc_h__ 
#define __tsmcmc_h__ 1

#include<plCurve.h>

/* 
   We will need some sort of framework for handling arbitrary
   triangulations.  The key insight here is that a triangulation is a
   tree. This means that each triangulation should be defined in terms
   of up to two daughter triangles. Embedding a triangulation is always
   a journey from root to branch. 

   Orientation is maintained in the following way: 

   1) Edges always have the orientation given by the vertex numbering
      on the polygon. 

   2) Diagonals are POSITIVELY oriented as daughters, and NEGATIVELY
      oriented as parents. 

   The list of daughter chords is maintained in orientation order on 
   the current triangle.

*/

typedef enum chord_enum {diagonal,edge} chordtype_t;

typedef struct chord_struct {

  int vt[2];    /* An chord joins two vertices. Edges ARE chords. */

} tsmcmc_chord_t;

typedef struct triangle_struct { 
  
  int         parent_chord;        /* The parent chord index.                              */
  chordtype_t parent_type;         /* Every chord is either a "diagonal" or "edge" chord.  */

  int         daughter_chord[2];   /* The other two chords in the triangle are daughters.  */
  chordtype_t daughter_type[2];    /* The daughters are either "diagonal" or "edge" chords */
  int         daughter_tri[2];     /* These are indices into the
				      "triangles" array for each
				      daughter chord of type
				      "diagonal", or (-1) for each daughter chord of type "edge". */

  /* The daughter chords are stored in order by the orientation of the triangle so that 
     the order goes "parent_diag->daughter_chord[0]->daughter_chord[1]". The chords are 
     stored oriented so that the daughter chords are oriented positively on this triangle,
     so the "free" vertex is the second vertex of daughter_chord[0] or the first vertex 
     of daughter_chord[1]. */

} tsmcmc_triangle_t;

typedef struct triangulation_struct { 

  int                ntri;      /* This is a master list of (n-2) triangles. */
  tsmcmc_triangle_t *triangles; /* They have an internal tree structure,
                                 but keeping them together in memory
                                 improves performance slightly and
                                 simplifies memory management
                                 greatly. */

  int                ndiags;    /* This is a master list of (n-3) "internal" diagonals. */
  tsmcmc_chord_t    *diags;     /* We use this to index into the array diagonal_lengths. */

  int                nedges;    /* This is a master list of polygon edges. */
  tsmcmc_chord_t    *edges;     /* We use this to index into the array edge_lengths. */

} tsmcmc_triangulation_t;

/*----------------------------------------*/

int tsmcmc_triangulation_chordnum(tsmcmc_triangulation_t T,int a,int b);
/* Finds the index of the chord with endpoints a and b, returns -1 if not found. */

/*----------------------------------------*/

/* We generate triangulations from a much simpler object called a
   chord system.  A chord system for an n-gon is a selection of (n-3)
   chords which do not cross each other. We can auto-generate a 
   triangulation from a given chord system. */

bool tsmcmc_chord_system_ok(int nchords,tsmcmc_chord_t *chord); 
/* Checks to make sure that no pair of chords cross. */

tsmcmc_triangulation_t tsmcmc_generate_triangulation(int nedges,tsmcmc_chord_t *chord);
/* Generate triangulation from a system of nedges-3 chords on a polygon with nedges. */

/*------------------------------------------*/

void tsmcmc_regular_ngon(tsmcmc_triangulation_t T,
			 double **edge_lengths,double **diagonal_lengths,double **dihedral_angles);
/* A standard regular n-gon, used to construct starting configurations. */

void tsmcmc_regular_failure_to_close_ngon(tsmcmc_triangulation_t T,
					  double ftc,
					  double **edge_lengths,double **diagonal_lengths,double **dihedral_angles);
/* Generate planar ngon inscribed in circle with first edgelength ftc and all other edgelengths 1. */

void tsmcmc_equilateral_ngon(gsl_rng *rng,tsmcmc_triangulation_t T,
			     double **edge_lengths,double **diagonal_lengths,double **dihedral_angles);
/* This equilateral unit edge length polygon is used as a default starting point. It is a regular
   n-gon after one moment polytope step and one dihedral step. */

void tsmcmc_failure_to_close_ngon(gsl_rng *rng,tsmcmc_triangulation_t T,double ftc,double **edge_lengths,double **diagonal_lengths,double **dihedral_angles);

void tsmcmc_confined_equilateral_ngon(gsl_rng *rng,tsmcmc_triangulation_t T,double confinement_radius,
				      double **edge_lengths,double **diagonal_lengths,double **dihedral_angles);
/* This is the standard starting configuration for confined sampling. It is based on a "folded triangle". */

/* Triangulation Management Code (in tsmcmc_triangulation.c) */
/*-----------------------------------------------------------*/

tsmcmc_triangulation_t tsmcmc_fan_triangulation(int nedges);
/* Constructs the fan triangulation for an n-edge polygon. */

tsmcmc_triangulation_t tsmcmc_spiral_triangulation(int nedges);
/* Constructs the spiral triangulation which recursively joins adjacent edges with a diagonal. */

tsmcmc_triangulation_t tsmcmc_random_triangulation(gsl_rng *rng,int nedges);
/* Constructs a random triangulation of an n-edge polygon. */

tsmcmc_triangulation_t tsmcmc_teeth_triangulation(int nedges);
/* Constructs a triangulation by "walking" across the polygon. */

void tsmcmc_triangulation_free(tsmcmc_triangulation_t T);

bool tsmcmc_triangulation_ok(tsmcmc_triangulation_t T);
/* Self-tests on a triangulation. */

void tsmcmc_triangulation_print(FILE *outfile,tsmcmc_triangulation_t T);
/* Prints to a text file to be picked up by Mathematica */

char *tsmcmc_triangulation_MathematicaForm(tsmcmc_triangulation_t T);
/* Allocates a string for a Mathematica representation of the triangulation */


/*-----------------------------------------------------------*/

plCurve *tsmcmc_embed_polygon(tsmcmc_triangulation_t T,double *edge_lengths,
			      double *diagonal_lengths,double *dihedral_angles);
/* Construct a space polygon from moment polytope and angle data. */

void tsmcmc_compute_diagonals(plCurve *L,tsmcmc_triangulation_t T,double *diagonal_lengths);

void tsmcmc_compute_edgelengths(plCurve *L,tsmcmc_triangulation_t T,double *edge_lengths);

void tsmcmc_compute_dihedral_angles(plCurve *L,tsmcmc_triangulation_t T,double *dihedral_angles,
				    bool *dihedral_defined);

bool tsmcmc_polygon_embedding_ok(plCurve *L, 
				 tsmcmc_triangulation_t T,double *edge_lengths,
				 double *diagonal_lengths, double *dihedral_angles);
/* Checks to see whether edgelengths, diagonallengths, and dihedrals are correct. */

/* Markov Chain Steps */

typedef enum step_enum {moment_polytope,dihedral,permute} tsmcmc_step_t;

void     tsmcmc_moment_polytope_step(gsl_rng *rng,tsmcmc_triangulation_t T,
				     double *edge_lengths,double *diagonal_lengths);
/*  Make a hit-and-run step in the moment polytope using the
    triangulation T, which alters diagonal_lengths but not edge_lengths. */

void     tsmcmc_confined_moment_polytope_step(gsl_rng *rng,tsmcmc_triangulation_t T,
					      double confinement_radius,
					      double *edge_lengths,double *diagonal_lengths);
/*  Make a hit-and-run step in the confined version of the moment
    polytope using the triangulation T (which is assumed to be the fan
    triangulation) which alters diagonal_lengths. */

bool     tsmcmc_edgelengths_equilateral(tsmcmc_triangulation_t T,double *edge_lengths);
/* Check to make sure all edges are length 1.0 */

bool     tsmcmc_diagonals_ok(tsmcmc_triangulation_t T,double *edge_lengths,double *diagonal_lengths);
/* Check to make sure that the diagonals obey the triangle inequalities from the given triangulation */

bool     tsmcmc_confined_diagonals_ok(tsmcmc_triangulation_t T,double confinement_radius,
				      double *edge_lengths,double *diagonal_lengths,int step);
/* Check to make sure that the diagonals obey the triangle AND CONFINEMENT inequalities for the 
   given triangulation. */

void     tsmcmc_dihedrals_step(gsl_rng *rng,tsmcmc_triangulation_t T,double *dihedral_angles);
/* Reset dihedral angles uniformly. */

void     tsmcmc_edgepermute_step(gsl_rng *rng,tsmcmc_triangulation_t T,
				 double *edge_lengths, double *diagonal_lengths, 
				 double *dihedral_angles);
/* Generate a polygon, permute edges, and then recompute diagonal_lengths and dihedral_angles. */
/* If edgelength 0 is different from the other edgelengths, permutes only edges 1,...,n-1. Doesn't */
/* check for other edgelength combinations at this point, though a very general algorithm would.*/

tsmcmc_step_t tsmcmc_dihedral_diagonal_step(gsl_rng *rng, tsmcmc_triangulation_t T,
					    double *edge_lengths, double *diagonal_lengths,double *dihedral_angles, 
					    double beta);
/* Randomly chooses dihedral or diagonal step. Returns the type of step taken. */

tsmcmc_step_t tsmcmc_confined_dihedral_diagonal_step(gsl_rng *rng, tsmcmc_triangulation_t T,
						     double confinement_radius,
						     double *edge_lengths, double *diagonal_lengths,
						     double *dihedral_angles, 
						     double beta, int moment_polytope_repeat);
/* Randomly chooses dihedral or confined moment polytope step. Returns type of step taken. 
   Here 

   (the fraction) beta of all steps are (confined) moment polytope steps
   all remaining steps are dihedral steps

   and each moment polytope step is actually a sequence of "moment_polytope_repeat" moment polytope steps.
*/

tsmcmc_step_t tsmcmc_dihedral_diagonal_permute_step(gsl_rng *rng, tsmcmc_triangulation_t T,
						    double *edge_lengths, double *diagonal_lengths,
						    double *dihedral_angles, 
						    double beta,double delta,int moment_polytope_repeat);
/* Randomly chooses dihedral, diagonal, or permute step. Returns type of step taken. Here 

   (the fraction) delta of all steps are permutations
   (the fraction) beta of non-permutation steps are moment polytope steps
   all remaining steps are dihedral steps.

   A "moment polytope step" is actually a sequence of "moment_polytope_repeat" moment polytope steps.
   A "permutation step" checks to see if edgelength 0 is different from the other edgelengths (the
   fixed failure-to-close case), and if so, permutes fixing that edge.
*/


/**********************************************************************************************/
/*
  API: These are the functions which are exposed to the end user. 

*/

typedef struct tsmcmc_run_stats_struct {

  int max_lagged_covariance_used; 
  int lagged_covariances_available;

  int dihedral_steps;
  int mp_steps;
  int permute_steps;
  
  double total_seconds;
  double geyer_ips_seconds;

} tsmcmc_run_stats;

typedef struct tsmcmc_run_parameters_struct { 

  int    burn_in;
  double delta;                  /* The fraction of permute steps. */
  double beta;                   /* The fraction of moment polytope steps among non-permute steps. */
  int    moment_polytope_repeat; /* Repeat moment polytope steps this many times. */

  int    log_interval;           /* Log data every this many steps */
  char   logfile_name[4096];     /* Filename to use for logging data (if any) */
  FILE   *logfile;               /* Filehandle for log. */

} tsmcmc_run_parameters;

char    *tsmcmc_run_stats_MathematicaForm(tsmcmc_run_stats run_stats);
char    *tsmcmc_run_params_MathematicaForm(tsmcmc_run_parameters run_params);

double   tsmcmc_equilateral_expectation(gsl_rng *rng,double integrand(plCurve *L),
					int max_steps,int max_seconds,
					tsmcmc_triangulation_t T,
					tsmcmc_run_parameters run_params,
					tsmcmc_run_stats *run_stats,
					double *error);
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


double   tsmcmc_confined_equilateral_expectation(gsl_rng *rng,double integrand(plCurve *L),
						 double confinement_radius, int nedges,
						 int max_steps,int max_seconds,
						 tsmcmc_run_parameters run_params,
						 tsmcmc_run_stats *run_stats,
						 double *error);


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

/*---*/

#endif
