/* plcRandomPolygon.h  */

/* This header file contains definitions and an API to create random
   open and closed polygons (and polygons with various kinds of constraints).
   It replaces the randompolygon section of the previous plCurve.h header.

   This file contains the exposed portion of the API. More is available
   in the (uninstalled) header file tsmcmc.h, which is distributed with
   plCurve. */

/* Copyright 2014 The University of Georgia */

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

/*@-exportlocal@*/
#ifndef PLCRANDOMPOLYGON_H
#define PLCRANDOMPOLYGON_H

#if (__cplusplus || c_plusplus)
extern "C" {
#endif

  /*
   * Take our chances that stdio.h and stdbool.h exist.  We need stdio to define
   * the FILE type and stdbool to define bool.  If stdbool doesn't exist, one
   * option is to just
   *   typedef int bool;
   * and
   *   #define true 1
   *   #define false 0
   */
#include <stdio.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>  /* We are going to need the gsl_rng type to be defined below. */
#include <plCurve.h>

  /*-

    This package implements several random polygon algorithms:

    1) The "toric symplectic Markov chain Monte Carlo (tsmcmc)"
       algorithm to generate Markov chains of equilateral (or fixed
       edgelength) polygons.

       This is described in some detail in our paper:

       Symplectic Geometry of Random Walks in 3-space by Cantarella and Shonkwiler.

       If you don't want to worry about it, "fan" is always an
       acceptable choice. If you want to play with special-purpose
       triangulations, here are the definitions.

    2) The direct sampling algorithm from 

       A fast direct sampling algorithm for closed equilateral polygons

       by Cantarella, Shonkwiler, and Uehara. This is slower (per sample),
       but produces samples that are known to be uncorrelated.



    We start with the extra stuff for the first algorithm. It depends
    on selecting a triangulation for a closed n-gon to determine which
    form of the moment polytope for the triangulation to use, and you 
    can control this here, if you want to.

   -*/


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

  void tsmcmc_triangulation_polymake(FILE *outfile,tsmcmc_triangulation_t T,double *edge_lengths,bool integral);
  /* Outputs the moment polytope of the triangulation in polymake form */

  /*-

    The equilateral polygon sampling algorithms are not direct sampling
    algorithms; but rather produce a Markov chain of polygons which is
    guaranteed to converge.  As a result, there's no direct API for
    accessing the individual polygons in the chain directly- you just
    compute a single (real-valued) function on each polygon and
    integrate it. The code below computes 95% confidence error bars
    for the result of the integral using the "Geyer IPS estimators".

    Notice that this explicitly works (including convergence and
    error bars) for functions like:

    double f(plCurve *L) {

    if (is_trefoil_knot(L)) { return 1.0; }
    else {return 0.0;}

    }

    Note: you can get to the chain itself using plcrwalk from the command line;
    since the chains are usually very large, they need to be written to
    disk and then post-processed from there.

    -*/


  /* These structures deal with internals of the algorithm. It's ok to
     use the default values. */

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
    FILE   *logfile;               /* Filehandle for log. Set to NULL to turn logging off. */

  } tsmcmc_run_parameters;

  char    *tsmcmc_run_stats_MathematicaForm(tsmcmc_run_stats run_stats);
  char    *tsmcmc_run_params_MathematicaForm(tsmcmc_run_parameters run_params);

  tsmcmc_run_parameters tsmcmc_default_unconfined_parameters(void);
  tsmcmc_run_parameters tsmcmc_default_confined_parameters(void);

  /* 
     Direct Sampler. 

  */

  plCurve *plc_random_equilateral(int n,gsl_rng *rng);
  /* Generate a random equilateral polygon using the CSU algorithm. */

  /*-
    Master functions for integrating over polygon space.
    -*/

    double   tsmcmc_equilateral_expectation(gsl_rng *rng,double integrand(plCurve *L, void *args),
                                            void *args,
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


    double   tsmcmc_confined_equilateral_expectation(gsl_rng *rng,double integrand(plCurve *L, void *args),
                                                     void *args,
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

    double   tsmcmc_fixed_ftc_expectation(gsl_rng *rng,double integrand(plCurve *L, void *args),
                                          void *args,
                                          double failure_to_close,
                                          int max_steps,int max_seconds,
                                          tsmcmc_triangulation_t T,
                                          tsmcmc_run_parameters run_params,
                                          tsmcmc_run_stats *run_stats,
                                          double *error);


  /* This is the "master" driver function for computing an expectation
     over equilateral polygons with a fixed failure to close (one long
     edge of length failure_to_close). These can be triangulated any way
     that is desired, but the edge from vertex 0 to vertex 1 is still
     the long edge. It uses the Geyer ips estimator to compute error
     bars.

     We set the run parameters with the usual run_params struct. Permutations
     only permute the equal edge-length steps.

     Detailed information about the run is returned run_stats (optional,
     set to NULL if you don't care).

     Note that the number of edges is set (implicitly) by the triangulation.

  */

  /*---*/

#endif
