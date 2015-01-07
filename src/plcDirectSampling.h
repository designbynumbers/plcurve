/* plcDirectSampling.h  */

/* This header file contains definitions and an API to create random
   closed equilateral polygons. This file contains the exposed portion
   of the API. More is available in the files kascentpermutation.h and
   plcDirectSampling.c. Note that you're going to have to include have
   <gsl/gsl_rng.h> around to use this package, since we depend on the
   GSL random number generators.
*/

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
#ifndef PLCDIRECTSAMPLING_H
#define PLCDIRECTSAMPLING_H

#if (__cplusplus || c_plusplus)
extern "C" {
#endif


#include <gsl/gsl_rng.h>  /* We are going to need the gsl_rng type to be defined below. */ 
#include <plCurve.h>

  /*-

    This package uses the rejection sampling algorithm from:

    "A fast direct sampling algorithm for equilateral random closed
    polygons" by Cantarella, Shonkwiler and Uehara

    to generate direct samples of equilateral closed polygons. This
    method is less general than the TSMCMC methods provided elsewhere
    in plCurve, but produces perfect samples.

    -*/

  double *plc_sample_hypercube_section(int n, gsl_rng *rng);
  /* R

  
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

  tsmcmc_run_parameters tsmcmc_default_unconfined_parameters();
  tsmcmc_run_parameters tsmcmc_default_confined_parameters();

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
