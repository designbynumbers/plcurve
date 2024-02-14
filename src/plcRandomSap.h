/* 

   plcRandomSap.h : 

   This is a header for the rejection sampler for self-avoiding equilateral
   polygons. 

*/

#ifndef __PLC_RANDOM_SAP_H__
#define __PLC_RANDOM_SAP_H__ 1

#ifndef PD_VERBOSE
#define PD_VERBOSE 0
#endif

#include <config.h>
#include"plCurve.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_STDARG_H
#include <stdarg.h>
#endif
#ifdef HAVE_CTYPE_H
#include <ctype.h>
#endif
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#ifdef HAVE_ASSERT_H
#include <assert.h>
#endif
#ifdef HAVE_FLOAT_H
#include <float.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

/* The plcRandomSap function is really inefficient, so implementation is important. */
/* It uses an internal random number generator (not GSL). We need to initialize the */
/* the rng with something like

   uint64_t *xos = plc_xoshiro_init((uint64_t)(time(0)));

*/
uint64_t     *plc_xoshiro_init(uint64_t init);

/* When we are ready to destroy the random number generator, we may use */
void          plc_xoshiro_free(uint64_t *xos);

/* Generate a uniform random number in [0,1]. */ 
inline double plc_xoshiro_uniform(uint64_t *xos);

/* Check if plCurve is self-avoiding in the "string of pearls" model (each
   vertex surrounded by a sphere of radius 1/2. */
bool          plc_is_sap(plCurve *L);

/* Generate a "pearl necklace" self-avoiding polygon by rejection sampling. */
/* The time is expected to be proportional to 1.72^n. */
plCurve      *plc_random_equilateral_closed_self_avoiding_polygon(uint64_t *xos,int n);

#endif
