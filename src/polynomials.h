/*****

	polynomials.h : Definitions for the polynomials source file. This is used for converting the output of lmpoly to a 
	more human-readable form. 


********/ 

#ifndef __POLYNOMIALS_H__
#define __POLYNOMIALS_H__ 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef HAVE_STDBOOL_H
#include <stdbool.h>
#endif

#include <ctype.h>

#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif

#include <assert.h>


typedef struct monomial_struct {

    int coeff;
    int l;  /* Power of l in this monomial */
    int m;  /* Power of m in this monomial */
    
  } monomial_t;

typedef struct polynomial_struct { 

    int nmonomials;
    monomial_t *mono;

} homfly_polynomial_t;

void homfly_polynomial_free(homfly_polynomial_t **homfly);

homfly_polynomial_t *lmpoly_to_polynomial(char *lmoutput);
char                *lmpoly_to_mathematica(char *lmpoly_output);
char                *lmpoly_to_latex(char *lmpoly_output);
homfly_polynomial_t *KnotTheory_to_polynomial(char *knottheoryform);

char *poly_to_mathematica(homfly_polynomial_t *poly);
char *poly_to_latex(homfly_polynomial_t *poly);

homfly_polynomial_t *product_polynomial(homfly_polynomial_t *pA,homfly_polynomial_t *pB);
char     *polynomial_to_lmpoly(monomial_t *poly,int nmonos);

char *lmpoly_check(char *lmpoly_output);

bool  polynomials_eq(homfly_polynomial_t *a, homfly_polynomial_t *b);
/* Compare two polynomials. */


#endif
