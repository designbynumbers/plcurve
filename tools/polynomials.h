/*****

	polynomials.h : Definitions for the polynomials source file.


********/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <ctype.h>
#include <limits.h>

#include"ordie.h"

typedef struct monostruct {

	double coeff;
	int l;
	int m;

} monomial;

monomial *lmpoly_to_polynomial(char *lmoutput,int *nmonomials);
char     *lmpoly_to_mathematica(char *lmpoly_output);
char     *lmpoly_to_latex(char *lmpoly_output);

char *poly_to_mathematica(monomial *poly, int *nmonos);
char *poly_to_latex(monomial *poly, int *nmonos);

monomial *product_polynomial(monomial *pA, int nA, monomial *pB, int nB, int *nProduct);
char     *polynomial_to_lmpoly(monomial *poly,int nmonos);

char *lmpoly_check(char *lmpoly_output);

