/*
 * @COPYRIGHT@
 * 
 * Routines for working with vectors.
 *
 * $Id: vector.c,v 1.3 2004-05-28 14:33:13 ashted Exp $
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "octrope.h"

/**************************************************************/
/*   Basic Linear Algebra Operations                          */
/**************************************************************/

/* Returns A + B. */
octrope_vector octrope_vplus(octrope_vector A,octrope_vector B) { 
  octrope_vector C;

  C.c[0] = A.c[0] + B.c[0];
  C.c[1] = A.c[1] + B.c[1];
  C.c[2] = A.c[2] + B.c[2];

  return C;
}
  
/* Returns A - B. */ 
octrope_vector octrope_vminus(octrope_vector A,octrope_vector B) { 
  octrope_vector C;

  C.c[0] = A.c[0] - B.c[0];
  C.c[1] = A.c[1] - B.c[1];
  C.c[2] = A.c[2] - B.c[2];

  return C;
}
  
/* Returns A x B. */
octrope_vector octrope_cross(octrope_vector A,octrope_vector B) { 
  octrope_vector C;

  C.c[0] = A.c[1] * B.c[2] - A.c[2] * B.c[1];
  C.c[1] = A.c[2] * B.c[0] - A.c[0] * B.c[2];
  C.c[2] = A.c[0] * B.c[1] - A.c[1] * B.c[0];

  return C;
}

/* Returns xA. */
octrope_vector octrope_scalarmult(double x,octrope_vector A) { 
  octrope_vector C;

  C.c[0] = x * A.c[0];
  C.c[1] = x * A.c[1];
  C.c[2] = x * A.c[2];

  return C;
}

