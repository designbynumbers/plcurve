/*
 * @COPYRIGHT@
 * 
 * Routines for working with vectors.
 *
 * $Id: vector.c,v 1.2 2004-03-02 20:53:42 ashted Exp $
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "octrope.h"

/**************************************************************/
/*   Basic Linear Algebra Operations                          */
/**************************************************************/

octrope_vector plus(octrope_vector A,octrope_vector B) { /* Returns A + B. */
  octrope_vector C;

  C.c[0] = A.c[0] + B.c[0];
  C.c[1] = A.c[1] + B.c[1];
  C.c[2] = A.c[2] + B.c[2];

  return C;
}
  
octrope_vector minus(octrope_vector A,octrope_vector B) { /* Returns A - B. */ 
  octrope_vector C;

  C.c[0] = A.c[0] - B.c[0];
  C.c[1] = A.c[1] - B.c[1];
  C.c[2] = A.c[2] - B.c[2];

  return C;
}
  
octrope_vector cross(octrope_vector A,octrope_vector B) { /* Returns A x B. */
  octrope_vector C;

  C.c[0] = A.c[1] * B.c[2] - A.c[2] * B.c[1];
  C.c[1] = A.c[2] * B.c[0] - A.c[0] * B.c[2];
  C.c[2] = A.c[0] * B.c[1] - A.c[1] * B.c[0];

  return C;
}

octrope_vector scalarmult(double x,octrope_vector A) { /* Returns xA. */
  octrope_vector C;

  C.c[0] = x * A.c[0];
  C.c[1] = x * A.c[1];
  C.c[2] = x * A.c[2];

  return C;
}

