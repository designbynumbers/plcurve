/*
 * @COPYRIGHT@
 * 
 * Routines for working with vectors.
 *
 * $Id: vector.c,v 1.1 2004-01-21 20:17:46 ashted Exp $
 *
 */

#include "config.h"
#include "ropelength.h"

/**************************************************************/
/*   Basic Linear Algebra Operations                          */
/**************************************************************/

vector plus(vector A,vector B) { /* Returns A + B. */
  vector C;

  C.c[0] = A.c[0] + B.c[0];
  C.c[1] = A.c[1] + B.c[1];
  C.c[2] = A.c[2] + B.c[2];

  return C;
}
  
vector minus(vector A,vector B) { /* Returns A - B. */ 
  vector C;

  C.c[0] = A.c[0] - B.c[0];
  C.c[1] = A.c[1] - B.c[1];
  C.c[2] = A.c[2] - B.c[2];

  return C;
}
  
vector cross(vector A,vector B) { /* Returns A x B. */
  vector C;

  C.c[0] = A.c[1] * B.c[2] - A.c[2] * B.c[1];
  C.c[1] = A.c[2] * B.c[0] - A.c[0] * B.c[2];
  C.c[2] = A.c[0] * B.c[1] - A.c[1] * B.c[0];

  return C;
}

vector scalarmult(double x,vector A) { /* Returns xA. */

  vector C;

  C.c[0] = x * A.c[0];
  C.c[1] = x * A.c[1];
  C.c[2] = x * A.c[2];

  return C;
}

