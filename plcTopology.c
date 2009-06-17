/*
 *  Routines to work with plCurves as topological objects
 *
 */

/* Copyright 2009 The University of Georgia. */

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

#include <config.h>
#include <plCurve.h>

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#ifdef HAVE_MATH_H
  #include <math.h>
#endif
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

/* Convert plCurve to Millett/Ewing crossing code (documented in plCurve.h) 

   This version is a cheat, writing the first component of the plCurve to 
   disk and calling the utility program compEric.f in order to compute the
   crossing code representation. We then read the representation from disk
   and parse it into the crossing code array. 

   This is scheduled to be replaced with a better version which handles links
   correctly, so this version should only be used for testing. */
			    
char *plc_ccode( plCurve *L )
{
  char *code;
  
  if (L->nc != 1) {

    printf("plc_ccode: L has %d components, but this version only works for 1 comp.\n",
	   L->nc);
    return NULL;

  }

  /* Step 0. Delete the filenames that we will need to create. */
  /*         We use the ISO C "remove" instead of "unlink" for portability. */

  remove("res.pts");
  remove("zkdata");
  remove("zmatrix");
  remove("zpoints");

  /* Step 1. Write a res.pts file in compEric (21.17f) floating point format. */

  int vt;
  FILE *resfile;

  resfile = fopen("res.pts","w");
  assert(resfile != NULL);


  for(vt=0;vt<L->cp[0].nv;vt++) {

    fprintf(resfile,"%21.17f  %21.17f  %21.17f \n",
	    plc_M_clist(L->cp[0].vt[vt]));

  }

  fclose(resfile);

  /* Step 2. Open a pipe to plcompEric and pass in answers 
     to the various questions to make it go. */

  int edges;
  FILE *CEpipe;

  edges = plc_num_edges(L);

  CEpipe = popen("plcompEric > /dev/null","w");
  assert(CEpipe != NULL);

  /* plcompEric expects 

 	# of edges
 	# of reps  (needs to be at least 1)
 	use restart? 1 yes, 0 no (we use 1 because we want to use res.pts)
 	save points in zpoints? 1/0 (no, so 0)
 	random seed, integer (useless, we set to 12345)
 	ropelength (useless here, we set to 30, but don't know why)
 	radius of perturbation (set to 0 to use this config)

	Note: comp516 makes the zmatrix file (also zkdata and zpoints, don't need) */

  fprintf(CEpipe,"%d\n1\n1\n1\n12345\n30\n0.01\n",edges); 
  pclose(CEpipe);

  /* Step 4. Allocate space for the crossing code. */

  int csize = 256*32;
  code = calloc( csize, sizeof(char));

  /* Step 3. Copy the zmatrix file generated by compEric. */

  FILE *zmatrix;
  int   spos = 0;

  zmatrix = fopen("zmatrix","r");
  assert(zmatrix != NULL);

  int zmchar;
  while((zmchar = fgetc(zmatrix)) != EOF) {

    code[spos++] = (char)(zmchar);

    if (spos > csize-10) { 

      csize *= 2; 
      code = realloc(code,csize*sizeof(char));
      assert(code != NULL);

    }

  }

  fclose(zmatrix);

  code[spos] = 0; /* Add the terminating 0 manually */

  /* Step 4. Check and return. */

  printf("Processed %d edge knot successfully. Crossing code is \n %s\n",
	 edges,code);

  return code;

}

  

  
  
