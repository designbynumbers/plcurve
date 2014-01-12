/*
 * Run the Hopf-link test (if Perl is present)
 *
 * $Id: hopf_test.c,v 1.5 2006-04-18 19:14:36 ashted Exp $
 *
 */

/* Copyright 2004 The University of Georgia. */

/* This file is part of liboctrope.
   
liboctrope is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

liboctrope is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with liboctrope; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef HAVE_SYS_WAIT_H
#include <sys/wait.h>
#endif
#ifdef HAVE_MATH_H
#include <math.h>
#endif
#ifdef HAVE_FLOAT_H
#include <float.h>
#endif
#include "octrope.h"


int main(int argc,char *argv[]) {

  int perl_return;
  int exit_stat;
  plCurve*L;
  int cnt;
  FILE *infile;
  plc_vector v[2];
  double length;

  int ns;
  double shortest;
  octrope_strut *sl;   /* strut list */
  int sl_size;

  if (argc < 2) { 
    /* No file given, assume this was called from the check system and 
       fire up perl, etc. */
    perl_return = system(
      "perl ./test_hopf 3 >/dev/null");
    exit_stat = WEXITSTATUS(perl_return);
  
    if (exit_stat == 127) {
      /* either /bin/sh or perl won't run */
      printf("Warning! Couldn't run hopf test (ok if running from make distcheck)\n");
      return 0;
    } else if (exit_stat != 0) {
      /* Something didn't work right */
      printf("Perl return: %d\n",exit_stat);
      return 0;
    }
    
    perl_return = system(
      "perl ./test_hopf 1234 >/dev/null");
    exit_stat = WEXITSTATUS(perl_return);
  
    if (exit_stat == 127) {
      /* either /bin/sh or perl won't run */
      return 0;
    } else if (exit_stat != 0) {
      /* Something didn't work right */
      printf("Perl return: %d\n",exit_stat);
      return 0;
    }
    return 0;
    
  } else {
    /* We got a filename, find the struts and print them */

    infile = fopen(argv[1],"r");
    if (infile == NULL) {
      printf("%s: Couldn't find input file %s.\n",argv[0],argv[1]);
      exit(-1);
    }

    octrope_error_num = 0;
    L = plc_read(infile,&octrope_error_num,octrope_error_str,80);
    fclose(infile);
    if (octrope_error_num != 0) {
      fprintf(stderr,"%s: Error in call to plc_read: %s\n",argv[0],
        octrope_error_str);
      exit(-1);
    }  

    if (L != NULL) {
      printf("%s: Loaded link from %s.\n",argv[0],argv[1]);

      octrope_set_debug(4);

      sl_size = 2*plc_num_edges(L);
      if ((sl = 
          (octrope_strut *)malloc(sl_size*sizeof(octrope_strut))) != NULL) {
        ns = octrope_struts(L,0,1e-12,sl,sl_size,&shortest,NULL,0);
        if (octrope_error_num != 0) {
          fprintf(stderr,"%s:Error in call to octrope_struts:%s\n",argv[0],
            octrope_error_str);
          free(sl);
          plc_free(L);
          exit(-1);
        }
        printf("Found %d struts, shortest of length %3.3f.\n",ns,shortest);
        if (ns == sl_size) {
          printf(" *** WARNING: may not have all struts, as ns == sl_size.\n");
          exit(-1); /* We know this shouldn't be the cases for the hopf link */
        }
      } else {
        fprintf(stderr,"%s:Unable to allocate space for %d struts.\n",argv[0],
          sl_size);
        sl_size = 0;
        ns = 0;
      }
  
      if (ns > 0) { 
        for (cnt = 0; cnt < ns; cnt++) {
          octrope_strut_ends(L,&sl[cnt],v);
          if (octrope_error_num != 0) {
            fprintf(stderr,"%s:Error in call to octrope_strut_ends:%s\n",
              argv[0], octrope_error_str);
            free(sl);
            plc_free(L);
            exit(-1);
          }
          printf("(%3.3f,%3.3f,%3.3f) - (%3.3f,%3.3f,%3.3f)",
            v[0].c[0],v[0].c[1],v[0].c[2], v[1].c[0],v[1].c[1],v[1].c[2]);
          length = plc_norm(plc_vect_diff(v[0],v[1]));
          printf(" %s %3.3f %d:%d--%d:%d\n",
             (sl[cnt].position[0] == 0 || sl[cnt].position[0] == 1) ?
            ((sl[cnt].position[1] == 0 || sl[cnt].position[1] == 1) ? "End-End" :
             "End-Mid") :
            ((sl[cnt].position[1] == 0 || sl[cnt].position[1] == 1) ? "Mid-End" :
             "Mid-Mid"),length,
              sl[cnt].component[0],sl[cnt].lead_vert[0],
              sl[cnt].component[1],sl[cnt].lead_vert[1]);
        }
      }

      plc_free(L);
      free(sl);
  
      exit(0);
    } else {
      printf("%s: Couldn't load link from %s.\n",argv[0],argv[1]);
      exit(1);
    }
  }
}
