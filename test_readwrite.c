/*
 * Sample program to show the use of liboctrope.a
 *
 * $Id: test_readwrite.c,v 1.2 2004-10-05 19:15:21 cantarel Exp $
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
#include "octrope.h"

int main(int argc,char *argv[]) {

  octrope_link *L;
  FILE *infile,*outfile;
  char outname[L_tmpnam],command[200 + 2*L_tmpnam],diffname[L_tmpnam];

  if (argc < 2) {

    printf("usage: test_readwrite <file.vect>\n");
    exit(2);

  }

  infile = fopen(argv[1],"r");

  if (infile == NULL) {

    printf("test_readwrite: Couldn't find input file %s.\n",argv[1]);
    exit(2);

  }

  L = octrope_link_read(infile);
  fclose(infile);

  if (L != NULL) {

    printf("test_readwrite: Loaded link from %s.\n",argv[1]);

  } else {

    printf("test_readwrite: Couldn't load link from %s. \n",argv[1]);
    exit(1);

  }

  tmpnam(outname);
  outfile = fopen(outname,"w");

  octrope_link_write(outfile,L);
  printf("test_readwrite: Wrote file to temporary file.\n");

  sprintf(diffname,"diff_results.txt");
  sprintf(command,"diff -bB %s %s > %s",argv[1],outname,diffname);
  
  fclose(outfile);

  return 0;
}