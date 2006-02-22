/*
 * Sample program to show the use of liboctrope.a
 *
 * $Id: test_readwrite.c,v 1.9 2006-02-22 22:54:11 ashted Exp $
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

#include <config.h>
#include <plCurve.h>

#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif

int main(int argc,char *argv[]) {
  plCurve *L;
  FILE *infile,*outfile;
  char outname[L_tmpnam],command[200 + 2*L_tmpnam];

  if (argc < 2) {
    printf("usage: test_readwrite <file.vect>\n");
    exit(2);
  }

  plCurve_version(NULL);
  printf("test_readwrite: $Revision: 1.9 $\n");

  infile = fopen(argv[1],"r");

  if (infile == NULL) {
    printf("test_readwrite: Couldn't find input file %s.\n",argv[1]);
    exit(2);
  }

  L = plCurve_read(infile);
  fclose(infile);

  if (plcl_error_num == 0) {
    printf("test_readwrite: Loaded plCurve from %s.\n",argv[1]);
  } else {
    printf("test_readwrite: Couldn't load plCurve from %s. \n",argv[1]);
    printf("  %s",plcl_error_str);
    exit(1);
  }

  tmpnam(outname);
  outfile = fopen(outname,"w");

  plCurve_write(outfile,L);
  plcl_status_check();
  printf("test_readwrite: Wrote file to temporary file %s.\n",outname);
  fclose(outfile);

  sprintf(command,"diff -u -bB %s %s",argv[1],outname);

  system(command);
  
  return 0;
}
