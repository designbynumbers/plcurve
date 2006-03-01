/*
 * Sample program to show the use of liboctrope.a
 *
 * $Id: test_readwrite.c,v 1.14 2006-03-01 15:51:05 ashted Exp $
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
  char err_str[80];
  int err_num;
  int chrs;

  if (argc < 2) {
    printf("usage: test_readwrite <file.vect>\n");
    exit(EXIT_FAILURE);
  }

  plCurve_version(NULL,0);
  printf("test_readwrite: $Revision: 1.14 $\n");

  infile = fopen(argv[1],"r");

  if (infile == NULL) {
    printf("test_readwrite: Couldn't find input file %s.\n",argv[1]);
    exit(EXIT_FAILURE);
  }

  L = plCurve_read(infile,&err_num,err_str,sizeof(err_str));
  (void)fclose(infile);

  if (err_num == 0 && L != NULL) {
    printf("test_readwrite: Loaded plCurve from %s.\n",argv[1]);
  } else {
    printf("test_readwrite: Couldn't load plCurve from %s. \n",argv[1]);
    printf("  %s",err_str);
    exit(EXIT_FAILURE);
  }

  (void)tmpnam(outname);
  outfile = fopen(outname,"w");
  if (outfile == NULL) {
    printf("test_readwrite: Couldn't open temporary file %s.\n",outname);
    plCurve_free(L);
    L = NULL;
    exit(EXIT_FAILURE);
  }

  plCurve_write(outfile,L);
  printf("test_readwrite: Wrote file to temporary file %s.\n",outname);
  (void)fclose(outfile);

  chrs = snprintf(command,sizeof(command),"diff -u -bB %s %s",argv[1],outname);
  if (chrs < (int)sizeof(command)) {
    (void)system(command);
  }

  plCurve_free(L);
  L = NULL;
  return 0;
}
