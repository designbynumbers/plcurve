/*
 * Sample program to show the use of liboctrope.a
 *
 * $Id: test_readwrite.c,v 1.1 2004-08-31 19:04:02 ashted Exp $
 *
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
