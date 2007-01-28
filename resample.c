/*
 * Author: Ted Ashton
 *
 * resample.c: Use the spline calls to change the number of vertices in a .vect
 *
 * $Id: resample.c,v 1.1 2007-01-28 21:37:05 ashted Exp $
 *
 */

#include <plCurve.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <argtable2.h>
#include <assert.h>

int main(int argc, char *argv[]) {
  FILE *vectfile;
  int err_num;
  char err_str[80];
  plCurve *L;
  plc_spline *S;
  int cmp, num_verts;
  int *nv;
  double *cmp_lengths;
  double curve_length;

  char revision[20] = "$Revision: 1.1 $";
  char *dollar;

  struct arg_file *filename = arg_file1(NULL,NULL,"<file>",
      "The input VECT file");
  struct arg_file *outname = arg_file0("o","outfile","<file>",
      "The resampled VECT file [stdout]");
  struct arg_int *num = arg_intn("n","num","<int>",0,100,
      "Number of ending vertices [no change]");
  struct arg_dbl *len = arg_dbl0("l","len","<dbl>",
      "Average ending edge length [no change]");
  struct arg_lit *help = arg_lit0("h","help","Print this help and exit");
  struct arg_end *end = arg_end(20);

  void *argtable[] = {help,num,len,outname,filename,end};
  int nerrors;

  dollar = strchr(&revision[1],'$');
  dollar[0] = '\0';
  plc_version(NULL,0);
  fprintf(stderr,"resample v%s\n",&revision[11]);

  if (arg_nullcheck(argtable) != 0) {
    /* NULL entries detected, allocations must have failed. */
    fprintf(stderr,"%s: insufficient memory\n",argv[0]);
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    exit(EXIT_FAILURE);
  }

  /* Parse the command line as defined by argtable[] */
  nerrors = arg_parse(argc,argv,argtable);

  /* special case: '--help' takes preceence over error reporting */
  if (help->count > 0) {
    fprintf(stderr,"Usage: %s ",argv[0]);
    arg_print_syntax(stderr,argtable,"\n");
    arg_print_glossary(stderr,argtable,"  %-25s %s\n");
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return 0;
  }

  /* If the parser returned any errors then display them and exit */
  if (nerrors > 0) {
    /* Display the error details contained in the arg_end struct.*/
    fprintf(stderr,"\n");
    arg_print_errors(stderr,end,argv[0]);
    fprintf(stderr,"\nUsage: %s ",argv[0]);
    arg_print_syntax(stderr,argtable,"\n");
    arg_print_glossary(stderr,argtable,"  %-25s %s\n");
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return 1;
  }

  vectfile = fopen(filename->filename[0],"r");
  assert(vectfile != NULL);
  L = plc_read(vectfile,&err_num,err_str,80);
  fclose(vectfile);
  if (err_num > 0) {
    fprintf(stderr,"Unable to read knot from file:\n%s\n",err_str);
    return(EXIT_FAILURE);
  }
  /* Convert to spline */
  S = plc_convert_to_spline(L,NULL);
  /* Decide how many vertices for each component */
  nv = calloc(L->nc,sizeof(int));
  if (nv == NULL) {
    fprintf(stderr,"%s: insufficient memory\n",argv[0]);
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    exit(EXIT_FAILURE);
  }
  num_verts = 0;
  for (cmp = 0; cmp < L->nc; cmp++) {
    num_verts += L->cp[cmp].nv;
    nv[cmp] = L->cp[cmp].nv;
  }
  if (num->count > 0) {
    if (len->count > 0) {
      fprintf(stderr,"Specified both num and len, ignoring len.\n");
    }
    if (num->count > 1) {
      if (num->count < L->nc) {
        fprintf(stderr,"Need to specify 0, 1 or %d counts for this knot.\n",
            L->nc);
        return(EXIT_FAILURE);
      }
    } else {
      for (cmp = 0; cmp < L->nc; cmp++) {
        nv[cmp] = num->ival[0]*L->cp[cmp].nv/num_verts;
        printf("%d * %d / %d = %d\n", num->ival[0], L->cp[cmp].nv, num_verts,
            nv[cmp]);
      }
    }
  } else if (len->count > 0) {
    cmp_lengths = calloc(L->nc,sizeof(double));
    if (cmp_lengths == NULL) {
      fprintf(stderr,"%s: insufficient memory\n",argv[0]);
      arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
      exit(EXIT_FAILURE);
    }
    curve_length = plc_arclength(L,cmp_lengths);
    for (cmp = 0; cmp < L->nc; cmp++) {
      nv[cmp] = cmp_lengths[cmp]/len->dval[0]+0.5;
      printf("%g / %g = %d\n",cmp_lengths[cmp],len->dval[0],nv[cmp]);
    }
  }
  /* Get rid of the previous one*/
  plc_free(L);
  /* Convert back */
  L = plc_convert_from_spline(S,nv);
  plc_spline_free(S);
  /* And write it out */
  if (outname->count > 0) {
    vectfile = fopen(outname->filename[0],"w");
    assert(vectfile != NULL);
    plc_write(vectfile,L);
    (void)fclose(vectfile);
  } else {
    plc_write(stdout,L);
  }
  return(EXIT_SUCCESS);
}
