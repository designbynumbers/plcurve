/*
 * Author: Jason Cantarella
 *
 * plcurvature.c: Use the plc_MR_curvature call to plot curvature for a .vect
 *
 * $Id: plcurvature.c,v 1.1 2007-03-21 21:07:32 cantarel Exp $
 *
 */

#include <config.h>
#include"plCurve.h"

#ifdef HAVE_STDIO_H
  #include <stdio.h>
#endif
#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif
#ifdef HAVE_FLOAT_H
  #include <float.h>
#endif
#ifdef HAVE_MATH_H
  #include <math.h>
#endif
#ifdef HAVE_STRING_H
  #include <string.h>
#endif
#ifdef HAVE_ASSERT_H
  #include <assert.h>
#endif
#ifdef HAVE_MALLOC_H
  #include <malloc.h>
#endif
#include"argtable2.h"


int main(int argc, char *argv[]) {

  FILE *vectfile;
  FILE *curvfile;
  int err_num;
  char err_str[80];
  plCurve *L;

  int cp, vt;

  char curvfilename[1024];
  char *vectloc;
  double s;

  char revision[20] = "$Revision: 1.1 $";
  char *dollar;

  struct arg_file *filename = arg_file1(NULL,NULL,"<file>",
      "The input VECT file");
  struct arg_file *outname = arg_file0("o","outfile","<file>",
      "The output file of MinRad curvature values [stdout]");
  struct arg_lit *gnuplot = arg_lit0("g","gnuplot","Plot curvature (requires GNUPLOT)");
  struct arg_lit *help = arg_lit0("h","help","Print this help and exit");
  struct arg_end *end = arg_end(20);

  void *argtable[] = {help,gnuplot,outname,filename,end};
  int nerrors;

  dollar = strchr(&revision[1],'$');
  dollar[0] = '\0';
  plc_version(NULL,0);
  fprintf(stderr,"plcurvature v%s\n",&revision[11]);

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
    fprintf(stderr,"computes an approximate curvature for each vertex of a polygon\n");
    fprintf(stderr,"Usage: %s ",argv[0]);
    arg_print_syntax(stderr,argtable,"\n");
    arg_print_glossary(stderr,argtable,"  %-25s %s\n");
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return 0;
  }

  /* If the parser returned any errors then display them and exit */
  if (nerrors > 0) {
    fprintf(stderr,"computes an approximate curvature for each vertex of a polygon\n");
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

  /* Now create a filename for the curvatures file. */

  if (outname->count > 0) { /* The user has given us the filename */

    strcpy(curvfilename,outname->basename[0]);

  } else {  /* We must make it up on our own. */

    strcpy(curvfilename,filename->basename[0]);
    
    if ((vectloc = strstr(curvfilename,".vect")) == NULL) {
      
      vectloc = &(curvfilename[strlen(curvfilename)]);
      
    }

    sprintf(vectloc,".curvature.dat");

  }
  
  /* Now we go ahead and fill the file. */
  
  curvfile = fopen(curvfilename,"w");
  assert(curvfile != NULL);

  fprintf(curvfile,
	  "# Curvature data file for %s output from plcurvature\n"	\
	  "# VECT file has %d components, %d verts \n" \
	  "# cp = component #, vt = vertex # \n" \
	  "# s = arclength from vertex 0 to vertex (vt+1), k = curvature, 1/k radius of curvature \n"
	  "# cp    vt     s     k       1/k     \n",
	  filename->filename[0],
	  L->nc,plc_num_verts(L));

  for(cp=0;cp<L->nc;cp++) {

    s = 0.0;
    
    for(vt=0;vt<L->cp[cp].nv;vt++) {

      s += plc_M_distance(L->cp[cp].vt[vt],L->cp[cp].vt[vt+1]);

      fprintf(curvfile,"  %5d %6d %5g %8g %8g \n",
	      cp,vt,s,plc_MR_curvature(L,cp,vt),1/plc_MR_curvature(L,cp,vt));

    }

    fprintf(curvfile,"#\n");

  }

  /* Now display total curvature information. */

  double tk,*tkc;

  tkc = calloc(sizeof(double),L->nc);
  tk = plc_totalcurvature(L,tkc);

  printf("total curvature (all components): %g (%g pi)\n",tk,tk/3.1415926);
  
  printf("total curvature (by component):");

  for(cp=0;cp<L->nc;cp++) {

    printf(" %g (%g pi) ",tkc[cp],tkc[cp]/3.1415926);

  }
  
  printf("\n");
  free(tkc);
  
  /* Done writing, now we close the files and exit. */

  (void) fclose(curvfile);

  /* Now, if we are supposed to call Gnuplot, please do */

  if (gnuplot->count > 0) {

    FILE *gpinstructions;
    gpinstructions = fopen("cplot.gnuplot","w");
  
    if (gpinstructions == NULL) {

      printf("Couldn't open cplot.gnuplot.\n"); 
      exit(1);

    }

    char plotfilename[4096];
    char pdffilename[4096];
    char texfilename[4096];

    strcpy(plotfilename,curvfilename);

    if ((vectloc = strstr(plotfilename,".dat")) == NULL) { 
      
      sprintf(plotfilename,"curvatureplot.eps");

    } else {

      sprintf(vectloc,".eps");

    }

    strcpy(pdffilename,plotfilename);

    if ((vectloc = strstr(pdffilename,".eps")) == NULL) { 
      
      sprintf(pdffilename,"curvatureplot.pdf");

    } else {

      sprintf(vectloc,".pdf");

    }

    strcpy(texfilename,pdffilename);

    if ((vectloc = strstr(texfilename,".pdf")) == NULL) { 
      
      sprintf(texfilename,"curvatureplot.tex");

    } else {

      sprintf(vectloc,".tex");

    }

    fprintf(gpinstructions,
	    "set output \"%s\" \n"
	    "set term postscript eps \n"

	    "set border 31 lw 2 \n"
	    "set grid back lw 1 \n"

	    "set title \"\" \n"
	    "unset label \n"
	    "unset key \n"

	    "set style line 1 lt -1 lw 1.5 pt 5 ps 1.5 \n"
            "set style data points  \n"
	    "set bars  4 \n",
	    plotfilename);

    fprintf(gpinstructions,"set xtics (");

    double arclengths[4096],totalal;
    totalal = plc_arclength(L,arclengths);

    int cp;
    double cumlength = 0;

    for(cp=0;cp<L->nc;cp++) {

      fprintf(gpinstructions," \"\" %g, \"\" %g, \"\" %g, \"\" %g, \"\" %g ",
	      cumlength, cumlength + 0.2*arclengths[cp], cumlength + 0.4*arclengths[cp], 
	      cumlength + 0.6*arclengths[cp], cumlength + 0.8*arclengths[cp] );

      cumlength += arclengths[cp];

    }

    fprintf(gpinstructions,

	    " )\n"
	    "set ytics ( \"\" 0, \"\" 0.5, \"\" 1, \"\" 1.5, \"\" 2); \n"
	    "set size ratio 0.4 \n"
	    "set xrange [0:%g] \n"
	    "set yrange [0:2] \n"
	    "plot \"%s\" using ($3):($4) with lines lt -1 lw 1 \n"
	    "set term x11 \n"
	    "replot",

	    totalal,curvfilename);

    fclose(gpinstructions);

    printf("Calling GNUPLOT and ps2pdf...\n");
    
    char sys[4096];
    sprintf(sys,"gnuplot cplot.gnuplot; ps2pdf -dEPSCrop %s %s",plotfilename,pdffilename);
    system(sys);


    printf("Writing basic TeX input file with overpic (will probably need to be adjusted by hand)...\n");

    FILE *texfile;
    texfile = fopen(texfilename,"w");

    fprintf(texfile,
	    "\\begin{overpic}{%s}\n"
	    "\\put(0,0){$0$}\n"
	    "\\put(0,25){$0.5$}\n"
	    "\\put(0,50){$1.0$}\n"
	    "\\put(0,75){$1.5$}\n"
	    "\\put(0,100){$2.0$}\n",
	    pdffilename);

    int tag;
    for(tag=0;tag <= 5*L->nc;tag++) {

      fprintf(texfile,
	      "\\put(%g,0){$%4.2f$}\n",
	      tag*(1.0/(5.0*L->nc)),
	      totalal*(tag*(1.0/(5.0*L->nc))));

    }

    fprintf(texfile,
	    "\\end{overpic}\n");

    fclose(texfile);

  }
  plc_free(L);
  return(EXIT_SUCCESS);
}
