/*

   test_homfly.c : Unit tests for the code in pdcode.c which converts pd codes to M/E codes and computes HOMFLY.


*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#ifdef HAVE_STDIO_H
   #include<stdio.h>
#endif

#ifdef HAVE_ASSERT_H
   #include<assert.h>
#endif

#ifdef HAVE_STRING_H
   #include<string.h>
#endif

#ifdef HAVE_STDINT_H
   #include<stdint.h>
#endif

#ifdef HAVE_STDLIB_H
   #include<stdlib.h>
#endif

#ifdef HAVE_STDBOOL_H
   #include<stdbool.h>
#endif

#ifdef HAVE_TIME_H
   #include<time.h>
#endif

#include<argtable2.h>
#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>

#include<pd_isomorphisms.h>
#include<pd_orientation.h>
#include<polynomials.h>

int PD_VERBOSE=50;

/* We now swap in an open-source replacement for getline. */

/* Always add at least this many bytes when extending the buffer.  */
#define MIN_CHUNK 64

/* Read up to (and including) a TERMINATOR from STREAM into *LINEPTR
   + OFFSET (and null-terminate it). *LINEPTR is a pointer returned from
   malloc (or NULL), pointing to *N characters of space.  It is realloc'd
   as necessary.  Return the number of characters read (not including the
   null terminator), or -1 on error or EOF.  On a -1 return, the caller
   should check feof(), if not then errno has been set to indicate
   the error.  */

int
mygetstr (lineptr, n, stream, terminator, offset)
     char **lineptr;
     size_t *n;
     FILE *stream;
     char terminator;
     int offset;
{
  int nchars_avail;		/* Allocated but unused chars in *LINEPTR.  */
  char *read_pos;		/* Where we're reading into *LINEPTR. */
  int ret;

  if (!lineptr || !n || !stream)
    {
      errno = EINVAL;
      return -1;
    }

  if (!*lineptr)
    {
      *n = MIN_CHUNK;
      *lineptr = malloc (*n);
      if (!*lineptr)
	{
	  errno = ENOMEM;
	  return -1;
	}
    }

  nchars_avail = *n - offset;
  read_pos = *lineptr + offset;

  for (;;)
    {
      int save_errno;
      register int c = getc (stream);

      save_errno = errno;

      /* We always want at least one char left in the buffer, since we
	 always (unless we get an error while reading the first char)
	 NUL-terminate the line buffer.  */

      assert((*lineptr + *n) == (read_pos + nchars_avail));
      if (nchars_avail < 2)
	{
	  if (*n > MIN_CHUNK)
	    *n *= 2;
	  else
	    *n += MIN_CHUNK;

	  nchars_avail = *n + *lineptr - read_pos;
	  *lineptr = realloc (*lineptr, *n);
	  if (!*lineptr)
	    {
	      errno = ENOMEM;
	      return -1;
	    }
	  read_pos = *n - nchars_avail + *lineptr;
	  assert((*lineptr + *n) == (read_pos + nchars_avail));
	}

      if (ferror (stream))
	{
	  /* Might like to return partial line, but there is no
	     place for us to store errno.  And we don't want to just
	     lose errno.  */
	  errno = save_errno;
	  return -1;
	}

      if (c == EOF)
	{
	  /* Return partial line, if any.  */
	  if (read_pos == *lineptr)
	    return -1;
	  else
	    break;
	}

      *read_pos++ = c;
      nchars_avail--;

      if (c == terminator)
	/* Return the line.  */
	break;
    }

  /* Done - NUL terminate and return the number of chars read.  */
  *read_pos = '\0';

  ret = read_pos - (*lineptr + offset);
  return ret;
}

int
mygetline (lineptr, n, stream)
     char **lineptr;
     size_t *n;
     FILE *stream;
{
  return mygetstr (lineptr, n, stream, '\n', 0);
}


char *pdcode_to_ccode(pd_code_t *pd);
char *plc_lmpoly(char *ccode,int timeout); // This is in pllmpoly02.c.
int   monomial_cmp(const void *a, const void *b);

bool millett_ewing_paper_test() { 

  printf("----------------------------------------------\n"
	 "testing crossing code from Millett/Ewing paper\n"
	 "----------------------------------------------\n");

  /* In the original "load-balanced homfly" paper, M/E give
     the example:

     ------>----d\  /a------>--------
     |             / 1 (-)          |
     |           c/  \b             |
     ^           |    |             V
     |          b^    Vc            |
     |            \  /              |
     |              / 2 (-)         |
     |             / \              |
     |            /   \             |
     |          av     ^d           |
     |           |     |            |
     |    b-<----+     +---<--a     |
     |    |                   |     |
     +<a-----c--<------<---d--|--b--+ 
     3 (+)|             4 (+) c
          d                   ^
	  +-------------------+   

    which they claim evaluates to the crossing code:
*/
  
  char correct_ccode[2048] =
    "1-4b2c2b3a\n"
    "2-3b1c1b4a\n"
    "3+1d2a4d4c\n"
    "4+2d1a3d3c\n"
    "\n\n";
  
  printf("testing homfly of code from paper...");
  char *lmpoly_homfly = plc_lmpoly(correct_ccode,300);
  printf("pass (completed computation)\n");

  printf("converting to homfly_polynomial_t...");
  homfly_polynomial_t *poly_homfly = lmpoly_to_polynomial(lmpoly_homfly);
  printf("done.\n");

  homfly_polynomial_t *comp_homfly = calloc(1,sizeof(homfly_polynomial_t));
  
  assert(comp_homfly != NULL);
  comp_homfly->nmonomials = 4;
  comp_homfly->mono = calloc(4,sizeof(monomial_t));
  assert(comp_homfly->mono != NULL);

  /* The expected polynomial (hand computation by Eric Lybrand):
     m^{2} - l^{2} - 1 - l^{-2}"
  */
  comp_homfly->mono[0].coeff = 1;
  comp_homfly->mono[0].l = 0;
  comp_homfly->mono[0].m = 2;

  comp_homfly->mono[1].coeff = -1;
  comp_homfly->mono[1].l = 2;
  comp_homfly->mono[1].m = 0;

  comp_homfly->mono[2].coeff = -1;
  comp_homfly->mono[2].l = 0;
  comp_homfly->mono[2].m = 0;

  comp_homfly->mono[3].coeff = -1;
  comp_homfly->mono[3].l = -2;
  comp_homfly->mono[3].m = 0;
  qsort(comp_homfly->mono,comp_homfly->nmonomials,sizeof(monomial_t),monomial_cmp);
  
  printf("testing against expected value for figure-8...");
  
  if (!polynomials_eq(comp_homfly,poly_homfly)) {

    printf("fail (got %s, expected %s)\n",
	   poly_to_latex(poly_homfly),poly_to_latex(comp_homfly));
    return false;

  } else {
    
    char *lm_latex = poly_to_latex(poly_homfly);
    char *cp_latex = poly_to_latex(comp_homfly);

    printf("pass (got %s, expected %s)\n",lm_latex,cp_latex);
    
    free(lm_latex);
    free(cp_latex);

  }

  homfly_polynomial_free(&poly_homfly);
  homfly_polynomial_free(&comp_homfly);
  free(lmpoly_homfly);

  printf("----------------------------------------------\n"
	 "crossing code from Millett/Ewing paper: pass  \n"
	 "----------------------------------------------\n");

  return true;
 
}

bool trefoil_ccode_test() {

  printf("-------------------------------------------\n"
	 "testing diagrams based on trefoil knot\n"
	 "-------------------------------------------\n");

  pd_code_t *trefoil_pd = pd_build_torus_knot(2,3);

  /*

    According to the docs for build torus knot, this should produce the pd code

   +-----------------------------------------------+
   |                                               |
   +-- 2  ----\    /---0---\   /-- 4 --\   /---2---+
               \0 /         \1/         \2/
	        \            \           \
	       / \          / \         / \
   +-- 5   ---/   \---3  --/   \---1---/   \---5--+
   |                                              |
   +----------------------------------------------+

   which, using the convention from the documentation for lmpoly

       a
       ^
       |
       |
   d<--|-->b (meaning that this can go either way, depending on the orientation of cr)
       |
       |
       c

   So a crossing code representation of a plCurve is a char buffer
   containing lines of the form

   17+2b10c11c31a

   meaning that crossing 17 is a positive crossing

   connected in the a position to the b position of crossing 2,
   connected in the b position to the c position of crossing 10,
   connected in the c position to the c position of crossing 11 and
   connected in the d position to the a position of crossing 31.

   means we should generate the crossing code determined by

   +--------<------------------------<-------------+
   |                                               |
   +------>--c\    /d---->c\   /d---->c\   /d>-----+
               \1 /         \2/         \3/
	        \            \           \
	       / \          / \         / \
   +------>--b/   \a---->-b/   \a>----b/   \a>----+
   |                                              |
   +------<--------------<-------------------<----+

   or

  */

  char correct_ccode[2048] =
    "1+2b3a3d2c\n"
    "2+3b1a1d3c\n"
    "3+1b2a2d1c\n"
    "\n\n";

  printf("testing +trefoil pd code\n");
  printf("computing ccode...");

  char *ccode = pdcode_to_ccode(trefoil_pd);

  printf("done\n");

  //printf("computed ccode:\n\n%s\n",ccode);
  //printf("expected ccode:\n\n%s\n",correct_ccode);

  printf("comparing computed to expected ccode...");

  if (strcmp(ccode,correct_ccode)) {

    printf("FAIL\n\n");
    pd_printf("from PD %PD",trefoil_pd);

    return false;

  }

  printf("pass\n");

  printf("computing HOMFLY using lmpoly...");
  char *lmpoly = plc_lmpoly(ccode,300);
  printf("pass (got %s).\n",lmpoly);

  printf("converting to standard form...");
  char *texform = lmpoly_to_latex(lmpoly);
  printf("pass (got %s)\n",texform);

  free(texform);
  free(lmpoly);
  free(ccode);

  printf("\ntesting -trefoil pd code\n");
  pd_idx_t i;
  for(i=0;i<trefoil_pd->ncross;i++) {

    trefoil_pd->cross[i].sign = PD_NEG_ORIENTATION;

  }

  printf("computing ccode...");
  ccode = pdcode_to_ccode(trefoil_pd);
  printf("done\n");

  /*
    If we switch all of the crossings to negative,
    we get:

   +--------<------------------------<-------------+
   |                                               |
   +------>--d\    /a---->d\   /a---->d\   /a>-----+
               \1 /         \2/         \3/
	        /            /           /
	       / \          / \         / \
   +------>--c/   \b---->-c/   \b>----c/   \b>----+
   |                                              |
   +------<--------------<-------------------<----+

   or
  */

  char minustref_expected[4096] =
    "1-2d2c3b3a\n"
    "2-3d3c1b1a\n"
    "3-1d1c2b2a\n"
    "\n\n";

  //printf("computed ccode:\n\n%s\n",ccode);

  //printf("expected ccode:\n\n%s\n",minustref_expected);

  printf("comparing computed to expected ccode...");

  if (strcmp(ccode,minustref_expected)) {

    printf("FAIL\n\n");
    pd_printf("from PD %PD",trefoil_pd);

    return false;

  }

  printf("pass\n");

  printf("computing HOMFLY using lmpoly...");
  lmpoly = plc_lmpoly(ccode,300);
  printf("pass (got %s).\n",lmpoly);

  printf("converting to standard form...");
  texform = lmpoly_to_latex(lmpoly);
  printf("pass (got %s)\n",texform);

  free(texform);
  free(lmpoly);
  free(ccode);

  printf("\ntesting +-+ unknot pd code\n");
  for(i=0;i<trefoil_pd->ncross;i++) {

    trefoil_pd->cross[i].sign = PD_POS_ORIENTATION;

  }

  trefoil_pd->cross[1].sign = PD_NEG_ORIENTATION;

  printf("computing ccode...");
  ccode = pdcode_to_ccode(trefoil_pd);
  printf("done\n");

  /*
    If we switch ONLY THE MIDDLE CROSSING to negative,
    we get:

   +--------<------------------------<-------------+
   |                                               |
   +------>--c\    /d---->d\   /a---->c\   /d>-----+
               \1 /         \2/         \3/
	        \            /           \
	       / \          / \         / \
   +------>--b/   \a---->-c/   \b>----b/   \a>----+
   |                                              |
   +------<--------------<-------------------<----+

   or
  */

  char pmp_expected[4096] =
    "1+2c3a3d2d\n"
    "2-3c3b1a1d\n"
    "3+1b2b2a1c\n"
    "\n\n";

  //printf("computed ccode:\n\n%s\n",ccode);
  //printf("expected ccode:\n\n%s\n",pmp_expected);

  printf("comparing computed to expected ccode...");

  if (strcmp(ccode,pmp_expected)) {

    printf("FAIL\n\n");
    pd_printf("from PD %PD",trefoil_pd);

    return false;

  }

  printf("pass\n");

  printf("computing HOMFLY using lmpoly...");
  lmpoly = plc_lmpoly(ccode,300);
  printf("pass (got %s).\n",lmpoly);

  printf("converting to standard form...");
  texform = lmpoly_to_latex(lmpoly);
  printf("pass (got %s)\n",texform);

  free(texform);
  free(lmpoly);

  free(ccode);

  pd_code_free(&trefoil_pd);
  printf("-----------------------------------------------\n"
	 "trefoil-based crossing code generation tests: PASS \n"
	 "-----------------------------------------------\n");

  return true;
}

bool test_all_signs(pd_code_t *pd) {

  pd_multidx_t *sign_iterator = pd_new_multidx(pd->ncross,NULL,orientation_ops);
  unsigned int max = pd_multidx_nvals(sign_iterator);
  unsigned int i;

  for(i=0;i<max;i++,pd_increment_multidx(sign_iterator)) {

    /* Actually set the crossings according to the information in the 
       iterator. */

    pd_idx_t k;
    for(k=0;k<sign_iterator->nobj;k++) { 

      pd->cross[k].sign = ((pd_orientation_t *)(sign_iterator->obj[k]))->orient;

    }

    pd_printf("\t testing crsigns %MULTIDX...",NULL,sign_iterator);

    if (!pd_ok(pd)) { 

      pd_printf("FAIL. pd %PD not ok after sign assignment.\n",pd);
      return false;

    }

    char *ccode = pdcode_to_ccode(pd);
    free(ccode); 

    printf("pass.\n");

  }
  
  pd_free_multidx(&sign_iterator);
  return true;

} 

bool unknot_generation_tests() {

  printf("-----------------------------------------------\n"
	 "testing diagrams based on unknots \n"
	 "-----------------------------------------------\n");

  printf("generating codes for unknots with 2-10 + crossings...\n"
	 "\t");

  pd_idx_t k;
  for(k=2;k<11;k++) {

    pd_code_t *pd = pd_build_unknot(k);
    char *ccode = pdcode_to_ccode(pd);
    pd_code_free(&pd);
    free(ccode);
    printf("%d ",k);

  }

  printf("\npass. Generated 2-10 crossing unknot codes.\n");
  printf("\ntesting all crossing signs for 5 crossing diagram...\n\n");

  pd_code_t *pd = pd_build_unknot(5);
  if (test_all_signs(pd)) { 

    printf("\npass (generated all ccodes w/o crashing)\n");


  } else {

    printf("FAIL.\n");
    return false;

  }
      
  pd_code_free(&pd);

  printf("\ntesting all crossing signs for 4 crossing diagram...\n\n");
  pd = pd_build_unknot(4);

  if (test_all_signs(pd)) { 

    printf("\npass (generated all ccodes w/o crashing)\n");

  } else {

    printf("FAIL.\n");
    return false;

  }
      
  pd_code_free(&pd);
  printf("\npass (survived)\n");

  printf("-----------------------------------------------\n"
	 "unknot-based crossing code generation tests: PASS \n"
	 "-----------------------------------------------\n");

  return true;
}

bool test_ccode_conversion() {

  printf("-----------------------------------------------\n"
	 "Conversion of pd_code_t to Millett/Ewing crossing code tests\n");
  printf("tests function \"char *pdcode_to_ccode(pd_code_t *pd)\",\n"
	 "(not part of exposed API)\n"
	 "-----------------------------------------------\n");

  if (!trefoil_ccode_test()) { return false; }
  if (!unknot_generation_tests()) { return false; }
  return true;
}


bool test_unknot_homfly_all_signs(pd_code_t *pd,char *desc) {

  pd_multidx_t *sign_iterator = pd_new_multidx(pd->ncross,NULL,orientation_ops);
  unsigned int max = pd_multidx_nvals(sign_iterator);
  unsigned int i;

  printf("\ttesting HOMFLY for %d crossing %s...",
	 pd->ncross,desc);
  fflush(stdout);

  for(i=0;i<max;i++,pd_increment_multidx(sign_iterator)) {

    /* Actually set the crossings according to the information in the 
       iterator. */

    pd_idx_t k;
    for(k=0;k<sign_iterator->nobj;k++) { 

      pd->cross[k].sign = ((pd_orientation_t *)(sign_iterator->obj[k]))->orient;

    }

    if (!pd_ok(pd)) { 

      pd_printf("FAIL at %MULTIDX. pd %PD not ok after sign assignment.\n",
		pd,sign_iterator);
      return false;

    }

    char *homfly = pd_homfly(pd);
    char unknot_homfly[32] = "1";

    if (strcmp(homfly,unknot_homfly)) {

      pd_printf("\t FAIL at crsigns %MULTIDX\n",NULL,
		sign_iterator);
      printf("\t generated HOMFLY: %s\n\t expected HOMFLY: %s\n",
	     homfly,unknot_homfly);
      return false;

    }

    free(homfly);

  }
  
  printf("pass. (%d crossing sign choices verified).\n",i);
  pd_free_multidx(&sign_iterator);
  return true;

} 

bool test_unknot_homflys() { 

  printf("--------------------------------------------------\n"
	 "test unknot HOMFLY\n"
	 "tests HOMFLY computation on various unknot diagrams\n"
	 "--------------------------------------------------\n");
  pd_idx_t i,j,k;

  printf("testing pd_build_unknot diagrams (simple twists)...\n\n");
  for (i=2;i<7;i++) {
    pd_code_t *pd = pd_build_unknot(i);
    if (!test_unknot_homfly_all_signs(pd,"simple twist")) {
      return false;
    }
    pd_code_free(&pd);
  }

  printf("\ntesting pd_build_unknot_wye (plectonemes)...\n\n");
  char desc[1024];
  for(i=2;i<4;i++) {
    for(j=2;j<4;j++) {
      for(k=2;k<4;k++) {
	sprintf(desc,"%d-%d-%d wye",i,j,k);
	pd_code_t *pd = pd_build_unknot_wye(i,j,k);
	if(!test_unknot_homfly_all_signs(pd,desc)) { 
	  return false;
	}
	pd_code_free(&pd);
      }
    }
  }

  printf("\n--------------------------------------------------\n"
	 "test unknot HOMFLY: PASS\n"
	 "--------------------------------------------------\n");
  return true;

}
    
/* We now write a test where we attempt to load a knot table and compute
   HOMFLYs for it. */

bool rolfsentabletest() 
{
 
  printf("\n--------------------------------------------------\n"
	 "Rolfsen Knot Table/Thistlethwaite Link Table test\n"
	 "--------------------------------------------------\n"); 
  
  printf("Trying to determine srcdir from environment...");
  char *srcdir = getenv("srcdir");
  bool free_srcdir = false;

  if (srcdir == NULL) { 

    printf("fail. (assuming data files local)\n");
    srcdir = calloc(4,sizeof(char));
    sprintf(srcdir,".");
    free_srcdir = true;

  } else {

    printf("pass (srcdir = %s)\n",srcdir);

  }

  char *rolfsentable = calloc(4096,sizeof(char));

  sprintf(rolfsentable,"%s/rolfsentable.txt",srcdir);
  printf("Opening data file %s...",rolfsentable);

  FILE *infile;

  infile = fopen(rolfsentable,"r");
  
  if (infile != NULL) {

    printf("pass\n");

  } else { 

    printf("fail\n");
    return false;

  }

  printf("Loading Mathematica format pd codes...");

  int pd_codes_expected;
  if (!(fscanf(infile,"pdcodes %d \n",&pd_codes_expected) == 1)) { 

    printf("fail. (Couldn't read # of codes in file)");
    return false;

  } 

  pd_code_t **pdbuf;
  pdbuf = calloc(pd_codes_expected,sizeof(pd_code_t *));
  assert(pdbuf != NULL);

  int loaded;
  
  for(loaded = 0;loaded < pd_codes_expected && !feof(infile); loaded++) {

    pdbuf[loaded] = pd_read_KnotTheory(infile);

    if (pdbuf[loaded] == NULL) {

      printf("fail (on pd code %d)\n",loaded);
      return false;

    }

  }

  if (loaded != pd_codes_expected) { 

    printf("fail. (expected %d pd codes, got %d)\n",
	   pd_codes_expected,loaded);
    return false;

  }

  printf("pass (%d pd codes loaded, %d expected).\n",loaded,pd_codes_expected);
  fclose(infile);


  printf("Computing HOMFLYs...");
  clock_t start, end;
  double cpu_time_used;

  char **homflybuf = calloc(loaded,sizeof(char *));
  assert(homflybuf != NULL);
  int computed;

  start = clock();

  for(computed = 0; computed < loaded; computed++) { 

    homflybuf[computed] = pd_homfly(pdbuf[computed]);
    if (homflybuf[computed] == NULL) {

      printf("fail. (pdcode %d doesn't produce HOMFLY)\n",computed);
      return false;

    }

    if (strlen(homflybuf[computed]) == 0) { 

      printf("fail. (produced homfly %s for pdcode %d)\n",
	     homflybuf[computed],computed);
      return false;

    }
    pd_code_free(&(pdbuf[computed]));

  }

  end = clock();
  cpu_time_used =  ((double)(end - start))/CLOCKS_PER_SEC;

  printf("pass (%d homfly polynomials computed in %2.2g sec).\n\n",computed,cpu_time_used);
  free(pdbuf);

  for(computed = 0; computed < loaded; computed++) { 

    free(homflybuf[computed]);

  }
  free(homflybuf);
  free(rolfsentable);

  char *thistlethwaitetable = calloc(4096,sizeof(char));

  sprintf(thistlethwaitetable,"%s/thistlethwaitetable.txt",srcdir);
  printf("Opening data file %s...",thistlethwaitetable);

  infile = fopen(thistlethwaitetable,"r");
  
  if (infile != NULL) {

    printf("pass\n");

  } else { 

    printf("fail\n");
    return false;

  }

  printf("Loading Mathematica format pd codes from file...");

  if (!(fscanf(infile,"pdcodes %d \n",&pd_codes_expected) == 1)) { 

    printf("fail. (Couldn't read # of codes in file)");
    return false;

  } 

  pdbuf = calloc(pd_codes_expected,sizeof(pd_code_t *));
  assert(pdbuf != NULL);

  for(loaded = 0;loaded < pd_codes_expected && !feof(infile); loaded++) {

    pdbuf[loaded] = pd_read_KnotTheory(infile);

    if (pdbuf[loaded] == NULL) {

      printf("fail (on pd code %d)\n",loaded);
      return false;

    }

  }

  if (loaded != pd_codes_expected) { 

    printf("fail. (expected %d pd codes, got %d)\n",
	   pd_codes_expected,loaded);
    return false;

  }

  printf("pass (%d pd codes loaded, %d expected).\n",loaded,pd_codes_expected);
  fclose(infile);


  printf("Computing HOMFLYs...");
  homflybuf = calloc(loaded,sizeof(char *));
  assert(homflybuf != NULL);
  
  start = clock();

  for(computed = 0; computed < loaded; computed++) { 

    homflybuf[computed] = pd_homfly(pdbuf[computed]);

    if (homflybuf[computed] == NULL) {

      printf("fail. (pdcode %d doesn't produce HOMFLY)\n",computed);
      return false;

    }

    if (strlen(homflybuf[computed]) == 0) { 

      printf("fail. (produced homfly %s for pdcode %d)\n",
	     homflybuf[computed],computed);
      return false;

    }
    pd_code_free(&(pdbuf[computed]));

  }

  end = clock();
  cpu_time_used =  ((double)(end - start))/CLOCKS_PER_SEC;

  printf("pass (%d homfly polynomials computed in %2.2g sec).\n",computed,cpu_time_used);

  for(computed = 0; computed < loaded; computed++) { 

    free(homflybuf[computed]);

  }
  free(homflybuf);
  free(pdbuf);
  free(thistlethwaitetable);
  if (free_srcdir) {  free(srcdir); }

  printf("\n"
	 "Note: This is not a test that the HOMFLYPTs are correct,\n"
	 "      just that they don't crash the system, since we\n"
	 "      we don't have an external verification for lmpoly's\n"
	 "      computation of the HOMFLYPT. (KnotTheory uses a \n"
	 "      different skein relation.)\n\n");
  
  printf("\n--------------------------------------------------\n"
	 "Rolfsen Knot Table test: pass\n"
	 "--------------------------------------------------\n\n"); 

  return true;

}
  

    

int main() {

  printf("test_homfly (%s)\n",PACKAGE_STRING);
  printf("--------------------------------------------------------\n"
	 "Unit tests for computing HOMFLY polynomial from pdcode. \n"
	 "========================================================\n");

  if (!millett_ewing_paper_test() ||!test_ccode_conversion() || !rolfsentabletest() || !test_unknot_homflys()) {

    printf("=======================================================\n");
    printf("test_homfly:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================================\n");
    printf("test_homfly:  PASS.\n");
    exit(0);

  }

  return 0;

}
