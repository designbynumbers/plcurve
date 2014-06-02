/*************

    polynomials.c : Procedures to convert the output of pllmpoly to 
    a standard polynomial form and manipulate such outputs. 

***************/

#include<plcTopology.h>
#include<polynomials.h>

void homfly_polynomial_free(homfly_polynomial_t **homfly)
{
  if (homfly != NULL) { 

    if (*homfly != NULL) {

      (*homfly)->nmonomials = 0;
      if ((*homfly)->mono != NULL) {

	free((*homfly)->mono);
	(*homfly)->mono = NULL;

      }

      free(*homfly);
      *homfly = NULL;

    } 

  }

}    

int monomial_cmp(const void *A,const void *B)
{
  monomial_t *a,*b;
  a = (monomial_t *)(A);
  b = (monomial_t *)(B);

  if (a->l != b->l) { 

    return a->l - b->l;

  } 

  if (a->m != b->m) {

    return a->m - b->m;

  }

  return a->coeff - b->coeff;

}
 

int initiallpower(char start[]) 

/* Reads forward in the string, counting blocks of stuff separated by whitespace until you get to a [ character. */

{
  int pwr=0;
  int i;
  
  for(i=0;start[i] != '[' && start[i] != 0;i++) { 
    
    if (isspace(start[i])) {for(;isspace(start[i]) && start[i] != 0;i++);} /* Skip whitespace */
    if (start[i] != '[' && start[i] != 0) { 
      for(;!isspace(start[i]) && start[i] != 0;i++) {} 
      pwr--;} /* Skip number, count it. */
    i--;

  }
  return pwr;
}

homfly_polynomial_t *lmpoly_to_polynomial(char *lmout)

/* We convert a polynomial in lmpoly's rather odd format 

   [-1 0 -2 0 [0]]N[0]N1 0 [0]N

   to an array of monomials. The way this works is that the numbers are coefficients
   of monomials in the form l^i m^j. We determine i and j in the following way.

   1) The row inside []'s ([1 0 [3] 0 1]) is the m^0 row.
   2) Rows with a common power of m are delimited by N's.
   3) Inside a row, the 0 power of n has the coefficient in brackets.

   This gives rise to the following algorithm. 

   a) Copy the string to a new buffer. 
   b) Count N's until we find something enclosed by [ [ ] ]. 
      This gives you the initial power of m.
   c) Delete the extra brackets from the copy.
   d) Now scan forward across the polynomial:

         Every time you hit an N, increment the power of m. 
	 Then scan forward to find the current power of l, 
	 which is - the number of spaces until the next [.
	 
   e) Every time you read a number, add it to the polynomial and increment the l power.
   
*/

{

  /* Step a */

  char *lmoutput;

  lmoutput = calloc(strlen(lmout)+10,sizeof(char));
  assert(lmoutput != NULL);
  strcpy(lmoutput,lmout);

  /* Step b */
 
  int i;
  int depth = 0;
  int start = 0,end = 0;
  int itflag = false;
  int initialmpower = 0;

  /* First, we find the start [ and end ] of the middle section. */

  for (i=0;lmoutput[i];i++) {

    if (lmoutput[i] == '[') {

      if (depth == 0) { start = i; depth++; }
      else { itflag = true; }

    }

    if (lmoutput[i] == ']') { 

      depth--; 

      if (itflag && depth == 0) { 

	end = i+1;
	break;  /* Break the enclosing for loop. */

      }

    }

    if (lmoutput[i] == 'N') { 

      initialmpower--;

    }

  }

  lmoutput[start] = ' ';
  lmoutput[end] = ' ';

  /* We now have the initialmpower set up. */
  /* We need to count the actual number of monomials in the strings. This 
     is certainly overestimated the number of spaces in the initial string. */

  int mono_count = 0;
  for(i=0;lmoutput[i];i++) { mono_count += (lmoutput[i] == ' '); }

  /* Now we can create a buffer for the monomials themselves. */

  homfly_polynomial_t *homfly = calloc(1,sizeof(homfly_polynomial_t));
  assert(homfly != NULL);

  homfly->nmonomials = 0;
  homfly->mono = (monomial_t *)(calloc(mono_count,sizeof(monomial_t)));
  
  /* Now we're actually ready to read out some monomials. */

  int mpower;
  int lpower;

  for(i=0,mpower = initialmpower,lpower=initiallpower(&(lmoutput[0])),mono_count=0;lmoutput[i];i++) {

    if(lmoutput[i] == 'N') {

      i++; /* Advance past the N to get initiallpower off on the right foot.*/
      lpower=initiallpower(&(lmoutput[i])); mpower++;
      i--; /* Return a step (we'll advance again in the loop increment) */

    } else if (isspace(lmoutput[i]) || lmoutput[i] == ']' || lmoutput[i] == '[') { 

    } else { /* We're at the lmoutput character of a coeff */

      /* Process the coefficient with sscanf */

      if (sscanf(&(lmoutput[i]),"%d",&(homfly->mono[mono_count].coeff)) != 1) {

	fprintf(stderr,"Couldn't process string %s.\n",&(lmoutput[i]));
	exit(1);

      }

      /* Advance the pointer to the next special char */

      for(;!isspace(lmoutput[i]) && lmoutput[i] != '[' && lmoutput[i] != ']' && lmoutput[i] != 'N';i++) {}
      i--;

      /* We only record this if the coefficient was nonzero. */

      if (homfly->mono[mono_count].coeff != 0.0) {

	/* And fill in the m and l powers */
	
	homfly->mono[mono_count].m = mpower;
	homfly->mono[mono_count].l = lpower;
	
	mono_count++;

      }
      
      /* Finally, increment the l power */

      lpower++;

    }

  }

  free(lmoutput);

  /* Now we put the monomials in standard form by sorting them by total degree, */
  /* then l degree, then m degree */

  qsort(homfly->mono,mono_count,sizeof(monomial_t),monomial_cmp);
  homfly->nmonomials=mono_count;
  return homfly;

}

char *poly_to_mathematica(homfly_polynomial_t *homfly)
{

  char *output;
  output = calloc(128*(homfly->nmonomials) + 1,sizeof(char));
  
  int i;
  for(i=0;i<(homfly->nmonomials);i++) {

    char buffer[256];
    
    if (i > 0) { sprintf(buffer," + "); strcat(output,buffer); }

    sprintf(buffer,"(%d)",homfly->mono[i].coeff); strcat(output,buffer);

    if (homfly->mono[i].l != 0) { 
      sprintf(buffer,"*a^(%d)",homfly->mono[i].l); strcat(output,buffer);
    }

    if (homfly->mono[i].m != 0) { 
      sprintf(buffer,"*z^(%d)",homfly->mono[i].m); strcat(output,buffer);
    }
    
    
  }

  return output;
}

char *poly_to_latex(homfly_polynomial_t *homfly)
{

  char *output;
  output = calloc(128*(homfly->nmonomials) + 1,sizeof(char));
  
  int i;
  for(i=0;i<homfly->nmonomials;i++) {

    char buffer[256];
    
    if (i > 0) { 

      if ((int)(homfly->mono[i].coeff) > 0) { sprintf(buffer," + "); strcat(output,buffer); }
      else { sprintf(buffer," - "); strcat(output,buffer); }

    } else {

      if (homfly->mono[i].coeff < 0) { sprintf(buffer,"-"); strcat(output,buffer); }

    }

    if (abs((int)(homfly->mono[i].coeff)) != 1 || (homfly->mono[i].l == 0 && homfly->mono[i].m == 0)) {
      
      sprintf(buffer,"%d",abs((int)(homfly->mono[i].coeff))); strcat(output,buffer);

    }

    if (homfly->mono[i].l != 0) { 
      sprintf(buffer,"a^{%d}",homfly->mono[i].l); strcat(output,buffer);
    }

    if (homfly->mono[i].m != 0) { 
      sprintf(buffer,"z^{%d}",homfly->mono[i].m); strcat(output,buffer);
    }
    
  }

  return output;
}

char *lmpoly_to_mathematica(char *lmpoly_output)
{
  homfly_polynomial_t *homfly;
  char *output;
  
  homfly = lmpoly_to_polynomial(lmpoly_output);
  output = poly_to_mathematica(homfly);
  homfly_polynomial_free(&homfly);

  return output;
}

char *lmpoly_to_latex(char *lmpoly_output)
{
  homfly_polynomial_t *homfly;
  char *output;
  
  homfly = lmpoly_to_polynomial(lmpoly_output);
  output = poly_to_latex(homfly);
  homfly_polynomial_free(&homfly);

  return output;
}

homfly_polynomial_t *product_polynomial(homfly_polynomial_t *pA, homfly_polynomial_t *pB) 

/* Compute the product of two polynomials. */

{
  monomial_t *rawProduct;
  int nRaw;

  nRaw = pA->nmonomials * pB->nmonomials;
  rawProduct = calloc(nRaw,sizeof(monomial_t));
  assert(rawProduct != 0);

  int i,j,k=0;

  for(i=0;i<pA->nmonomials;i++) {
    for(j=0;j<pB->nmonomials;j++,k++) {
      rawProduct[k].coeff = pA->mono[i].coeff * pB->mono[j].coeff;
      rawProduct[k].m = pA->mono[i].m + pB->mono[j].m;
      rawProduct[k].l = pA->mono[i].l + pB->mono[j].l;
    }
  }

  /* The product is now complete, but may have some duplicate entries. */
  /* We put the rawProduct in standard form by sorting the monomials by total degree, */
  /* then l degree, then m degree */

  qsort(rawProduct,nRaw,sizeof(monomial_t),monomial_cmp);

  /* Now we can scan the rawProduct entries; any entries with the same powers must be */
  /* adjacent after sorting, so we need only read into a new array, advancing to a new */
  /* monomial if the powers of m and n are different and adding coefficients if they are */
  /* the same. */

  homfly_polynomial_t *finishedProduct = calloc(1,sizeof(homfly_polynomial_t));
  assert(finishedProduct != NULL);
  finishedProduct->mono = calloc(nRaw,sizeof(monomial_t));
  assert(finishedProduct->mono != NULL);

  int nFin = 0;

  finishedProduct->mono[0] = rawProduct[0];

  for(nFin=0,i=1;i<nRaw;i++) {

    if (rawProduct[i].m == finishedProduct->mono[nFin].m && rawProduct[i].l == finishedProduct->mono[nFin].l) {

      finishedProduct->mono[nFin].coeff += rawProduct[i].coeff;

    } else {

      nFin++;
      finishedProduct->mono[nFin] = rawProduct[i];

    }

  }

  nFin++;
  free(rawProduct);

  finishedProduct->nmonomials = nFin;
  return finishedProduct;

}

int lmpoly_cmp(const void *A,const void *B)
{
  monomial_t *a,*b;
  a = (monomial_t *)(A);
  b = (monomial_t *)(B);

  if (a->m == b->m) {

    return a->l - b->l;

  } else {

    return a->m - b->m;

  }
}

int lspan(monomial_t *poly,int nmonos) {

  int lowl = INT_MAX, hil = INT_MIN;
  int i;

  for(i=0;i<nmonos;i++) {

    lowl = (poly[i].l < lowl) ? poly[i].l : lowl;
    hil  = (poly[i].l > hil)  ? poly[i].l : hil;

  }

  return hil - lowl;
}

int mspan(monomial_t *poly,int nmonos) {

  int lowm = INT_MAX, him = INT_MIN;
  int i;

  for(i=0;i<nmonos;i++) {

    lowm = (poly[i].m < lowm) ? poly[i].m : lowm;
    him  = (poly[i].m > him)  ? poly[i].m : him;

  }

  return him - lowm;
}


void write_lmpoly_lcoeff(char *buf,int coeff,int l)

{
  char writebuf[256];

  if (l == 0) {
    
    sprintf(writebuf,"[%d] ",coeff); strcat(buf,writebuf);

  } else {

    sprintf(writebuf,"%d ",coeff); strcat(buf,writebuf);

  }

}
    

void write_lmpoly_mline(char *buf,monomial_t *mline,int nmonos)
  
/* Writes the line of coefficients to the buffer buf, surrounding 
   the entire line with brackets if m == 0, the l == 0 coefficient
   with a bracket, and terminating the line with an N. */

{

  /* First, if we end before the 0 power of l, add a last monomial
     for the line to make sure that we seem to end at this power. */

  if (mline[nmonos-1].l < 0) {

    mline[nmonos].coeff = 0;
    mline[nmonos].m = mline[0].m;
    mline[nmonos].l = 0;

    nmonos++;

  }

  /* Similarly, if the first l power is > 0, we need to add a fake starting entry. */

  if (mline[0].l > 0) {

    int i; for(i=nmonos;i>0;i--) { mline[i] = mline[i-1]; }
    
    mline[0].coeff = 0;
    mline[0].m = mline[1].m;
    mline[0].l = 0;

    nmonos++;

  }

  /* Now initialize the line. */

  char writebuf[256];

  if (mline[0].m == 0) {

    sprintf(writebuf,"["); strcat(buf,writebuf);

  }

  /* Now write the first coefficient */

  write_lmpoly_lcoeff(buf,(int)(mline[0].coeff),(int)(mline[0].l));

  int i;

  for(i=1;i<nmonos;i++) {

    /* Fill in any placeholders between the last power of l and us. */

    int lpow = mline[i-1].l+1;
    for(;lpow<mline[i].l;lpow++) { write_lmpoly_lcoeff(buf,0,lpow); }

    /* Now actually write the coefficient */
    
    write_lmpoly_lcoeff(buf,(int)(mline[i].coeff),(int)(mline[i].l));

  }    

  /* This is a total hack, but we just wrote an extra space to buf. Kill it. */
  buf[strlen(buf)-1] = 0;

  if (mline[0].m == 0) {

    sprintf(writebuf,"]"); strcat(buf,writebuf);

  }
  
  sprintf(writebuf,"N"); strcat(buf,writebuf);

}

char *polynomial_to_lmpoly(monomial_t *poly,int nmonos)

/* Returns a buffer giving the polynomial in the (slightly weird) lmpoly format */

{
  monomial_t *sortPoly;
  int i;

  sortPoly = calloc(nmonos+1,sizeof(monomial_t));
  assert(sortPoly != NULL);
  for(i=0;i<nmonos;i++) { sortPoly[i] = poly[i]; }
  

  qsort(sortPoly,nmonos,sizeof(monomial_t),lmpoly_cmp);

  /* We now have a buffer in the sortpoly format. The size of the output string is funny;
     it actually depends on the l and m power values. We can estimate it by computing the 
     span in m and the span in l. */

  char *outbuf;
  outbuf = calloc(lspan(sortPoly,nmonos)*mspan(sortPoly,nmonos)*256,sizeof(char));
  assert(outbuf != 0);

  monomial_t *mline;
  int lsize;
  char writebuf[256];
  
  mline = calloc(nmonos+1,sizeof(monomial_t));
  assert(mline != NULL);
  lsize = 1;
  mline[0] = sortPoly[0];
  
  for(i=1;i<nmonos;i++) {

    if (mline[lsize-1].m == sortPoly[i].m) { /* Part of the same line, push onto it. */

      mline[lsize] = sortPoly[i];
      lsize++;

    } else { /* We have finished the line. */

      /* First, we dispose of the past line. */

      write_lmpoly_mline(outbuf,mline,lsize);

      int mcount = mline[0].m+1;
      for(;mcount<sortPoly[i].m;mcount++) { /* Add placeholder lines. */

	sprintf(writebuf,"[0]N");strcat(outbuf,writebuf);

      }

      /* Now we start a new line */

      mline[0] = sortPoly[i];
      lsize = 1;

    }

  }

  /* We have a last line, so write it. */
  write_lmpoly_mline(outbuf,mline,lsize);

  /* Now add a trailing space */
  sprintf(writebuf," "); strcat(outbuf,writebuf);

  free(mline);
  free(sortPoly);

  return outbuf;
}


char *lmpoly_check(char *lmpoly_output)

/* This procedure, used for debugging, converts lmpoly to polynomial form and back. */

{
  char *cvtback;
  homfly_polynomial_t *homfly;

  homfly = lmpoly_to_polynomial(lmpoly_output);
  cvtback = polynomial_to_lmpoly(homfly->mono,homfly->nmonomials);
  homfly_polynomial_free(&homfly);
  
  return cvtback;
}


  
homfly_polynomial_t *KnotTheory_to_polynomial(char *knottheoryform) 

/* Convert a polynomial in our output form from Mathematica, eg:

   1 a^0 z^0 + -1 a^2 z^0 + 2 a^4 z^0 + -1 a^6 z^0 + 1 a^0 z^2 + -2 a^2 z^2 + 2 a^4 z^2 + -1 a^2 z^4

   to a homfly polynomial data structure. 
*/

{
  /* The first step is to count terms: this is exactly the same as counting pluses and adding 1 */

  int nterms = 1;
  char *pos = knottheoryform;

  for(;*pos != 0;pos++) {
    
    if (*pos == '+') { nterms++; }

  }

  homfly_polynomial_t *homfly = calloc(1,sizeof(homfly_polynomial_t *));
  assert(homfly != NULL);
  homfly->mono = calloc(nterms,sizeof(monomial_t));
  assert(homfly->mono != NULL);
  homfly->nmonomials = nterms;
  
  /* Now we tokenize on the pluses */

  char *monostring;
  int i;
  for(monostring = strtok(knottheoryform,"+"),i=0;i<nterms;i++,monostring = strtok(NULL,"+")) { 

    if (monostring == NULL) {

      pd_error( SRCLOC ,"ran out of monomials early: expected %d, got as far as %d\n",
		NULL, nterms, i);
      return NULL;
    
    }

    if (sscanf(monostring," %d a^%d z^%d ",
	       &(homfly->mono[i].coeff),&(homfly->mono[i].l), &(homfly->mono[i].m)) != 3) { 
      
      pd_error( SRCLOC ,"failed to decode monomial %s", NULL, monostring);
      return NULL;
	
    }

  }

  /* Finally, we sort terms. */

  qsort(homfly->mono,homfly->nmonomials,sizeof(monomial_t),monomial_cmp);

  /* Now we're done. */

  return homfly;

}

bool  polynomials_eq(homfly_polynomial_t *a, homfly_polynomial_t *b) 
/* Compare two polynomials. */
{
  if (a->nmonomials != b->nmonomials) { return false; }

  qsort(a->mono,a->nmonomials,sizeof(monomial_t),monomial_cmp);
  qsort(b->mono,b->nmonomials,sizeof(monomial_t),monomial_cmp);
  
  pd_idx_t i;

  for(i=0;i<a->nmonomials;i++) {

    if (monomial_cmp(&(a->mono[i]),&(b->mono[i])) != 0) { 

      return false;

    }

  }

  return true;

}
