/* 

   ccode_test.c: Make sure the code for converting plCurves to pd_codes is working.

*/

#include<config.h>
#include<plCurve.h>
#include"plcTopology.h"

#ifdef HAVE_MATH_H
  #include<math.h>
#endif
#ifdef HAVE_STDLIB_H
  #include<stdlib.h>
#endif
#ifdef HAVE_TIME_H
  #include<time.h>
#endif
#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

#define VERBOSE 1

int gcd(int x, int y)
{
  int t;

  while (y) {
    t = x;
    x = y;
    y = t % y;
  }
  return(x);
}


plCurve *torusknot(int verts,int p,int q,double major_radius,double minor_radius) 

{

  int cp; 
  int vt;
  
  int *nv;
  bool *open;
  int *cc;

  plCurve *tk;

  cp = gcd(p,q);

  int i;

  nv = calloc(cp,sizeof(int));
  open = calloc(cp,sizeof(bool));
  cc = calloc(cp,sizeof(int));

  for (i=0;i<cp;i++) { 

    nv[i] = verts;
    open[i] = false;
    cc[i] = 0;

  }

  tk = plc_new(cp,nv,open,cc);

  /* We are now prepared to build the actual knot */

  double pofs;
  double theta; 
  double tstep;
  double pi = 3.1415926535897932384626433;
  double pangle,qangle;

  for (i=0;i<cp;i++) {

    for(vt=0,theta=0,pofs=i*(2*pi/cp),tstep = 2*pi/(double)(verts);
	vt<verts;
	vt++,theta+=tstep) {

      plc_vector loc;

      /* Theta is going to run from 0 to 2 pi for this component. If this is 
	 a torus knot, we should circle q times around the major axis while 
	 circling p times around the minor axis. 

	 However, if gcd(p,q) != 1, we should divide both of these by the gcd
	 (or the number of components). So, for instance, in the (2,4) torus
	 link, each component circles 2/2 = 1 time in the q direction and 4/2 = 2
	 times in the p direction. */

      pangle = pofs + (p/cp)*theta;
      qangle = (q/cp)*theta;

      loc = plc_build_vect(
			   major_radius*cos(qangle)*(1+(minor_radius/major_radius)*cos(pangle)),
			   major_radius*sin(qangle)*(1+(minor_radius/major_radius)*cos(pangle)),
			   minor_radius*sin(pangle)
			   );

      tk->cp[i].vt[vt] = loc;

    }
    
  }

  plc_fix_wrap(tk);

  free(nv);
  free(open);
  free(cc);

  return tk;

}

void set_pd_code_from_plCurve_debug(bool val);
  
bool torus_knot_test(gsl_rng *rng,int verts,int q) {

  clock_t start,end;
  double  cpu_time_used;

  printf("(2,%d) torus knot test \n",q);
  printf("------------------------------------------\n");

  printf("generating (2,%d) torus knot at %d verts...",q,verts);  
  plCurve *L = torusknot(verts,q,2,5.0,2.0);
  printf("done\n");

 /*  printf("writing to file torusknot_test.vect...."); */
/*   FILE *outfile = fopen("/Users/cantarel/plcurve/test/torusknot_test.vect","w"); */
/*   plc_write(outfile,L); */
/*   fclose(outfile); */
/*   printf("done\n"); */

  start=clock();
  printf("computing pd_code (random rotation disabled)...");
  
  set_pd_code_from_plCurve_debug(true);
  pd_code_t *projected_pd = pd_code_from_plCurve(rng,L);
  
  end=clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("done (%2.2g sec)\n",cpu_time_used);
  
  printf("building pd_code for (2,%d) torus knot...",q);
  pd_code_t *expected_pd = pd_build_torus_knot(2,q);
  printf("done\n");

  printf("comparing expected and actual pd code...");

  if (!pd_isomorphic(projected_pd,expected_pd)) {

    printf("FAIL (not pd_isomorphic) \n");
    printf("pd code generated from plCurve was:\n");
    printf("-----------------------------------\n");
    pd_write(stdout,projected_pd);
    printf("-----------------------------------\n");
    printf("pd code expected was:\n");
    pd_write(stdout,expected_pd);
    printf("------------------------------------\n");
    printf("torusknot_test: FAIL\n");
    exit(1);

  }

  if (!pd_diagram_isotopic(projected_pd,expected_pd)) {

    printf("FAIL (not pd_diagram_isotopic) \n");
    printf("pd code generated from plCurve was:\n");
    printf("-----------------------------------\n");
    pd_write(stdout,projected_pd);
    printf("-----------------------------------\n");
    printf("pd code expected was:\n");
    pd_write(stdout,expected_pd);
    printf("------------------------------------\n");
    printf("torusknot_test: FAIL\n");
    exit(1);

  }

  plc_free(L);
  pd_code_free(&projected_pd);
  pd_code_free(&expected_pd);

  printf("PASS\n");
  printf("---------------------------------------------\n");
  return true;

}

bool torus_knot_rotation_test(gsl_rng *rng,int verts,int q) {

  clock_t start,end;
  double  cpu_time_used;

  printf("(2,%d) torus knot rotation test \n",q);
  printf("------------------------------------------\n");

  printf("generating (2,%d) torus knot at %d verts...",q,verts);  
  plCurve *L = torusknot(verts,q,2,5.0,2.0);
  printf("done\n");

  printf("looking for fencepost errors by checking each renumbering of L...");

  start=clock();

  int ofs; 
  for(ofs=0;ofs < verts;ofs++) { 

    /* We now shift the vertices in each component of L by one */

    int cp, vt;
    for(cp = 0;cp < L->nc;cp++) { 
      for(vt=0;vt < L->cp[cp].nv;vt++) { 
	L->cp[cp].vt[vt] = L->cp[cp].vt[vt+1];
      }
    }
    plc_fix_wrap(L);
      
    //printf("computing pd_code (random rotation disabled)...");
  
    set_pd_code_from_plCurve_debug(true);
    pd_code_t *projected_pd = pd_code_from_plCurve(rng,L);
 
    //printf("done (%2.2g sec)\n",cpu_time_used);
  
    //printf("building pd_code for (2,%d) torus knot...",q);
    pd_code_t *expected_pd = pd_build_torus_knot(2,q);
    //printf("done\n");

    //printf("comparing expected and actual pd code...");

    if (!pd_isomorphic(projected_pd,expected_pd)) {

      printf("FAIL at offset %d\n",ofs);
      printf("pd code generated from plCurve was:\n");
      printf("-----------------------------------\n");
      pd_write(stdout,projected_pd);
      printf("-----------------------------------\n");
      printf("pd code expected was:\n");
      pd_write(stdout,expected_pd);
      printf("------------------------------------\n");
      printf("torusknot_test: FAIL\n");
      exit(1);
      
    }

    

    pd_code_free(&projected_pd);
    pd_code_free(&expected_pd);
  
  }


  end=clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;
  plc_free(L);

  printf("pass\n   (%2.2g sec total - %2.4g sec per pdcode)\n",
	 cpu_time_used,cpu_time_used/(double)(verts));
  printf("---------------------------------------------\n");
  return true;

}

#define MAXARCS 32

typedef struct arc_presentation_cmp_struct {

  int narcs;
  int page[MAXARCS];
  int startlevel[MAXARCS];

  /* Arc i occupies page[i] of the book, and 
    starts at level startlevel[i].
  */

} arc_presentation_cmp;

#define MAXCMP 12

typedef struct arc_presentation_struct {

  int nc;
  int pages;
  arc_presentation_cmp cp[MAXCMP];

} arc_presentation;

plCurve *plCurve_from_arcpres(arc_presentation *AP) {

  double PI = 3.14159265359;
  int *cc,*nv;
  bool *open;

  cc = calloc(AP->nc,sizeof(int));
  nv = calloc(AP->nc,sizeof(int));
  open = calloc(AP->nc,sizeof(bool));

  int cp;

  for(cp=0;cp<AP->nc;cp++) { 

    nv[cp] = 3*AP->cp[cp].narcs;
    open[cp] = false;
    cc[cp] = 0;

  }

  plCurve *L = plc_new(AP->nc,nv,open,cc);
  
  for(cp=0;cp<L->nc;cp++) { 

    int arc = 0;
    for(arc=0;arc < AP->cp[cp].narcs;arc++) {

      int pageNum = AP->cp[cp].page[arc];

      L->cp[cp].vt[3*arc] = plc_build_vect(5*cos(pageNum*2*PI/((double)(AP->pages))),
					   5*sin(pageNum*2*PI/((double)(AP->pages))),
					   AP->cp[cp].startlevel[arc]);

      L->cp[cp].vt[3*arc+1] = plc_build_vect(5*cos(pageNum*2*PI/((double)(AP->pages))),
					     5*sin(pageNum*2*PI/((double)(AP->pages))),
					     AP->cp[cp].startlevel[(arc+1)%AP->cp[cp].narcs]);

      L->cp[cp].vt[3*arc+2] = plc_build_vect(0,0,AP->cp[cp].startlevel[(arc+1)%AP->cp[cp].narcs]);

    }
      
  }

  plc_fix_wrap(L);
  free(nv);
  free(cc);
  free(open);
  
  return L;

}

plCurve *basic_hopf() {

  arc_presentation AP;

  AP.pages = 4;
  AP.nc = 2;

  AP.cp[0].narcs = 2;
  AP.cp[0].page[0] = 0;       AP.cp[0].page[1] = 2;
  AP.cp[0].startlevel[0] = 0; AP.cp[0].startlevel[1] = 2;

  AP.cp[1].narcs = 2;
  AP.cp[1].page[0] = 1;       AP.cp[1].page[1] = 3;
  AP.cp[1].startlevel[0] = 1; AP.cp[1].startlevel[1] = 3;

  return plCurve_from_arcpres(&AP);

}

plCurve *basic_trefoil() { 

  arc_presentation AP;

  AP.pages = 5;
  AP.nc = 1;
  AP.cp[0].narcs = 5;

  AP.cp[0].page[0] = 1; AP.cp[0].startlevel[0] = 1;
  AP.cp[0].page[1] = 3; AP.cp[0].startlevel[1] = 3;
  AP.cp[0].page[2] = 0; AP.cp[0].startlevel[2] = 5;
  AP.cp[0].page[3] = 2; AP.cp[0].startlevel[3] = 2;
  AP.cp[0].page[4] = 4; AP.cp[0].startlevel[4] = 4;
  
  return plCurve_from_arcpres(&AP);

}

bool arcpresentation_tests(gsl_rng *rng) {

  clock_t start,end;
  double cpu_time_used;

  printf("===============================================\n"
	 "arc presentation tests (nongeneric projections)\n"
	 "===============================================\n");

  printf("generating plCurve for basic hopf link...");
  plCurve *hopf = basic_hopf();
  printf("done\n");

  /* printf("writing test file hopftest.vect..."); */
/*   FILE *outfile; */
/*   outfile = fopen("hopftest.vect","w"); */
/*   plc_write(outfile,hopf); */
/*   fclose(outfile); */
/*   printf("done\n"); */

  start=clock();
  printf("computing pd_code (random rotation disabled)...");
  
  set_pd_code_from_plCurve_debug(true);
  pd_code_t *projected_pd = pd_code_from_plCurve(rng,hopf);
  
  end=clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("done (%2.2g sec)\n",cpu_time_used);

  printf("checking pd_ok...");
  
  if (pd_ok(projected_pd)) { 
    
    printf("pass\n");

  } else { 

    printf("FAIL\n");
    printf("pd code generated from plCurve was:\n");
    printf("-----------------------------------\n");
    pd_write(stdout,projected_pd);
    pd_code_free(&projected_pd);

    exit(1);
  }

  pd_code_free(&projected_pd);

  start=clock();
  printf("computing pd_code (random rotation enabled)...");
  
  set_pd_code_from_plCurve_debug(false);
  projected_pd = pd_code_from_plCurve(rng,hopf);
  
  end=clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("done (%2.2g sec)\n",cpu_time_used);

  printf("checking pd_ok...");
  
  if (pd_ok(projected_pd)) { 
    
    printf("pass\n");
    pd_code_free(&projected_pd);

  } else { 

    printf("FAIL\n");
    printf("pd code generated from plCurve was:\n");
    printf("-----------------------------------\n");
    pd_write(stdout,projected_pd);
    pd_code_free(&projected_pd);

    exit(1);
  }
  plc_free(hopf);

  printf("---------------------------------------\n"
	 "generating plCurve for trefoil...");

  plCurve *tref = basic_trefoil();
  printf("done\n");

 /*  printf("writing test file treftest.vect..."); */
/*   outfile = fopen("treftest.vect","w"); */
/*   plc_write(outfile,tref); */
/*   fclose(outfile); */
/*   printf("done\n"); */

  start=clock();
  printf("computing pd_code (random rotation disabled)...");
  
  set_pd_code_from_plCurve_debug(true);
  projected_pd = pd_code_from_plCurve(rng,tref);
  
  end=clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("done (%2.2g sec)\n",cpu_time_used);

  printf("checking pd_ok...");
  
  if (pd_ok(projected_pd)) { 
    
    printf("pass\n");
    pd_code_free(&projected_pd);

  } else { 

    printf("FAIL\n");
    printf("pd code generated from plCurve was:\n");
    printf("-----------------------------------\n");
    pd_write(stdout,projected_pd);
    pd_code_free(&projected_pd);

    exit(1);
  }

 
  start=clock();
  printf("computing pd_code (random rotation enabled)...");
  
  set_pd_code_from_plCurve_debug(false);
  projected_pd = pd_code_from_plCurve(rng,tref);
  
  end=clock();
  cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

  printf("done (%2.2g sec)\n",cpu_time_used);

  printf("checking pd_ok...");
  
  if (pd_ok(projected_pd)) { 
    
    printf("pass\n");
    pd_code_free(&projected_pd);

  } else { 

    printf("FAIL\n");
    printf("pd code generated from plCurve was:\n");
    printf("-----------------------------------\n");
    pd_write(stdout,projected_pd);
    pd_code_free(&projected_pd);

    exit(1);
  }

  plc_free(tref);
  return true;

}

void set_pd_code_from_plCurve_verbose(bool val);

bool randomwalk_test(gsl_rng *rng,int nedges) {

  clock_t start,end;
  double cpu_time_used;

  printf("------------------------------------------------\n"
	 "random equilateral n-gon test\n"
	 "------------------------------------------------\n");

  start = clock();
  printf("generating equilateral %d gon...",nedges);
  plCurve *L = plc_random_equilateral_closed_polygon(rng,nedges);
  end = clock();
  cpu_time_used = ((double)(end - start)/CLOCKS_PER_SEC);
  printf("done (%2.4g sec)\n",cpu_time_used);

  printf("computing pd_code (random rotation enabled, verbose enabled)...\n\n");
  set_pd_code_from_plCurve_verbose(true);
  pd_code_t *projected_pd = pd_code_from_plCurve(rng,L);
  printf("done\n");

  printf("crossings of produced pd code...%d\n",projected_pd->ncross);
  printf("checking pd_ok...");
  if (pd_ok(projected_pd)) { 
    
    printf("pass\n");
    pd_code_free(&projected_pd);

  } else { 

    printf("FAIL\n");
    printf("pd code generated from plCurve was:\n");
    printf("-----------------------------------\n");
    pd_write(stdout,projected_pd);
    pd_code_free(&projected_pd);

    exit(1);
  }

  plc_free(L);
  return true;
 
}
  

plCurve *equilateral_ngon(int nedges) {

  int n = nedges;
  double ngon_radius;
  double TWO_PI = 6.2831853071795864769;
  double phi,theta;
  int i;

  plCurve *L;
  bool open = false;
  int  nv = (int)(n);
  int  cc = 0;

  L = plc_new(1,&nv,&open,&cc);

  /* For the regular n-gon with side length 1, we have that the radius */
  /* of the circumscribed circle obeys:

     (1/2) / r = sin (theta/2), 

     or r = (1/2) / sin (theta/2), where theta = TWO_PI/n_edges. 
  */

  theta = TWO_PI/(double)(nv);
  ngon_radius = (1.0/2.0) / sin(theta/2.0);

  for(i=0,phi=0;i<L->cp[0].nv;i++,phi += theta) { 

    L->cp[0].vt[i].c[0] = ngon_radius*cos(phi);
    L->cp[0].vt[i].c[1] = ngon_radius*sin(phi);
    L->cp[0].vt[i].c[2] = 0;

  }

  plc_fix_wrap(L);
  return L;
}

bool unknot_and_split_component_test(gsl_rng *rng) {

  clock_t start,end;
  double cpu_time_used;

  printf("------------------------------------------------\n"
	 "unknot and split component test\n"
	 "------------------------------------------------\n");
  
  int nedges = 50;

  start = clock();
  printf("generating regular planar %d gon (0-crossing unknot)...",nedges);
  plCurve *L = equilateral_ngon(nedges);
  end = clock();
  cpu_time_used = ((double)(end - start)/CLOCKS_PER_SEC);
  printf("done (%2.4g sec)\n",cpu_time_used);

  printf("computing pd_code (random rotation enabled, verbose enabled)...\n\n");
  set_pd_code_from_plCurve_verbose(true);
  pd_code_t *projected_pd = pd_code_from_plCurve(rng,L);
  printf("done\n");

  printf("crossings of produced pd code == 1 (expect 1 virtual crossing)...");

  if (projected_pd->ncross == 1) { 
    
    printf("ok (%d crossings actual == 1)\n",projected_pd->ncross);

  } else {

    printf("FAIL (%d != 1)\n",projected_pd->ncross);
    printf("pd code generated from plCurve was:\n");
    printf("-----------------------------------\n");
    pd_write(stdout,projected_pd);
    pd_code_free(&projected_pd);

    exit(1);
  }

  printf("checking pd_ok...");
  if (pd_ok(projected_pd)) { 
    
    printf("pass\n");
    pd_code_free(&projected_pd);

  } else { 

    printf("FAIL\n");
    printf("pd code generated from plCurve was:\n");
    printf("-----------------------------------\n");
    pd_write(stdout,projected_pd);
    pd_code_free(&projected_pd);

    exit(1);
  }

  plc_free(L);

  /***********************/ 
  
 /*  start = clock(); */
/*   printf("generating (2,4) torus link + split component gon..."); */
/*   plCurve *split = equilateral_ngon(nedges); */
/*   plCurve *link = torusknot(250,4,2,5.0,2.0); */
/*   plc_add_component(link,2,split->cp[0].nv,split->cp[0].open,split->cp[0].cc,split->cp[0].vt,split->cp[0].clr); */
/*   end = clock(); */
/*   cpu_time_used = ((double)(end - start)/CLOCKS_PER_SEC); */
/*   printf("done (%2.4g sec)\n",cpu_time_used); */

/*   printf("writing to file splitlink_test.vect...."); */
/*   FILE *outfile = fopen("/Users/cantarel/plcurve/test/splitlink_test.vect","w"); */
/*   plc_write(outfile,link); */
/*   fclose(outfile); */
/*   printf("done\n"); */

/*   printf("computing pd_code (random rotation disabled, verbose enabled)...\n\n"); */

/*   set_pd_code_from_plCurve_verbose(true); */
/*   set_pd_code_from_plCurve_debug(true); */

/*   projected_pd = pd_code_from_plCurve(rng,link); */
/*   printf("done\n"); */

/*   printf("crossings of produced pd code == 5 (expect 4 link + 1 virtual crossing)..."); */

/*   if (projected_pd->ncross == 5) {  */
    
/*     printf("ok (%d crossings actual == 5)\n",projected_pd->ncross); */

/*   } else { */

/*     printf("FAIL (%d != 5)\n",projected_pd->ncross); */
/*     printf("pd code generated from plCurve was:\n"); */
/*     printf("-----------------------------------\n"); */
/*     pd_write(stdout,projected_pd); */
/*     free(projected_pd); */

/*     exit(1); */
/*   } */

/*   printf("checking pd_ok..."); */
/*   if (pd_ok(projected_pd)) {  */
    
/*     printf("pass\n"); */
/*     free(projected_pd); */

/*   } else {  */

/*     printf("FAIL\n"); */
/*     printf("pd code generated from plCurve was:\n"); */
/*     printf("-----------------------------------\n"); */
/*     pd_write(stdout,projected_pd); */
/*     free(projected_pd); */

/*     exit(1); */
/*   } */

/*   plc_free(split); */
/*   plc_free(link); */

  printf("disabling verbose mode...");
  set_pd_code_from_plCurve_verbose(false);
  printf("done\n");
  printf("---------------------------------------------\n");
  return true;
 
}


int main () {

  const gsl_rng_type * T;
  gsl_rng *rng;
   
  gsl_rng_env_setup();
  
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  
  int seedi;
  
  seedi = time(0); 
  gsl_rng_set(rng,seedi);
 
  printf("pd_code_from_plCurve test suite\n");   
  printf("with %s random number gen, seeded with %d.\n",gsl_rng_name(rng),seedi);
  printf("==========================================\n");

  torus_knot_test(rng,150,3);
  torus_knot_test(rng,150,4);
  torus_knot_test(rng,550,8);
  
  arcpresentation_tests(rng);
  unknot_and_split_component_test(rng);
  
  //torus_knot_test(rng,250,3);
  //torus_knot_test(rng,550,3);
  //torus_knot_test(rng,1050,3);
  //torus_knot_test(rng,150,7);
  //torus_knot_test(rng,550,9);

  torus_knot_rotation_test(rng,150,3);
  torus_knot_rotation_test(rng,150,4);

  randomwalk_test(rng,10);
  //randomwalk_test(rng,101);
  //randomwalk_test(rng,250);

  printf("=======================================\n"
	 "pd_code_from_plCurve test suite PASSED.\n");

  gsl_rng_free(rng);
  exit(0);

}
