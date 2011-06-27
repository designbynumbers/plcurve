/* 

   symmetry_test.c : Test code for the symmetry functions in plCurve. 

*/

#include<plCurve.h>
#include<config.h>

#ifdef HAVE_MATH_H
 #include<math.h>
#endif

#ifdef HAVE_STDLIB_H
 #include<stdlib.h>
#endif

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

    nv[i] = p*q*(verts/(p*q));  /* We need to make sure that verts is a multiple of pq to have symmetry. */
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

    for(vt=0,theta=0,pofs=i*(2*pi/q),tstep = 2*pi/(cp*(double)(nv[i]));
	vt<nv[i];
	vt++,theta+=tstep) {

      plc_vector loc;

      pangle = pofs + p*theta;
      qangle = q*theta;

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

void torus_tests(int verts, int p, int q)

/* Tests for the existence of the Z/pZ rotational symmetry in a torus knot or link. */

{
  plc_vector Zaxis = {{0,0,1}}; 
  plCurve *L;
  L = torusknot(verts,p,q,2.0,1.0); 

  printf("\nGenerated (%d,%d) torus knot with %d verts. Testing for Z/%dZ rotational symmetry.\n",p,q,plc_num_verts(L),p);
  L->G = plc_rotation_group(L,Zaxis,p);
  
  if (L->G != NULL) {

    printf("Symmetry construction passed.\n");

  } else {

    printf("FAILED to build %d-fold rotational symmetry.\n",p);
    exit(1);

  }

 /*  printf("Now checking for Z/%dZ rotational symmetry (should NOT be able to construct it).\n",p+1); */

  plc_symmetry_group *check;
  check = plc_rotation_group(L,Zaxis,p+1);

  if (check != NULL) {

    printf("FAILED. Detected spurious %d-fold symmetry in (%d,%d) torus knot.\n",p+1,p,q);
    exit(1);

  } else {

    printf("Passed. Could not construct symmetry.\n");
    
  }

  plc_symmetry_group_free(&check);

  double rad = 0.01;
  printf("Now perturbing with radius %g.\n",rad);
  plc_perturb(L,rad);

  printf("Perturbed version has error %g.\n",plc_symmetry_group_check(L));
  printf("Symmetrizing perturbed link.\n");
  plc_symmetrize(L);

  printf("Symmetrized version has error %g.\n",plc_symmetry_group_check(L));

  printf("Trying to rebuild symmetry group for symmetrized version.\n");
  check = plc_rotation_group(L,Zaxis,p);
  
  if (check != NULL) {

    printf("Passed. Symmetrized version has symmetry.\n");

  } else {

    printf("FAILED. Symmetrized version does NOT have desired symmetry.\n");
    exit(1);

  }

  printf("Checking plc_symmetry_group_free.\n");
  plc_symmetry_group_free(&check);

  /* Now we need to test the variation symmetrization. */

  printf("Now building random variation of this curve.\n");

  plc_vector *testvariation;
  testvariation = calloc(plc_num_verts(L),sizeof(plc_vector));

  int numverts,i;
  for(i=0,numverts=plc_num_verts(L);i<numverts;i++) {
    testvariation[i] = plc_build_vect((6.0)*rand()/RAND_MAX - 3,(6.0)*rand()/RAND_MAX - 3,(6.0)*rand()/RAND_MAX - 3);
  }

  printf("Generated random variation with symmetry error %g.\n",plc_symmetry_group_variation_check(L,testvariation));
  printf("Now symmetrizing variation.\n");
  
  plc_symmetrize_variation(L,testvariation);
  printf("After symmetrizing symmetry error is %g.\n",plc_symmetry_group_variation_check(L,testvariation));
  
  free(testvariation);
  printf("%d vertex (%d,%d) torus knot/link test PASSED.\n\n",plc_num_verts(L),p,q);

  plc_free(L);

}
  
  
int main() {

  torus_tests(841,3,2);
  torus_tests(41,3,5);
  torus_tests(247,2,4);
  torus_tests(341,4,2);
  torus_tests(243,6,9);

  return 0;

}

