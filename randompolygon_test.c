/* 

   randompolygon_test.c : Test code for the symmetry functions in plCurve. 

*/

#include<plCurve.h>
#include<config.h>

#ifdef HAVE_MATH_H
 #include<math.h>
#endif

#ifdef HAVE_STDLIB_H
 #include<stdlib.h>
#endif

#ifdef HAVE_TIME_H
#include<time.h>
#endif

int main() {

  printf("Random Polygon Generation Tests \n"
	 "------------------------------- \n"
	 "plCurve generates random closed polygons using the direct sampling from the symmetric measure algorithm\n"
	 "of Cantarella, Deguchi, and Shonkwiler. The call plc_random_polygon(n) generates a polygon directly \n"
	 "sampled from this distribution. The polygon is guaranteed to be closed and have length 2.\n\n");

  plCurve *test;
  int i,j;
  int verts = 200;

  clock_t start,end;
  double cpu_time_used;

  printf("Verts        Length      Time to generate \n");

  for(i=0;i<4;i++,verts*=10) {

    for(j=0;j<5;j++) {

      start = clock();
      test = plc_random_polygon(verts);
      end = clock();
      cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;

      printf("%07d  %05g  %05g \n",verts,plc_arclength(test,NULL),cpu_time_used);

      plc_free(test);

    }

    printf("\n");

  }

}
