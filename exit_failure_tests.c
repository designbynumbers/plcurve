#include <config.h>
#include <plCurve.h>

#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif
#ifdef HAVE_STDIO_H
  #include <stdio.h>
#endif

/* 
 */

int main(int argc,char **argv) {

  if (argc == 2) {
    switch (argv[1][0]) {
      case '1' :
        (void)plc_normalize_vect(plc_build_vect(0.0,0.0,0.0),NULL);
        break;
      case '2' :
        (void)plc_component_div(plc_build_vect(1.0,1.0,1.0),
                                plc_build_vect(1.0,1.0,0.0),NULL);
        break;
      case '3' :
        exit(EXIT_SUCCESS);
    }
  }

  /* Because we are expecting EXIT_FAILURE, we exit with something completely
     different */
  /*@-exitarg@*/ /* So splint doesn't complain */
  exit(17);
  /*@=exitarg@*/
}
