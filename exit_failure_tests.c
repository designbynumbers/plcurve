#include <config.h>
#include <plCurve.h>

#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif

int main(int argc,char **argv) {

  if (argc == 1) {
    switch (argv[0][0]) {
      case '1' :
        (void)plcl_normalize_vect(plcl_build_vect(0.0,0.0,0.0),NULL);
        break;
      case '2' :
        (void)plcl_component_div(plcl_build_vect(1.0,1.0,1.0),
                                 plcl_build_vect(1.0,1.0,0.0),NULL);
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
