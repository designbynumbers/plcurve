#include <plCurve.h>
#include <plcTopology.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>

int PD_VERBOSE=0;

int main(int argc, char** argv) {
    plCurve* K;
    gsl_rng *r;
    const gsl_rng_type *T;
    //    int i =0;
    int j =0;
    //double theta;

    //f = fopen("r2_testG_before.pdstor", "r");
    //L = pd_read(f);
    //fclose(f);

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    for (j=1; j<500000; j++) {
        //theta = gsl_rng_uniform(r);
        //K = plc_random_closed_polygon(r, 50);
        K = plc_random_equilateral_closed_polygon(r, 50);

        plc_free(K);
    }

    return 0;
}
