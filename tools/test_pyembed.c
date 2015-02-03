#include <plCurve.h>
#include <plcTopology.h>
#include <pdcode/diagram_api.h>
#include <Python.h>
#include <stdio.h>

int main(int argc, char** argv) {
    pd_code_t *K, *L;
    pd_code_t **results;
    pd_code_t ***neighbors;
    int nn = 0;
    int *nd;
    int ndias = 0;
    int i =0;
    int j =0;
    FILE *f;
    //f = fopen("r2_testG_before.pdstor", "r");
    //L = pd_read(f);
    //fclose(f);

    Py_Initialize();
    import_libpl__pdcode__diagram();
    for (j=3; j<8; j++) {

        K = pd_build_unknot(j);
        int nd = 0;

        printf("Testing pd_simplify on %d-crossing unknot...\n", K->ncross);
        printf("Original diagram had %d crossings\n", K->ncross);

        results = pd_simplify(K, &nd);
        for (i=0; i<nd; i++) {
            printf("Result diagram %d has %d crossings\n",
                   i+1, results[i]->ncross);
            pd_code_free(&results[i]);
        }
        free(results);
        pd_code_free(&K);
        //pd_code_free(&L);
    }
    Py_Finalize();

    /* printf("\nTesting pd_neighbors on 7-crossing unknot...\n"); */
    /* neighbors = pd_neighbors(K, &nn, &nd); */
    /* for (i=0; i<nn; i++) { */
    /*     printf("Considering neighbor number %d\n",i+1); */
    /*     for (j=0; j<nd[i]; j++) { */
    /*         printf("Piece %d has %d crossings\n", */
    /*                j+1, neighbors[i][j]->ncross); */

    /*         assert(!strcmp(pd_homfly(K), pd_homfly(neighbors[i][j]))); */
    /*         pd_code_free(&neighbors[i][j]); */
    /*     } */
    /*     free(neighbors[i]); */
    /* } */
    /* free(neighbors); */
    /* free(nd); */

    /* printf("\nTesting pd_simplify on a R2 candidate which splits into 2...\n"); */
    /* printf("Original diagram had %d crossings\n", L->ncross); */

    /* results = pd_simplify(L, &ndias); */
    /* for (i=0; i<ndias; i++) { */
    /*     printf("Result diagram %d has %d crossings\n", */
    /*            i+1, results[i]->ncross); */
    /*     pd_code_free(&results[i]); */
    /* } */
    /* free(results); */

    /* printf("\nTesting pd_neighbors on R2 candidate...\n"); */
    /* neighbors = pd_neighbors(L, &nn, &nd); */
    /* for (i=0; i<nn; i++) { */
    /*     printf("Considering neighbor number %d\n",i+1); */
    /*     for (j=0; j<nd[i]; j++) { */
    /*         printf("Piece %d has %d crossings\n", */
    /*                j+1, neighbors[i][j]->ncross); */
    /*         pd_code_free(&neighbors[i][j]); */
    /*     } */
    /*     free(neighbors[i]); */
    /* } */
    /* free(neighbors); */
    /* free(nd); */


    return 0;
}
