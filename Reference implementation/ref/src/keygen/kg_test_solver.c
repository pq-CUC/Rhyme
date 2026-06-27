#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "kg_solver.h"

static uint64_t rngs = 1;
static uint64_t xs(void){ rngs^=rngs<<13; rngs^=rngs>>7; rngs^=rngs<<17; return rngs; }
static int32_t cbd(unsigned eta){
    int32_t a=0,b=0;
    for(unsigned i=0;i<eta;i++){ uint64_t r=xs(); a+=r&1; b+=(r>>1)&1; }
    return a-b;
}

int main(int argc, char **argv) {
    unsigned n = argc > 1 ? (unsigned)atoi(argv[1]) : 64;
    unsigned d = argc > 2 ? (unsigned)atoi(argv[2]) : 3;
    unsigned eta = argc > 3 ? (unsigned)atoi(argv[3]) : 2;
    unsigned seed = argc > 4 ? (unsigned)atoi(argv[4]) : 7;
    rngs = seed * 2654435761u + 1;

    int32_t *cols = malloc((size_t)(d-1)*d*n*4);
    int32_t *F = malloc((size_t)d*n*4);
    clock_t t0 = clock();
    int tries = 0, rc;
    do {
        for (unsigned j = 0; j < d-1; j++)
            for (unsigned i = 0; i < d; i++)
                for (unsigned t = 0; t < n; t++)
                    cols[((size_t)j*d+i)*n+t] = cbd(eta);
        rc = kg_solve_last_column(F, cols, n, d);
        tries++;
        if (tries > 40) { printf("too many tries (last rc=%d)\n", rc); return 1; }
    } while (rc != 0);
    double el = (double)(clock()-t0)/CLOCKS_PER_SEC;
    if (kg_det_check(cols, F, n, d) != 0) { printf("DET CHECK FAILED\n"); return 1; }
    int mx = 0; double ss = 0;
    for (unsigned i = 0; i < d; i++)
        for (unsigned t = 0; t < n; t++) {
            int v = abs(F[(size_t)i*n+t]);
            if (v > mx) mx = v;
            ss += (double)F[(size_t)i*n+t]*F[(size_t)i*n+t];
        }
    printf("OK n=%u d=%u eta=%u: tries=%d det=1, max|F|=%d std=%.2f, %.1fs\n",
           n, d, eta, tries, mx, __builtin_sqrt(ss/(d*n)), el);
    return 0;
}
