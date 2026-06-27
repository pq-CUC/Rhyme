#include <stdio.h>
#include <string.h>
#include "params.h"
#include "ntt.h"
#include "reduce.h"
#include "ntt_vectors.h"

#if RHYME_MODE == 128
#define TV(x) tv_##x##_3329_256
#elif RHYME_MODE == 256
#define TV(x) tv_##x##_9473_512
#elif RHYME_MODE == 384
#define TV(x) tv_##x##_11777_512
#else
#define TV(x) tv_##x##_18433_1024
#endif

int main(void) {
    int32_t a[N], b[N], c[N];
    memcpy(a, TV(a), sizeof a);
    memcpy(b, TV(b), sizeof b);
    /* End-to-end check against schoolbook negacyclic ground truth:
     * ntt -> basemul -> invntt_tomont must reproduce a*b mod (x^n+1, q). */
    ntt(a);
    ntt(b);
    basemul(c, a, b);
    invntt_tomont(c);
    for (int i = 0; i < N; i++) {
        if (freeze(c[i]) != freeze(TV(prod)[i])) {
            printf("PROD mismatch at %d: got %d want %d\n", i, freeze(c[i]), TV(prod)[i]);
            return 1;
        }
    }
    printf("NTT OK (mode %d, q=%d, n=%d)\n", RHYME_MODE, Q, N);
    return 0;
}
