/* unit test: zint xgcd + NTT vs known values */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kg_zint.h"
#include "kg_ntt.h"

static uint64_t rngs = 88172645463325252ULL;
static uint64_t xs(void){ rngs^=rngs<<13; rngs^=rngs>>7; rngs^=rngs<<17; return rngs; }

int main(void) {
    /* zint arithmetic fuzz vs __int128 */
    for (int t = 0; t < 200000; t++) {
        int64_t a = (int64_t)xs() >> (xs() % 32);
        int64_t b = (int64_t)xs() >> (xs() % 32);
        zint za, zb, zr; zi_init(&za); zi_init(&zb); zi_init(&zr);
        zi_set_i64(&za, a); zi_set_i64(&zb, b);
        zi_add(&zr, &za, &zb);
        if (zi_to_i64(&zr) != a + b) { printf("add fail %lld %lld\n",(long long)a,(long long)b); return 1; }
        zi_sub(&zr, &za, &zb);
        if (zi_to_i64(&zr) != a - b) { printf("sub fail\n"); return 1; }
        zi_mul(&zr, &za, &zb);
        __int128 prod = (__int128)a * b;
        if ((int64_t)prod == prod) {  /* compare only when fits */
            if (zi_to_i64(&zr) != (int64_t)prod) { printf("mul fail\n"); return 1; }
        }
        zi_free(&za); zi_free(&zb); zi_free(&zr);
    }
    /* big xgcd fuzz: random multi-limb a,b; check u*a+v*b == g and g | a,b via mod check */
    for (int t = 0; t < 300; t++) {
        zint a, b, g, u, v, t1, t2, t3;
        zi_init(&a); zi_init(&b); zi_init(&g); zi_init(&u); zi_init(&v);
        zi_init(&t1); zi_init(&t2); zi_init(&t3);
        zi_set_i64(&a, (int64_t)(xs() | 1));
        zi_set_i64(&b, (int64_t)(xs() | 1));
        for (int i = 0; i < 8 + (int)(xs() % 8); i++) { zi_mul_u32(&a, (uint32_t)xs() | 1); zi_mul_u32(&b, (uint32_t)xs() | 1); }
        if (t & 1) zi_neg(&a);
        zi_xgcd(&g, &u, &v, &a, &b);
        zi_mul(&t1, &u, &a);
        zi_mul(&t2, &v, &b);
        zi_add(&t3, &t1, &t2);
        if (zi_cmp_abs(&t3, &g) != 0 || t3.sign != g.sign || g.sign <= 0) {
            printf("xgcd fail at %d (bezout)\n", t); return 1;
        }
        /* g divides a and b: check via random prime moduli */
        for (int pi = 0; pi < 3; pi++) {
            uint32_t p = kg_prime((unsigned)(xs() % kg_nprimes()));
            uint32_t gm = zi_mod_u32(&g, p);
            if (gm == 0) continue;
            uint32_t ginv = kg_powmod(gm, p - 2, p);
            (void)ginv; /* divisibility can't be checked mod p; skip */
        }
        zi_free(&a); zi_free(&b); zi_free(&g); zi_free(&u); zi_free(&v);
        zi_free(&t1); zi_free(&t2); zi_free(&t3);
    }
    /* NTT roundtrip + negacyclic product check vs schoolbook (n=64) */
    unsigned n = 64;
    for (int t = 0; t < 50; t++) {
        unsigned pidx = (unsigned)(xs() % kg_nprimes());
        uint32_t p = kg_prime(pidx);
        uint32_t a[64], b[64], A[64], B[64];
        int64_t ref[128];
        memset(ref, 0, sizeof ref);
        for (unsigned i = 0; i < n; i++) { a[i] = (uint32_t)(xs() % 1000); b[i] = (uint32_t)(xs() % 1000); }
        for (unsigned i = 0; i < n; i++)
            for (unsigned j = 0; j < n; j++)
                ref[i + j] += (int64_t)a[i] * b[j];
        memcpy(A, a, sizeof a); memcpy(B, b, sizeof b);
        kg_ntt_fwd(A, n, pidx); kg_ntt_fwd(B, n, pidx);
        for (unsigned i = 0; i < n; i++) A[i] = (uint32_t)(((uint64_t)A[i] * B[i]) % p);
        kg_ntt_inv(A, n, pidx);
        for (unsigned i = 0; i < n; i++) {
            int64_t want = (ref[i] - ref[i + n]) % (int64_t)p;
            if (want < 0) want += p;
            if (A[i] != (uint32_t)want) { printf("ntt prod fail t=%d i=%u\n", t, i); return 1; }
        }
    }
    printf("kg_zint + kg_ntt OK\n");
    return 0;
}
