#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "params.h"
#include "sign.h"
#include "poly.h"
#include "ntt.h"
#include "sampler.h"
#include "packing.h"
#include "symmetric.h"
#include "cpucycles.h"

#define NWARM 3
#define NSIGN 100
#define NPRIM 2000   /* iterations for per-primitive micro-benchmarks */

static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : x > y;
}
static uint64_t median(uint64_t *t, int n) {
    qsort(t, n, sizeof(uint64_t), cmp_u64);
    return t[n/2];
}
static uint64_t average(const uint64_t *t, int n) {
    uint64_t acc = 0;
    for (int i = 0; i < n; i++) acc += t[i];
    return acc / n;
}
/* print median + average for a sample array (matches the reference
 * speed_print format: two lines per operation). average() must run before
 * median() sorts the array — report() handles the order. */
static void report(const char *name, uint64_t *t, int n) {
    uint64_t avg = average(t, n);
    uint64_t med = median(t, n);          /* sorts t in place */
    printf("%s\n", name);
    printf("  median:  %12llu cycles\n", (unsigned long long)med);
    printf("  average: %12llu cycles\n", (unsigned long long)avg);
}

/* keep results from being optimised away */
static volatile uint64_t g_sink;

int main(void) {
    static uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    static uint8_t sig[CRYPTO_BYTES];
    static uint64_t tsign[NSIGN], tverify[NSIGN], tprim[NPRIM];
    size_t siglen;
    uint8_t msg[59] = {1,2,3};

    printf("=== %s  (pk %d B, sk %d B, sig cap %d B) ===\n",
           CRYPTO_ALGNAME, CRYPTO_PUBLICKEYBYTES, CRYPTO_SECRETKEYBYTES, CRYPTO_BYTES);

    /* ---------------- key generation ---------------- */
    clock_t c0 = clock();
    uint64_t k0 = cpucycles();
    if (crypto_sign_keypair(pk, sk)) { printf("keypair fail\n"); return 1; }
    uint64_t k1 = cpucycles();
    printf("\n-- full operations --\n");
    printf("keypair:        %10.2f s  (%llu cycles)\n",
           (double)(clock() - c0) / CLOCKS_PER_SEC, (unsigned long long)(k1 - k0));

    /* ---------------- full sign / verify ---------------- */
    size_t slen_sum = 0;
    for (int i = 0; i < NWARM; i++) {
        crypto_sign_signature(sig, &siglen, msg, sizeof msg, sk);
        crypto_sign_verify(sig, siglen, msg, sizeof msg, pk);
    }
    for (int i = 0; i < NSIGN; i++) {
        msg[0] = (uint8_t)i;
        uint64_t a = cpucycles();
        if (crypto_sign_signature(sig, &siglen, msg, sizeof msg, sk)) { printf("sign fail\n"); return 1; }
        uint64_t b = cpucycles();
        if (crypto_sign_verify(sig, siglen, msg, sizeof msg, pk)) { printf("verify fail\n"); return 1; }
        uint64_t c = cpucycles();
        tsign[i] = b - a; tverify[i] = c - b; slen_sum += siglen;
    }
    report("sign", tsign, NSIGN);
    report("verify", tverify, NSIGN);
    printf("avg sig size:   %10zu B\n", slen_sum / NSIGN);

    /* ---------------- per-primitive micro-benchmarks ---------------- */
    printf("\n-- core primitives (median + average over %d) --\n", NPRIM);

    /* poly NTT / iNTT */
    poly a; for (int j = 0; j < N; j++) a.coeffs[j] = (int32_t)(j % Q);
    for (int i = 0; i < NPRIM; i++) {
        poly t = a; uint64_t s = cpucycles();
        poly_ntt(&t); uint64_t e = cpucycles();
        tprim[i] = e - s; g_sink ^= t.coeffs[0];
    }
    report("poly_ntt", tprim, NPRIM);

    poly an = a; poly_ntt(&an);
    for (int i = 0; i < NPRIM; i++) {
        poly t = an; uint64_t s = cpucycles();
        poly_invntt_tomont(&t); uint64_t e = cpucycles();
        tprim[i] = e - s; g_sink ^= t.coeffs[0];
    }
    report("poly_invntt", tprim, NPRIM);

    /* pointwise montgomery */
    poly b2 = an;
    for (int i = 0; i < NPRIM; i++) {
        poly t; uint64_t s = cpucycles();
        poly_basemul(&t, &an, &b2); uint64_t e = cpucycles();
        tprim[i] = e - s; g_sink ^= t.coeffs[0];
    }
    report("poly_basemul", tprim, NPRIM);

    /* challenge sampling */
    uint8_t cseed[CTILDEBYTES] = {0};
    for (int i = 0; i < NPRIM; i++) {
        poly c; cseed[0] = (uint8_t)i;
        uint64_t s = cpucycles();
        SampleChallenge(&c, cseed); uint64_t e = cpucycles();
        tprim[i] = e - s; g_sink ^= c.coeffs[0];
    }
    report("SampleChallenge", tprim, NPRIM);

    /* y1 gaussian sampling (CDT) */
    uint8_t yseed[CRHBYTES] = {0};
    for (int i = 0; i < NPRIM; i++) {
        poly y; yseed[0] = (uint8_t)i;
        uint64_t s = cpucycles();
        SampleY1(&y, yseed, (uint16_t)i); uint64_t e = cpucycles();
        tprim[i] = e - s; g_sink ^= y.coeffs[0];
    }
    report("SampleY1(CDT)", tprim, NPRIM);

    /* SHAKE256 absorb+squeeze of one block (symmetric primitive) */
    uint8_t hbuf[CRHBYTES]; uint8_t hin[64] = {0};
    for (int i = 0; i < NPRIM; i++) {
        hin[0] = (uint8_t)i;
        uint64_t s = cpucycles();
        shake256(hbuf, CRHBYTES, hin, sizeof hin); uint64_t e = cpucycles();
        tprim[i] = e - s; g_sink ^= hbuf[0];
    }
    report("shake256(64B)", tprim, NPRIM);

    (void)g_sink;
    return 0;
}
