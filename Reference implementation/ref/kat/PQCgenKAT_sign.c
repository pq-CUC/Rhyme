/* Standard NIST PQC KAT generator for the Rhyme signature scheme.
 * Produces PQCsignKAT_<sklen>.req / .rsp.  KAT count defaults to 100
 * (NIST standard); set the RHYME_KAT_COUNT environment variable to
 * generate fewer vectors (useful for the slow keygen of high modes). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "api.h"
#include "nist_rng.h"

#define KAT_SUCCESS 0
#define KAT_FILE_OPEN_ERROR -1
#define KAT_CRYPTO_FAILURE -4

static void fprintBstr(FILE *fp, const char *S, const unsigned char *A, size_t LEN) {
    fprintf(fp, "%s", S);
    if (LEN == 0) fprintf(fp, "00");
    for (size_t i = 0; i < LEN; i++) fprintf(fp, "%02X", A[i]);
    fprintf(fp, "\n");
}

int main(void) {
    char fn_req[64], fn_rsp[64];
    FILE *fp_req, *fp_rsp;
    static unsigned char seed[100][48];
    static unsigned char msgs[100][3300];
    unsigned char entropy_input[48];
    unsigned char *m, *sm, *m1;
    size_t mlen, smlen, mlen1;
    int count = 100;
    static unsigned char pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];

    const char *env = getenv("RHYME_KAT_COUNT");
    if (env) {
        int c = atoi(env);
        if (c > 0 && c <= 100) count = c;
    }

    sprintf(fn_req, "PQCsignKAT_%s_%d.req", CRYPTO_ALGNAME, CRYPTO_SECRETKEYBYTES);
    sprintf(fn_rsp, "PQCsignKAT_%s_%d.rsp", CRYPTO_ALGNAME, CRYPTO_SECRETKEYBYTES);
    if (!(fp_req = fopen(fn_req, "w"))) return KAT_FILE_OPEN_ERROR;
    if (!(fp_rsp = fopen(fn_rsp, "w"))) return KAT_FILE_OPEN_ERROR;

    for (int i = 0; i < 48; i++) entropy_input[i] = (unsigned char)i;
    randombytes_init(entropy_input, NULL, 256);
    for (int i = 0; i < count; i++) {
        fprintf(fp_req, "count = %d\n", i);
        randombytes(seed[i], 48);
        fprintBstr(fp_req, "seed = ", seed[i], 48);
        mlen = 33 * (size_t)(i + 1);
        fprintf(fp_req, "mlen = %zu\n", mlen);
        randombytes(msgs[i], mlen);
        fprintBstr(fp_req, "msg = ", msgs[i], mlen);
        fprintf(fp_req, "pk =\nsk =\nsmlen =\nsm =\n\n");
    }
    fclose(fp_req);

    fprintf(fp_rsp, "# %s\n\n", CRYPTO_ALGNAME);
    for (int i = 0; i < count; i++) {
        mlen = 33 * (size_t)(i + 1);
        m = malloc(mlen);
        m1 = malloc(mlen + CRYPTO_BYTES);
        sm = malloc(mlen + CRYPTO_BYTES);
        memcpy(m, msgs[i], mlen);            /* per NIST: msg comes from the req */
        randombytes_init(seed[i], NULL, 256);
        fprintf(fp_rsp, "count = %d\n", i);
        fprintBstr(fp_rsp, "seed = ", seed[i], 48);
        fprintf(fp_rsp, "mlen = %zu\n", mlen);
        fprintBstr(fp_rsp, "msg = ", m, mlen);
        if (crypto_sign_keypair(pk, sk)) { printf("keypair failure\n"); return KAT_CRYPTO_FAILURE; }
        fprintBstr(fp_rsp, "pk = ", pk, CRYPTO_PUBLICKEYBYTES);
        fprintBstr(fp_rsp, "sk = ", sk, CRYPTO_SECRETKEYBYTES);
        if (crypto_sign(sm, &smlen, m, mlen, sk)) { printf("sign failure\n"); return KAT_CRYPTO_FAILURE; }
        fprintf(fp_rsp, "smlen = %zu\n", smlen);
        fprintBstr(fp_rsp, "sm = ", sm, smlen);
        fprintf(fp_rsp, "\n");
        if (crypto_sign_open(m1, &mlen1, sm, smlen, pk)) { printf("open failure\n"); return KAT_CRYPTO_FAILURE; }
        if (mlen != mlen1 || memcmp(m, m1, mlen)) { printf("message mismatch\n"); return KAT_CRYPTO_FAILURE; }
        free(m); free(m1); free(sm);
        printf("kat %d/%d ok\n", i + 1, count);
        fflush(stdout);
    }
    fclose(fp_rsp);
    printf("KAT files written: %s %s\n", fn_req, fn_rsp);
    return KAT_SUCCESS;
}
