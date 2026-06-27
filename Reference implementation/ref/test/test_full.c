#include <stdio.h>
#include <string.h>
#include <time.h>
#include "params.h"
#include "sign.h"

int main(void) {
    static uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    static uint8_t sig[CRYPTO_BYTES];
    size_t siglen;
    clock_t t0 = clock();
    if (crypto_sign_keypair(pk, sk)) { printf("keypair FAIL\n"); return 1; }
    double kg = (double)(clock()-t0)/CLOCKS_PER_SEC;
    uint8_t msg[33];
    t0 = clock();
    int nsig = 10;
    for (int t = 0; t < nsig; t++) {
        for (int i = 0; i < (int)sizeof msg; i++) msg[i] = (uint8_t)(7*t + i);
        if (crypto_sign_signature(sig, &siglen, msg, sizeof msg, sk)) { printf("sign FAIL %d\n", t); return 1; }
        if (crypto_sign_verify(sig, siglen, msg, sizeof msg, pk)) { printf("verify FAIL %d\n", t); return 1; }
        msg[3] ^= 2;
        if (!crypto_sign_verify(sig, siglen, msg, sizeof msg, pk)) { printf("forgery accepted %d\n", t); return 1; }
    }
    double sg = (double)(clock()-t0)/CLOCKS_PER_SEC;
    printf("MODE %d FULL OK: keygen %.1fs, %d sign+verify in %.2fs (pk=%d sk=%d sig=%zu actual, cap %d)\n",
           RHYME_MODE, kg, nsig, sg, CRYPTO_PUBLICKEYBYTES, CRYPTO_SECRETKEYBYTES, siglen, CRYPTO_BYTES);
    return 0;
}
