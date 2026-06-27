#include <stdio.h>
#include <string.h>
#include "params.h"
#include "sign.h"
#include "packing.h"
#include "randombytes.h"

/* Construction 4 ("cut-F") sign/verify smoke test with negatives.
 * Uses the real keygen (the previous hardcoded test_basis_256.h was for the
 * old square-B geometry and no longer matches the cut-F secret_basis). */
int main(void) {
    static uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    static uint8_t sig[CRYPTO_BYTES];
    size_t siglen;

    if (crypto_sign_keypair(pk, sk)) { printf("keypair fail\n"); return 1; }
    printf("pk %d bytes, sk %d bytes, sigmax %d bytes\n",
           CRYPTO_PUBLICKEYBYTES, CRYPTO_SECRETKEYBYTES, CRYPTO_BYTES);

    uint8_t msg[59];
    for (int t = 0; t < 20; t++) {
        for (int i = 0; i < (int)sizeof msg; i++) msg[i] = (uint8_t)(t + i);
        if (crypto_sign_signature(sig, &siglen, msg, sizeof msg, sk)) {
            printf("sign fail at %d\n", t); return 1;
        }
        if (crypto_sign_verify(sig, siglen, msg, sizeof msg, pk)) {
            printf("verify fail at %d\n", t); return 1;
        }
        /* negative: flip message */
        msg[0] ^= 1;
        if (!crypto_sign_verify(sig, siglen, msg, sizeof msg, pk)) {
            printf("forged msg accepted at %d!\n", t); return 1;
        }
        msg[0] ^= 1;
        /* negative: corrupt sig (challenge bytes start at offset 2) */
        sig[2 + 5] ^= 0x40;
        if (!crypto_sign_verify(sig, siglen, msg, sizeof msg, pk)) {
            printf("corrupt sig accepted at %d!\n", t); return 1;
        }
        sig[2 + 5] ^= 0x40;
    }
    printf("sign/verify OK (20 messages, with negatives)\n");
    return 0;
}
