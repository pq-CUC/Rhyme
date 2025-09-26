#ifndef ccc_SIGN_H
#define ccc_SIGN_H

#include "params.h"
#include "poly.h"
#include "polymat.h"
#include "polyvec.h"
#include <stddef.h>
#include <stdint.h>

#define crypto_sign_keypair ccc_NAMESPACE(keypair)
int crypto_sign_keypair(uint8_t *pk, size_t *pklen, uint8_t *sk, size_t *sklen);

#define crypto_sign_signature ccc_NAMESPACE(signature)
int crypto_sign_signature(uint8_t *sig, size_t *siglen, const uint8_t *m,
                          size_t mlen, const uint8_t *sk);

#define crypto_sign_sign ccc_NAMESPACE(sign)
int crypto_sign_sign(uint8_t *sm, size_t *smlen, const uint8_t *m, size_t mlen,
                     const uint8_t *sk);

#define crypto_sign_verify ccc_NAMESPACE(verify)
int crypto_sign_verify(const uint8_t *sig, size_t siglen, const uint8_t *m, size_t mlen, 
                       const uint8_t *pk, size_t pklen);

#define crypto_sign_open ccc_NAMESPACE(open)
int crypto_sign_open(uint8_t *m, size_t *mlen, const uint8_t *sm, size_t smlen,
                     const uint8_t *pk, const size_t pklen);

#endif // ccc_SIGN_H