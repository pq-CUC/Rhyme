#ifndef CCC_PACKING_H
#define CCC_PACKING_H

#include "params.h"
#include "polyvec.h"
#include "poly.h"
#include <stdint.h>
#include <stddef.h> 
#include "config.h"
#include "encoding.h" // Include encoding header

// Public Key Packing/Unpacking
#define pack_pk ccc_NAMESPACE(pack_pk)
size_t pack_pk(uint8_t *pk, const uint8_t seedA[SEEDBYTES], const polyveck *b);
#define unpack_pk ccc_NAMESPACE(unpack_pk)
int unpack_pk(uint8_t seedA[SEEDBYTES], polyveck *b, const uint8_t *pk, size_t pklen);
 
// Secret Key Packing/Unpacking
#define pack_sk ccc_NAMESPACE(pack_sk)
size_t pack_sk(uint8_t *sk, const uint8_t *pk, size_t pklen, const polyvecm *s_gen, const polyveck *e_gen, const polyvecl *s_prime, const poly *s_bar_prime_0, const uint8_t key[SEEDBYTES]);

#define unpack_sk ccc_NAMESPACE(unpack_sk)
int unpack_sk(uint8_t *pk, size_t *pklen, polyvecm *s_gen, polyveck *e_gen, polyvecl *s_prime, poly *s_bar_prime_0, uint8_t key[SEEDBYTES], const uint8_t *sk, size_t sklen);
// Signature Packing/Unpacking
#define pack_sig ccc_NAMESPACE(pack_sig)

size_t pack_sig(uint8_t sig[CRYPTO_BYTES], const polyvecl *z, const poly *c);

#define unpack_sig ccc_NAMESPACE(unpack_sig)
int unpack_sig(polyvecl *z, poly *c, const uint8_t *sig, size_t siglen);



#define pack_polyveck_sigma1 ccc_NAMESPACE(pack_polyveck_sigma1)
void pack_polyveck_sigma1(uint8_t *buf, const polyveck *s_k);
#define unpack_polyveck_sigma1 ccc_NAMESPACE(unpack_polyveck_sigma1)
void unpack_polyveck_sigma1(polyveck *s_k, const uint8_t *buf);

#define pack_polyvecm_sigma1 ccc_NAMESPACE(pack_polyvecm_sigma1)
void pack_polyvecm_sigma1(uint8_t *buf, const polyvecm *s_m);
#define unpack_polyvecm_sigma1 ccc_NAMESPACE(unpack_polyvecm_sigma1)
void unpack_polyvecm_sigma1(polyvecm *s_m, const uint8_t *buf);

#define pack_polyvecl_sigma1 ccc_NAMESPACE(pack_polyvecl_sigma1)
void pack_polyvecl_sigma1(uint8_t *buf, const polyvecl *s);
#define unpack_polyvecl_sigma1 ccc_NAMESPACE(unpack_polyvecl_sigma1)
void unpack_polyvecl_sigma1(polyvecl *s, const uint8_t *buf);

#define pack_poly_challenge ccc_NAMESPACE(pack_poly_challenge)
void pack_poly_challenge(uint8_t *buf, const poly *c);
#define unpack_poly_challenge ccc_NAMESPACE(unpack_poly_challenge)
void unpack_poly_challenge(poly *c, const uint8_t *buf);

#define pack_polyveck_2q ccc_NAMESPACE(pack_polyveck_2q)
void pack_polyveck_2q(uint8_t *buf, const polyveck *v);


#define pack_polyvecl_s_prime ccc_NAMESPACE(pack_polyvecl_s_prime)
void pack_polyvecl_s_prime(uint8_t *buf, const polyvecl *s_prime);

#define unpack_polyvecl_s_prime ccc_NAMESPACE(unpack_polyvecl_s_prime)
void unpack_polyvecl_s_prime(polyvecl *s_prime, const uint8_t *buf);
#endif // CCC_PACKING_H