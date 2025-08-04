#ifndef CCC_SAMPLER_H
#define CCC_SAMPLER_H

#include "params.h"
#include "poly.h"
#include "polyvec.h" // For polyvecl type
#include <stdint.h>
#include "fips202.h"   // For stream state if needed directly, or via symmetric.h
#include "symmetric.h" // For stream states and functions

// Samples polynomial with coefficients {-1, 1}
void SampleSigma1(poly *p, const uint8_t seed[CRHBYTES], uint16_t nonce);

void SampleX_poly(poly *x, const poly *c, const uint8_t seed[CRHBYTES], uint16_t nonce);

void SampleY1(polyvecl *y1, const uint8_t seed[CRHBYTES], uint16_t nonce_base);

void SampleSPrime(polyvecl *s_prime, const uint8_t seed[CRHBYTES], uint16_t nonce_base);

void SampleY(polyvecl *y,
             const polyvecl *s_prime,
             const poly *s_bar_prime_0,         
             const uint8_t seed_y1[CRHBYTES],
             const uint8_t seed_c_prime[CRHBYTES],
             const uint8_t seed_x_prime[CRHBYTES],
             const uint8_t seed_c_prime2[CRHBYTES], 
             const uint8_t seed_x_prime2[CRHBYTES], 
             uint16_t nonce);

void SampleChallenge(poly *c, const uint8_t seed[CRHBYTES]);

unsigned int rej_uniform(int16_t *a, unsigned int len, const uint8_t *buf, unsigned int buflen, unsigned int *pos);

void poly_uniform_ntt(poly *a, const uint8_t seed[SEEDBYTES], uint16_t nonce);
#endif // CCC_SAMPLER_H