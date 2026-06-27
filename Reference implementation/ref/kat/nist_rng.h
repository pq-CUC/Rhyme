#ifndef NIST_RNG_H
#define NIST_RNG_H
#include <stdint.h>
#include <stddef.h>
void randombytes_init(unsigned char *entropy_input,
                      unsigned char *personalization_string,
                      int security_strength);
int randombytes(uint8_t *x, size_t xlen);
void nist_aes256_ecb(const uint8_t key[32], const uint8_t in[16], uint8_t out[16]);
#endif
