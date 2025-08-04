#ifndef NTT_H
#define NTT_H

#include "params.h" 
#include "config.h" 
#include <stdint.h>

#define ntt ccc_NAMESPACE(ntt)
void ntt(int16_t a[N]);

#define invntt_tomont ccc_NAMESPACE(invntt_tomont)
void invntt_tomont(int16_t a[N]);

#define fqmul ccc_NAMESPACE(fqmul)
int16_t fqmul(int16_t a, int16_t b);

#define zetas ccc_NAMESPACE(zetas)
extern  int16_t zetas[128];

#define basemul ccc_NAMESPACE(basemul)
void basemul(int16_t r[2], const int16_t a[2], const int16_t b[2], int16_t zeta);

#endif