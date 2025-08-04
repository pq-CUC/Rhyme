#include <stdio.h> 
#include <stdint.h>
#include "params.h"
#include "ntt.h"
#include "reduce.h"


int16_t zetas[128] = {
       1,   -16,    64,     4,    -8,   128,     2,   -32,
    -121,  -120,   -34,    30,   -60,   -68,    15,    17,
      81,   -11,    44,    67,   123,    88,   -95,   -22,
     -35,    46,    73,   117,    23,  -111,   -70,    92,
       9,   113,    62,    36,   -72,   124,    18,   -31,
     -61,   -52,   -49,    13,   -26,   -98,  -122,  -104,
     -42,   -99,  -118,    89,    79,    21,   -84,    59,
     -58,  -100,  -114,    25,   -50,    29,  -116,    57,
       3,   -48,   -65,    12,   -24,   127,     6,   -96,
    -106,  -103,  -102,    90,    77,    53,    45,    51,
     -14,   -33,  -125,   -56,   112,     7,   -28,   -66,
    -105,  -119,   -38,    94,    69,   -76,    47,    19,
      27,    82,   -71,   108,    41,   115,    54,   -93,
      74,   101,   110,    39,   -78,   -37,  -109,   -55,
    -126,   -40,   -97,    10,   -20,    63,     5,   -80,
      83,   -43,   -85,    75,   107,    87,   -91,   -86,
};
/*************************************************
* Name:        fqmul
*
* Description: Multiplication followed by Montgomery reduction
*
* Arguments:   - int16_t a: first factor
*              - int16_t b: second factor
*
* Returns 16-bit integer congruent to a*b*R^{-1} mod q
**************************************************/
int16_t fqmul(int16_t a, int16_t b) {
  return (int16_t)montgomery_reduce((int64_t)a * b);
}

/*************************************************
* Name:        ntt
*
* Description: Inplace number-theoretic transform (NTT) in Rq.
*              input is in standard order, output is in bitreversed order
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements of Zq
**************************************************/
void ntt(int16_t r[256]) {
  unsigned int len, start, j, k;
  int16_t t, zeta;

  k = 1;
  for(len = 128; len >= 2; len >>= 1) {
    for(start = 0; start < 256; start = j + len) {
      zeta = zetas[k++];
      for(j = start; j < start + len; j++) {
        t = fqmul(zeta, r[j + len]);
        r[j + len] = r[j] - t;
        r[j] = r[j] + t;
      }
    }
  }
}

/*************************************************
* Name:        invntt_tomont
*
* Description: Inplace inverse number-theoretic transform in Rq and
*              multiplication by Montgomery factor 2^16.
*              Input is in bitreversed order, output is in standard order
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements of Zq
**************************************************/
void invntt_tomont(int16_t r[256]) {
  unsigned int start, len, j, k;
  int16_t t, zeta;
  const int16_t f = 255; // mont^2/128

  k = 127;
  for(len = 2; len <= 128; len <<= 1) {
    for(start = 0; start < 256; start = j + len) {
      zeta = zetas[k--];
      for(j = start; j < start + len; j++) {
        t = r[j];
        r[j] = barrett_reduce(t + r[j + len]);
        r[j + len] = r[j + len] - t;
        r[j + len] = fqmul(zeta, r[j + len]);
      }
    }
  }

  for(j = 0; j < 256; j++)
    r[j] = fqmul(r[j], f);
}

/*************************************************
* Name:        basemul
*
* Description: Multiplication of polynomials in Zq[X]/(X^2-zeta)
* used for multiplication of elements in Rq in NTT domain.
*
* Arguments:   - int16_t r[2]: pointer to the output polynomial
* - const int16_t a[2]: pointer to the first factor
* - const int16_t b[2]: pointer to the second factor
* - int16_t zeta: integer defining the reduction polynomial
**************************************************/
void basemul(int16_t r[2], const int16_t a[2], const int16_t b[2], int16_t zeta)
{
  r[0]  = fqmul(a[1], b[1]);
  r[0]  = fqmul(r[0], zeta);
  r[0] += fqmul(a[0], b[0]);
  r[1]  = fqmul(a[0], b[1]);
  r[1] += fqmul(a[1], b[0]);
}


