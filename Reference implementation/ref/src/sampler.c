#include "sampler.h"
#include "symmetric.h"
#include "params.h"
#include "reduce.h"
#include "polyvec.h"
#include "poly.h"
#include "fips202.h"
#include "aes256ctr.h"
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h> 

#ifndef STREAM256_BLOCKBYTES
#define STREAM256_BLOCKBYTES SHAKE256_RATE
#endif
#ifndef STREAM128_BLOCKBYTES
#define STREAM128_BLOCKBYTES SHAKE128_RATE
#endif

static uint8_t get_random_bit(stream256_state *state, uint8_t *buf, size_t *buf_len, size_t *buf_pos);
static void get_random_bytes(uint8_t *out, size_t num_bytes, stream256_state *state, uint8_t *buf, size_t *buf_len, size_t *buf_pos);


/*************************************************
* Name:        sample_chi_k
*
* Description: Samples a value from the difference of two Hamming weights
* of k random bits. Used for sampling noise.
*
* Arguments:   - int32_t k: the multiplier factor (and number of bits)
* - stream256_state *state: pointer to PRNG state
* - uint8_t *buf: pointer to random buffer
* - size_t *buf_len: pointer to current buffer length
* - size_t *buf_pos: pointer to current buffer position
*
* Returns:     Sampled integer value.
**************************************************/
static int32_t sample_chi_k(int32_t k, stream256_state *state, uint8_t *buf, size_t *buf_len, size_t *buf_pos) {
    uint8_t b1 = get_random_bit(state, buf, buf_len, buf_pos);
    uint8_t b2 = get_random_bit(state, buf, buf_len, buf_pos);
    return k * ((int32_t)b1 - (int32_t)b2);
}

/*************************************************
* Name:        hw_bytes
*
* Description: Computes the Hamming weight (population count) of a byte array.
*
* Arguments:   - const uint8_t *buf: pointer to the input buffer
* - size_t nbytes: number of bytes to count
*
* Returns:     Total Hamming weight as int32_t.
**************************************************/
static inline int32_t hw_bytes(const uint8_t *buf, size_t nbytes) {
    int32_t r = 0;
    size_t i = 0;
    uint64_t val64;

    // Processes 8 bytes (64 bits) at a time.
    while (i + 8 <= nbytes) {
        memcpy(&val64, buf + i, 8);
        r += __builtin_popcountll(val64);
        i += 8;
    }

    while (i < nbytes) {
        r += __builtin_popcount(buf[i]);
        i++;
    }
    return r;
}

/*************************************************
* Name:        SampleY0
*
* Description: Samples the first polynomial y0 of the masking vector y.
* It uses the difference of Hamming weights of two random bit streams.
*
* Arguments:   - poly *y: pointer to output polynomial
* - const uint8_t seed[]: pointer to input seed (CRHBYTES)
* - uint16_t nonce: nonce for the sampling process
**************************************************/
void SampleY0(poly *y, const uint8_t seed[CRHBYTES], uint16_t nonce) {

    const uint8_t *key_part1 = seed;       
    const uint8_t *key_part2 = seed + 32;  

    const size_t bytes_per_sign = TAU0 / 8;
    const size_t total_bytes = 2 * bytes_per_sign;
    const size_t nblocks = total_bytes / 64; 

    uint8_t buf[total_bytes + 64]; 
    aes256ctr_ctx state;

    // Only one initialization is performed
    aes256ctr_init(&state, key_part1, 0); 

    for (int i = 0; i < N / 2; ++i) {
        uint64_t combined_nonce = ((uint64_t)nonce << 16) | (uint16_t)i;
        uint64_t retry_ctr = 0;
        int32_t val;

        do {
            uint64_t final_nonce = combined_nonce ^ (retry_ctr << 48);
            aes256ctr_set_nonce(&state, final_nonce);

            aes256ctr_squeezeblocks(buf, nblocks, &state);

            int32_t w_pos = hw_bytes(buf, bytes_per_sign);
            int32_t w_neg = hw_bytes(buf + bytes_per_sign, bytes_per_sign);
            val = w_pos - w_neg;
            
            retry_ctr++;
        } while (abs(val) > B0 + 1);
        
        y->coeffs[i] = val;
    }
    aes256ctr_free(&state);
    aes256ctr_init(&state, key_part2, 0);

    for (int i = N / 2; i < N; ++i) {
        uint64_t combined_nonce = ((uint64_t)nonce << 16) | (uint16_t)i;
        uint64_t retry_ctr = 0;
        int32_t val;

        do {
            uint64_t final_nonce = combined_nonce ^ (retry_ctr << 48);
            aes256ctr_set_nonce(&state, final_nonce);
            
            aes256ctr_squeezeblocks(buf, nblocks, &state);

            int32_t w_pos = hw_bytes(buf, bytes_per_sign);
            int32_t w_neg = hw_bytes(buf + bytes_per_sign, bytes_per_sign);
            val = w_pos - w_neg;
            
            retry_ctr++;
        } while (abs(val) > B0 + 1);
        
        y->coeffs[i] = val;
    }
    aes256ctr_free(&state);
}
/*************************************************
* Name:        SampleY_rest
*
* Description: Samples the remaining part of the masking vector y (y_rest).
* Each coefficient is sampled from Chi(0, 2).
*
* Arguments:   - polyvecd_rest *y_rest: pointer to output polynomial vector
* - const uint8_t seed[]: pointer to input seed (CRHBYTES)
* - uint16_t nonce: nonce for the sampling process
**************************************************/
void SampleY_rest(polyvecd_rest *y_rest, const uint8_t seed[CRHBYTES], uint16_t nonce) {
    uint8_t buf[STREAM256_BLOCKBYTES];
    size_t buf_len = 0;
    size_t buf_pos = STREAM256_BLOCKBYTES * 8;
    stream256_state state;
    stream256_init(&state, seed, nonce);

    for (int i = 0; i < D_REST; ++i) {
        for (int j = 0; j < N; ++j) {
             y_rest->vec[i].coeffs[j] = sample_chi_k(2, &state, buf, &buf_len, &buf_pos);
        }
    }
}

/*************************************************
* Name:        get_random_double
*
* Description: Generates a random double in [0, 1) using 32 random bits.
* Optimized to read 4 bytes at once.
*
* Arguments:   - stream256_state *state: pointer to stream state
* - uint8_t *buf: pointer to buffer
* - size_t *buf_len: pointer to buffer length
* - size_t *buf_pos: pointer to current bit position
*
* Returns:     Random double value.
**************************************************/

// static double get_random_double(stream256_state *state, uint8_t *buf, size_t *buf_len, size_t *buf_pos) {
//     uint32_t r;

//     /* 1. Buffer check (byte-wise) */
//     if (*buf_pos + 32 > (*buf_len) * 8) {
//         stream256_squeezeblocks(buf, 1, state);
//         *buf_len = STREAM256_BLOCKBYTES;
//         *buf_pos = 0;
//     }

//     /* 2. Read 4 bytes at once */
//     size_t idx = *buf_pos / 8;
//     r = (uint32_t)buf[idx] 
//       | ((uint32_t)buf[idx + 1] << 8) 
//       | ((uint32_t)buf[idx + 2] << 16) 
//       | ((uint32_t)buf[idx + 3] << 24);

//     /* 3. Update pointer */
//     *buf_pos += 32;

//     /* 4. Return double */
//     return (double)r / 4294967296.0;
// }


/*************************************************
* Name:        get_random_uint32
*
* Description: Extracts a 32-bit unsigned integer from a buffered stream.
* Refills the buffer from the stream when necessary.
* Note: The buffer position (buf_pos) is tracked in bits.
*
* Arguments:   - stream256_state *state: pointer to the stream state
* - uint8_t *buf: pointer to the buffer
* - size_t *buf_len: pointer to the buffer length (in bytes)
* - size_t *buf_pos: pointer to the current bit position in the buffer
*
* Returns:     A 32-bit unsigned integer.
**************************************************/
static inline uint32_t get_random_uint32(stream256_state *state, uint8_t *buf, size_t *buf_len, size_t *buf_pos) {

    if (*buf_pos + 32 > (*buf_len) * 8) { 
        stream256_squeezeblocks(buf, 1, state);
        *buf_len = STREAM256_BLOCKBYTES;
        *buf_pos = 0;
    }
    
    uint32_t r;
    memcpy(&r, &buf[*buf_pos / 8], 4);
    
    *buf_pos += 32;
    return r; 
}


/*************************************************
* Name:        Sample_Z
*
* Description: Performs the 3C sampling rejection step for the first component.
* Determines z0 and x based on the challenge c and masking y0.
*
* Arguments:   - poly *z0: pointer to output polynomial z0
* - poly *x: pointer to output polynomial x
* - const poly *c: pointer to challenge polynomial
* - const poly *y0: pointer to masking polynomial
* - stream256_state *state: pointer to RNG state
*
* Returns:     1 if sampling is successful (accept), 0 otherwise (reject).
**************************************************/
int Sample_Z(poly *z0, poly *x, const poly *c, const poly *y0, stream256_state *state) {
    uint8_t buf[STREAM256_BLOCKBYTES];
    size_t buf_len = 0;
    size_t buf_pos = STREAM256_BLOCKBYTES * 8; 

   /* Pre-compute pointers for optimization */
    const int16_t *c_coeffs = c->coeffs;
    const int16_t *y0_coeffs = y0->coeffs;
    int16_t *z0_coeffs = z0->coeffs;
    int16_t *x_coeffs = x->coeffs;

    for (int i = 0; i < N; ++i) {
        int32_t y_val = y0_coeffs[i];
        
        /* Main loop */
        if (c_coeffs[i] == 1) {
            /* 1. Get random number */
            uint32_t rnd32 = get_random_uint32(state, buf, &buf_len, &buf_pos);
            
            /* 2. Integer arithmetic (no division) */
            int64_t num_a = TAU0 + y_val;
            int64_t den_a = TAU0 - y_val + 1;
            
            // Check 1
            uint64_t lhs_a = (uint64_t)rnd32 * den_a;
            uint64_t rhs_a = (uint64_t)SCALED_M_INV * num_a;
            
            if (lhs_a < rhs_a) {
                z0_coeffs[i] = y_val - 1;
                x_coeffs[i] = -2;
            } else {
                // Check 2
                int64_t num_b = TAU0 - y_val;
                int64_t den_b = TAU0 + y_val + 1;
                
                uint64_t lhs_total = lhs_a * den_b; 
                uint64_t rhs_total = (uint64_t)SCALED_M_INV * (num_a * den_b + num_b * den_a);
                
                if (lhs_total < rhs_total) {
                    z0_coeffs[i] = y_val + 1;
                    x_coeffs[i] = 0;
                } else {
                    return 0; // Reject
                }
            }
        } else {
            z0_coeffs[i] = y_val;
            x_coeffs[i] = 0;
        }
        
        /* Final check: Norm validation */
        if (abs(z0_coeffs[i]) > B0) return 0;
    }
    return 1;
}
/*************************************************
* Name:        get_random_bit
*
* Description: Extracts a single random bit from a buffered stream.
* Refills the buffer from the stream when necessary.
*
* Arguments:   - stream256_state *state: pointer to the stream state
* - uint8_t *buf: pointer to the buffer
* - size_t *buf_len: pointer to the buffer length
* - size_t *buf_pos: pointer to the current bit position in the buffer
*
* Returns:     A single random bit (0 or 1).
**************************************************/
static uint8_t get_random_bit(stream256_state *state, uint8_t *buf, size_t *buf_len, size_t *buf_pos) {
    if (*buf_pos >= (*buf_len) * 8) {
        stream256_squeezeblocks(buf, 1, state);
        *buf_len = STREAM256_BLOCKBYTES;
        *buf_pos = 0;
    }
    uint8_t bit = (buf[*buf_pos / 8] >> (*buf_pos % 8)) & 1;
    (*buf_pos)++;
    return bit;
}

/*************************************************
* Name:        get_random_bytes
*
* Description: Extracts a specified number of random bytes from a buffered stream.
*
* Arguments:   - uint8_t *out: pointer to the output byte array
* - size_t num_bytes: the number of bytes to extract
* - stream256_state *state: pointer to the stream state
* - uint8_t *buf: pointer to the buffer
* - size_t *buf_len: pointer to the buffer length
* - size_t *buf_pos: pointer to the current bit position in the buffer
**************************************************/
static void get_random_bytes(uint8_t *out, size_t num_bytes, stream256_state *state, uint8_t *buf, size_t *buf_len, size_t *buf_pos) {
    for(size_t i = 0; i < num_bytes; ++i) {
        uint8_t byte = 0;
        for(int bit_idx = 0; bit_idx < 8; ++bit_idx) {
             if (*buf_pos >= (*buf_len) * 8) {
                stream256_squeezeblocks(buf, 1, state);
                *buf_len = STREAM256_BLOCKBYTES;
                *buf_pos = 0;
            }
            byte |= ( (buf[*buf_pos / 8] >> (*buf_pos % 8)) & 1 ) << bit_idx;
            (*buf_pos)++;
        }
        out[i] = byte;
    }
}


/*************************************************
* Name:        sample_uniform_index_secure
*
* Description: Samples a uniform random index in the range [0, bound-1]
* using rejection sampling to avoid bias.
*
* Arguments:   - uint32_t bound: the upper bound (exclusive) for the index
* - stream256_state *state: pointer to the stream state
* - uint8_t *buf: pointer to the buffer for random bytes
* - size_t *buf_len: pointer to the buffer length
* - size_t *buf_pos: pointer to the current bit position
*
* Returns:     A uniformly random 32-bit integer in [0, bound-1].
**************************************************/
static uint32_t sample_uniform_index_secure(uint32_t bound, stream256_state *state, uint8_t *buf, size_t *buf_len, size_t *buf_pos) {
    if (bound <= 1) return 0;
    uint32_t threshold = (bound == 0) ? 0 : (0xFFFFFFFF / bound) * bound;
    if (threshold > 0xFFFFFFFF) threshold = 0;
    uint32_t random_u32;
    uint8_t temp_bytes[4];
    do {
        get_random_bytes(temp_bytes, 4, state, buf, buf_len, buf_pos);
        random_u32 = (uint32_t)temp_bytes[0] | ((uint32_t)temp_bytes[1] << 8) |
                     ((uint32_t)temp_bytes[2] << 16) | ((uint32_t)temp_bytes[3] << 24);
    } while (random_u32 >= threshold);
    return random_u32 % bound;
}


/*************************************************
* Name:        SampleSigma1
*
* Description: Samples a polynomial with coefficients uniformly from {-1, 1}.
*
* Arguments:   - poly *p: pointer to output polynomial
* - const uint8_t seed[]: pointer to input seed (of length CRHBYTES)
* - uint16_t nonce: 16-bit nonce
**************************************************/
void SampleSigma1(poly *p, const uint8_t seed[CRHBYTES], uint16_t nonce) {
    uint8_t buf[STREAM256_BLOCKBYTES];
    size_t buf_len = 0;
    size_t buf_pos = STREAM256_BLOCKBYTES * 8;
    stream256_state state;
    stream256_init(&state, seed, nonce);
    for (int i = 0; i < N; ++i) {
        uint8_t random_bit = get_random_bit(&state, buf, &buf_len, &buf_pos);
        p->coeffs[i] = 1 - 2 * (int32_t)random_bit;
    }
}


/*************************************************
* Name:        SampleX_poly
*
* Description: Samples the auxiliary polynomial x for the 3C sampling process.
* For each coefficient c_i of the challenge c, if c_i is non-zero,
* the corresponding coefficient x_i is sampled from {-2, 0}.
* This implementation achieves the same distribution for (c_i + x_i) as
* described in the paper for the case where c_i=1.
*
* Arguments:   - poly *x: pointer to output polynomial
* - const poly *c: pointer to the input challenge polynomial
* - const uint8_t seed[]: pointer to input seed (of length CRHBYTES)
* - uint16_t nonce: 16-bit nonce
**************************************************/
void SampleX_poly(poly *x, const poly *c, const uint8_t seed[CRHBYTES], uint16_t nonce) {
    uint8_t buf[STREAM256_BLOCKBYTES];
    size_t buf_len = 0;
    size_t buf_pos = STREAM256_BLOCKBYTES * 8;
    stream256_state state;
    stream256_init(&state, seed, nonce);

    for (int i = 0; i < N; ++i) {

        uint8_t r = get_random_bit(&state, buf, &buf_len, &buf_pos);

        // r=0 -> value=-2; r=1 -> value=0
        int32_t value = (1 - r) * -2;

        // ci=0 -> mask=0; ci=1 -> mask=-1 
        int32_t mask = -c->coeffs[i];

        // ci=0 -> xi = 0 & value = 0
        // ci=1 -> xi = -1 & value = value
        x->coeffs[i] = mask & value;
    }
}

/*************************************************
* Name:        sample_chi_0_1
*
* Description: Samples a single coefficient from the distribution Chi(0,1),
* which outputs -1, 0, 1 with probabilities 1/4, 1/2, 1/4.
*
* Arguments:   - stream256_state *state: pointer to the stream state
* - uint8_t *buf: pointer to the buffer for random bits
* - size_t *buf_len: pointer to the buffer length
* - size_t *buf_pos: pointer to the current bit position
*
* Returns:     A single coefficient sampled from Chi(0,1).
**************************************************/
static int32_t sample_chi_0_1(stream256_state *state, uint8_t *buf, size_t *buf_len, size_t *buf_pos) {
    uint8_t bit1 = get_random_bit(state, buf, buf_len, buf_pos);
    uint8_t bit2 = get_random_bit(state, buf, buf_len, buf_pos);

    int32_t selection_mask = (bit1 ^ bit2) - 1;

    int32_t value = (bit1 * 2) - 1;

    return selection_mask & value;
}

/*************************************************
* Name:        SampleY1
*
* Description: Samples the first part of the masking vector, y1. Each coefficient
* is sampled from the distribution Chi(0,1) + Chi(0,1).
*
* Arguments:   - polyvecl *y1: pointer to output vector
* - const uint8_t seed[]: pointer to input seed (of length CRHBYTES)
* - uint16_t nonce_base: base nonce for sampling
**************************************************/
void SampleY1(polyvecl *y1, const uint8_t seed[CRHBYTES], uint16_t nonce_base) {
     uint8_t buf[STREAM256_BLOCKBYTES];
     size_t buf_len = 0;
     size_t buf_pos = STREAM256_BLOCKBYTES * 8;
     stream256_state state;
     stream256_init(&state, seed, nonce_base);
     for (int i = 0; i < K + L; ++i) {
         for (int j = 0; j < N; ++j) {
             int32_t term1 = sample_chi_0_1(&state, buf, &buf_len, &buf_pos);
             int32_t term2 = sample_chi_0_1(&state, buf, &buf_len, &buf_pos);
             y1->vec[i].coeffs[j] = term1 + term2;
         }
     }
}

/*************************************************
* Name:        SampleSPrime
*
* Description: Samples the sparse vector s', where s'_0 is the zero polynomial and
* each subsequent polynomial s'_i has Hamming weight 1 with a
* coefficient of either 1 or -1.
*
* Arguments:   - polyvecl *s_prime: pointer to output vector
* - const uint8_t seed[]: pointer to input seed (of length CRHBYTES)
* - uint16_t nonce_base: base nonce for sampling
**************************************************/
void SampleSPrime(polyvecl *s_prime, const uint8_t seed[CRHBYTES], uint16_t nonce_base) {
    uint8_t buf[STREAM256_BLOCKBYTES];
    size_t buf_len = 0;
    size_t buf_pos = STREAM256_BLOCKBYTES * 8;
    stream256_state state;
    stream256_init(&state, seed, nonce_base);
    poly_zero(&s_prime->vec[0]);
    for (int i = 1; i < K + L; ++i) {
        poly_zero(&s_prime->vec[i]);
        uint32_t j = sample_uniform_index_secure(N, &state, buf, &buf_len, &buf_pos);
        uint8_t sign_bit = get_random_bit(&state, buf, &buf_len, &buf_pos);
        int32_t value = 1 - 2 * (int32_t)sign_bit;
        if (j < N) {
            s_prime->vec[i].coeffs[j] = value;
        }
    }
}


/*************************************************
* Name:        SampleY
*
* Description: Samples the full masking vector y using the 3C sampling method
* as described in the Rhyme paper.
* Computes y = y1 + s'*(c'+x') + s_bar'*(c''+x'').
*
* Arguments:   - polyvecl *y: pointer to output vector
* - const polyvecl *s_prime: pointer to the sparse vector s'
* - const poly *s_bar_prime_0: pointer to the s_bar_prime_0 polynomial
* - const uint8_t seed_y1[]: seed for y1
* - const uint8_t seed_c_prime[]: seed for c'
* - const uint8_t seed_x_prime[]: seed for x'
* - const uint8_t seed_c_prime2[]: seed for c''
* - const uint8_t seed_x_prime2[]: seed for x''
* - uint16_t nonce: the nonce for this iteration
**************************************************/
void SampleY(polyvecl *y,
             const polyvecl *s_prime,
             const poly *s_bar_prime_0,         
             const uint8_t seed_y1[CRHBYTES],
             const uint8_t seed_c_prime[CRHBYTES],
             const uint8_t seed_x_prime[CRHBYTES],
             const uint8_t seed_c_prime2[CRHBYTES], 
             const uint8_t seed_x_prime2[CRHBYTES], 
             uint16_t nonce)
{
    
    polyvecl y1, y2_part1, y2_part2;
    poly c_prime, x_prime, c_plus_x_prime;
    poly c_prime2, x_prime2, c_plus_x_prime2;

    // 1. Sample y1
    SampleY1(&y1, seed_y1, nonce);

    
    SampleChallenge(&c_prime, seed_c_prime);
    SampleX_poly(&x_prime, &c_prime, seed_x_prime, nonce);
    poly_add(&c_plus_x_prime, &c_prime, &x_prime);
    for (int i = 0; i < K + L; ++i) {
        
        //poly_mul_integer(&y2_part1.vec[i], &s_prime->vec[i], &c_plus_x_prime);
        poly_mul_sparse_ternary(&y2_part1.vec[i], &s_prime->vec[i], &c_prime, &x_prime);
    }

    
    SampleChallenge(&c_prime2, seed_c_prime2);
    SampleX_poly(&x_prime2, &c_prime2, seed_x_prime2, nonce);
    poly_add(&c_plus_x_prime2, &c_prime2, &x_prime2);

    
    //poly_mul_integer(&y2_part2.vec[0], s_bar_prime_0, &c_plus_x_prime2);
    poly_mul_sparse_ternary(&y2_part2.vec[0], s_bar_prime_0, &c_prime2, &x_prime2);
    for (int i = 1; i < K + L; ++i) {
        poly_zero(&y2_part2.vec[i]);
    }

    
    polyvecl_add(y, &y1, &y2_part1);
    polyvecl_add(y, y, &y2_part2);
}


/*************************************************
* Name:        SampleChallenge
*
* Description: Samples a challenge polynomial c with a fixed number of
* non-zero coefficients (TAU), which are all 1.
*
* Arguments:   - poly *c: pointer to output polynomial
* - const uint8_t seed[]: pointer to input seed (of length CRHBYTES)
**************************************************/
void SampleChallenge(poly *c, const uint8_t seed[CRHBYTES]) {
    uint8_t buf[STREAM256_BLOCKBYTES];
    size_t buf_len = 0;
    size_t buf_pos = STREAM256_BLOCKBYTES * 8;
    stream256_state state;
    uint32_t list[N];
    int32_t i;
    uint32_t j_idx;
    uint32_t tmp;
    stream256_init(&state, seed, 0);
    for(i=0; i<N; ++i) list[i] = i;
    for(i=0; i < TAU; ++i) {
        j_idx = i + sample_uniform_index_secure(N - i, &state, buf, &buf_len, &buf_pos);
        if (j_idx < N) {
            tmp = list[i];
            list[i] = list[j_idx];
            list[j_idx] = tmp;
        }
    }
    poly_zero(c);
    for(i=0; i<TAU; ++i) {
        if (list[i] < N) {
           c->coeffs[list[i]] = 1;
        }
    }
}


/*************************************************
* Name:        rej_uniform
*
* Description: Uniformly sample coefficients in [0, q-1] using rejection sampling
* on a byte stream.
*
* Arguments:   - int16_t *a: pointer to output array of coefficients
* - unsigned int len: number of coefficients to be sampled
* - const uint8_t *buf: pointer to input byte array
* - unsigned int buflen: length of input byte array
* - unsigned int *pos: pointer to the current position in the byte array
*
* Returns:     The number of coefficients that have been successfully sampled.
**************************************************/
unsigned int rej_uniform(int16_t *a, unsigned int len, const uint8_t *buf, unsigned int buflen, unsigned int *pos) { 
    unsigned int ctr = 0;
    uint16_t val16;
    
    //const uint32_t threshold = 19 * Q; 
    const uint32_t threshold = (0xFFFF / Q) * Q;
    while (ctr < len && *pos + 2 <= buflen) {
        val16 = (uint16_t)buf[*pos] | ((uint16_t)buf[*pos + 1] << 8);
        *pos += 2;
        if (val16 < threshold) {
            
            a[ctr++] = (int16_t)(val16 % Q); 
        }
    }
    return ctr;
}

#define POLY_UNIFORM_NBLOCKS (((N * 16 / 8) + STREAM128_BLOCKBYTES - 1) / STREAM128_BLOCKBYTES + 1) // Estimate based on 16 bits per coeff needed for rej_uniform
#define POLY_UNIFORM_BUF_SIZE (POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES + 2) // +2 for safety margin for rej_uniform read

/*************************************************
* Name:        poly_uniform_ntt
*
* Description: Samples a polynomial with uniformly random coefficients
* in [0, q-1]. This function is intended for generating
* polynomials that will be used directly in the NTT domain.
*
* Arguments:   - poly *a: pointer to output polynomial
* - const uint8_t seed[]: pointer to input seed (of length SEEDBYTES)
* - uint16_t nonce: 16-bit nonce
**************************************************/
void poly_uniform_ntt(poly *a, const uint8_t seed[SEEDBYTES], uint16_t nonce) {

    unsigned int ctr = 0;
    unsigned int pos = 0;
    unsigned int buflen = 0;
    uint8_t buf[POLY_UNIFORM_BUF_SIZE];
    stream128_state state;

    stream128_init(&state, seed, nonce);

    while (ctr < N) {
        unsigned int bytes_available = buflen - pos;
        // Need 2 bytes for rej_uniform
        if (bytes_available < 2) {
            // Move remaining bytes (0 or 1) to the beginning
             if (bytes_available > 0) {
                 memmove(buf, buf + pos, bytes_available);
            }
            buflen = bytes_available;
            pos = 0;

            // Squeeze more blocks
            unsigned int nblocks_to_squeeze = POLY_UNIFORM_NBLOCKS;
            if (buflen + nblocks_to_squeeze * STREAM128_BLOCKBYTES > POLY_UNIFORM_BUF_SIZE) {
                // Adjust if buffer would overflow (shouldn't happen with calculation above)
                nblocks_to_squeeze = (POLY_UNIFORM_BUF_SIZE - buflen) / STREAM128_BLOCKBYTES;
            }

            if (nblocks_to_squeeze > 0) {
                 stream128_squeezeblocks(buf + buflen, nblocks_to_squeeze, &state);
                 buflen += nblocks_to_squeeze * STREAM128_BLOCKBYTES;
            } else if (buflen < 2) {
                // Error: Cannot get enough bytes
                fprintf(stderr, "Error: poly_uniform cannot squeeze enough bytes.\n");
                poly_zero(a); // Zero out polynomial on error
                return;
            }
        }

        // Call rej_uniform (expects buffer, buffer length, pointer to position)
        ctr += rej_uniform(a->coeffs + ctr, N - ctr, buf, buflen, &pos);
    }
    stream128_free(&state);
}