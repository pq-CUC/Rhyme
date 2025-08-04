#include "params.h" // Include params first
#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "sampler.h"   // Need rej_uniform for poly_uniform
#include "symmetric.h"
#include <stdint.h>
#include <string.h> // For memset, memcmp
#include <stdlib.h> // For calloc/free if used in poly_mul
#include <stdio.h>  // For fprintf in error cases

#ifndef STREAM128_BLOCKBYTES // Define if not already defined
#define STREAM128_BLOCKBYTES SHAKE128_RATE

#endif

/*************************************************
* Name:        poly_add
*
* Description: Add two polynomials. No modular reduction is performed.
*
* Arguments:   - poly *c: pointer to output polynomial
* - const poly *a: pointer to first operand
* - const poly *b: pointer to second operand
**************************************************/
void poly_add(poly *c, const poly *a, const poly *b) {
    for (unsigned int i = 0; i < N; ++i)
        c->coeffs[i] = a->coeffs[i] + b->coeffs[i];
    // Note: No reduction performed here. Caller must reduce if needed.
}

/*************************************************
* Name:        poly_sub
*
* Description: Subtract one polynomial from another. No modular reduction is performed.
*
* Arguments:   - poly *c: pointer to output polynomial
* - const poly *a: pointer to first operand
* - const poly *b: pointer to second operand
**************************************************/
void poly_sub(poly *c, const poly *a, const poly *b) {
    for (unsigned int i = 0; i < N; ++i)
        c->coeffs[i] = a->coeffs[i] - b->coeffs[i];
    // Note: No reduction performed here. Caller must reduce if needed.
}

/*************************************************
 * Name:        poly_pointwise_montgomery
 * Description: Pointwise multiplication in NTT domain (Montgomery form).
 * Assumes inputs a, b are in Montgomery form. Output c is too.
 **************************************************/
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b) {
    // for (unsigned int i = 0; i < N; ++i)
    //     //c->coeffs[i] = (int32_t)product;
    
    // // Result is in [-Q+1, Q-1]
    unsigned int i;
    for(i=0;i<N/4;i++) {
      basemul(&c->coeffs[4*i], &a->coeffs[4*i], &b->coeffs[4*i], zetas[64+i]);
      basemul(&c->coeffs[4*i+2], &a->coeffs[4*i+2], &b->coeffs[4*i+2], -zetas[64+i]);
    }
}

/*************************************************
* Name:        poly_reduce2q
*
* Description: Apply centered reduction to all coefficients of a polynomial
* for modulus 2q. The result is in the range (-q, q].
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_reduce2q(poly *a) {
    for (unsigned int i = 0; i < N; ++i)
        a->coeffs[i] = reduce32_2q(a->coeffs[i]); // Macro points to centered version
}

/*************************************************
* Name:        poly_freeze2q
*
* Description: Apply non-negative reduction to all coefficients of a polynomial
* for modulus 2q. The result is in the range [0, 2q-1].
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_freeze2q(poly *a) {
    for (unsigned int i = 0; i < N; ++i)
        a->coeffs[i] = freeze2q(a->coeffs[i]); // Macro points to non-negative version
}

/*************************************************
* Name:        poly_reduce
*
* Description: Apply reduction to all coefficients of a polynomial
* for modulus q. The result is in the range [0, q-1].
*
* Arguments:   - poly *p: pointer to input/output polynomial
**************************************************/
void poly_reduce(poly *p) {
    for (int i = 0; i < N; ++i) {
        p->coeffs[i] = freeze(p->coeffs[i]); // freeze results in [0, Q-1]
    }
}

/*************************************************
* Name:        poly_freeze
*
* Description: Apply reduction to all coefficients of a polynomial
* for modulus q. The result is in the range [0, q-1].
* (Alias for poly_reduce)
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_freeze(poly *a) {
    poly_reduce(a);
}

// Define buffer size calculation based on N and rate (adjust if needed)
#define POLY_UNIFORM_NBLOCKS (((N * 16 / 8) + STREAM128_BLOCKBYTES - 1) / STREAM128_BLOCKBYTES + 1) // Estimate based on 16 bits per coeff needed for rej_uniform
#define POLY_UNIFORM_BUF_SIZE (POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES + 2) // +2 for safety margin for rej_uniform read
/*************************************************
* Name:        poly_uniform
*
* Description: Sample a polynomial with uniformly random coefficients
* in [0, q-1] using rejection sampling on the output of
* a XOF.
*
* Arguments:   - poly *a: pointer to output polynomial
* - const uint8_t seed[]: pointer to input seed
* (of length SEEDBYTES)
* - uint16_t nonce: 16-bit nonce
**************************************************/
void poly_uniform(poly *a, const uint8_t seed[SEEDBYTES], uint16_t nonce) {
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
}

/*************************************************
* Name:        poly_ntt
*
* Description: Inplace forward Number-Theoretic Transform (NTT).
* Output coefficients are in bitreversed order and reduced
* in the range [0, q-1].
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_ntt(poly *a) {
    ntt(a->coeffs); // Assumes ntt modifies coeffs in-place
    poly_reduce(a); // Reduce result mod q -> [0, Q-1] (standard practice)
}

/*************************************************
* Name:        poly_invntt_tomont
*
* Description: Inplace inverse Number-Theoretic Transform (NTT) and
* multiplication by Montgomery factor.
* Input coefficients are in bitreversed order.
* Output coefficients are in standard order and in range [-q+1, q-1].
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_invntt_tomont(poly *a) {
    invntt_tomont(a->coeffs); // Assumes invntt_tomont modifies in-place
    // Result is in normal domain, coefficients likely in [-Q+1, Q-1]
}

/*************************************************
* Name:        poly_zero
*
* Description: Set all coefficients of a polynomial to zero.
*
* Arguments:   - poly *p: pointer to polynomial
**************************************************/
void poly_zero(poly *p) {
    memset(p->coeffs, 0, N * sizeof(p->coeffs[0]));
}

/*************************************************
* Name:        poly_mul_const
*
* Description: Multiply a polynomial by a constant integer factor.
* The result coefficients are reduced centered modulo 2q.
*
* Arguments:   - poly *r: pointer to output polynomial
* - const poly *a: pointer to input polynomial
* - int32_t factor: constant integer factor
**************************************************/
void poly_mul_const(poly *r, const poly *a, int32_t factor) {
     for (int i = 0; i < N; ++i) {
         int64_t temp = (int64_t)a->coeffs[i] * factor;
         // Use centered reduction mod 2q
         r->coeffs[i] = reduce32_2q((int32_t)(temp % DQ));
     }
 }

/*************************************************
* Name:        poly_compare
*
* Description: Compare two polynomials for equality.
*
* Arguments:   - const poly *a: pointer to first polynomial
* - const poly *b: pointer to second polynomial
*
* Returns:     0 if equal, non-zero otherwise.
**************************************************/
int poly_compare(const poly *a, const poly *b) {
    return memcmp(a->coeffs, b->coeffs, N * sizeof(a->coeffs[0]));
}


/*************************************************
* Name:        poly_mul
*
* Description: Schoolbook multiplication of two polynomials in R_2q.
* c = a * b mod (x^n + 1).
* The result coefficients are reduced centered modulo 2q.
*
* Arguments:   - poly *c: pointer to output polynomial
* - const poly *a: pointer to first factor
* - const poly *b: pointer to second factor
**************************************************/
void poly_mul(poly *c, const poly *a, const poly *b) {
    // Use stack allocation for result buffer if 2N-1 is not too large
    int64_t res[2 * N - 1]; // For N=256, this is 511*8 = 4088 bytes, likely OK on stack
    memset(res, 0, sizeof(res));
    int i, j;

    // Schoolbook multiplication
    for (i = 0; i < N; ++i) {
        if (a->coeffs[i] == 0) continue; // Optimization
        for (j = 0; j < N; ++j) {
            res[i + j] += (int64_t)a->coeffs[i] * b->coeffs[j];
        }
    }

    // Reduction using x^N = -1
    for (i = N; i < 2 * N - 1; ++i) {
        res[i - N] -= res[i]; // res[0] to res[N-1] now hold result before mod 2q
    }

    // Final reduction mod 2q (centered)
    // Final reduction mod 2q (centered) - CORRECTED VERSION
    for (i = 0; i < N; ++i) {
        int64_t val = res[i]; 

        int64_t r64 = val % DQ;
        
        if (r64 < 0) {
            r64 += DQ;
        }
     
        if (r64 > Q) {
            r64 -= DQ; 
        }

        c->coeffs[i] = (int32_t)r64;
    }
}

/*************************************************
* Name:        polyq_pack
*
* Description: Pack a polynomial with coefficients in [0, q-1].
*
* Arguments:   - uint8_t *r: pointer to output byte array
* (of length POLYQ_PACKEDBYTES)
* - const poly *a: pointer to input polynomial
**************************************************/
void polyq_pack(uint8_t *r, const poly *a) {
    size_t out_idx = 0;
    uint64_t bit_buffer = 0;
    int bits_in_buffer = 0;
    const int bits_per_coeff = POLYQ_BITS; // 12

    for (int i = 0; i < N; ++i) {
        // Assume a->coeffs[i] is already in [0, Q-1] by poly_reduce/freeze
        uint32_t coeff = (uint32_t)a->coeffs[i];
        if (coeff >= Q) { // Safety check
             // coeff = Q - 1; // Or handle error
        }
        bit_buffer |= ((uint64_t)coeff << bits_in_buffer);
        bits_in_buffer += bits_per_coeff;
        while (bits_in_buffer >= 8) {
            if (out_idx >= POLYQ_PACKEDBYTES) { /*fprintf(stderr, "Overflow in polyq_pack\n");*/ return; }
            r[out_idx++] = (uint8_t)(bit_buffer & 0xFF);
            bit_buffer >>= 8;
            bits_in_buffer -= 8;
        }
    }
    if (bits_in_buffer > 0) {
        if (out_idx >= POLYQ_PACKEDBYTES) { /*fprintf(stderr, "Overflow in polyq_pack (final byte)\n");*/ return; }
        r[out_idx++] = (uint8_t)(bit_buffer & 0xFF);
    }
    // Zero padding? Usually not needed if POLYQ_PACKEDBYTES is calculated correctly
    // while (out_idx < POLYQ_PACKEDBYTES) r[out_idx++] = 0;
}

/*************************************************
* Name:        polyq_unpack
*
* Description: Unpack a polynomial with coefficients in [0, q-1].
*
* Arguments:   - poly *r: pointer to output polynomial
* - const uint8_t *a: pointer to input byte array
* (of length POLYQ_PACKEDBYTES)
**************************************************/
void polyq_unpack(poly *r, const uint8_t *a) {
    size_t in_idx = 0;
    uint64_t bit_buffer = 0;
    int bits_in_buffer = 0;
    const int bits_per_coeff = POLYQ_BITS; // 12
    const uint32_t mask = (1U << bits_per_coeff) - 1; // 0xFFF

    for (int i = 0; i < N; ++i) {
        while (bits_in_buffer < bits_per_coeff) {
            if (in_idx >= POLYQ_PACKEDBYTES) { /* Error handling */ poly_zero(r); return; }
            bit_buffer |= ((uint64_t)a[in_idx++] << bits_in_buffer);
            bits_in_buffer += 8;
        }
        uint32_t coeff = (uint32_t)(bit_buffer & mask);
        bit_buffer >>= bits_per_coeff;
        bits_in_buffer -= bits_per_coeff;

        if (coeff < Q) {
            r->coeffs[i] = (int32_t)coeff;
        } else {
            // Invalid coefficient found during unpack
            // Set to 0 or Q-1? Or handle error? Let's set to 0 for now.
            r->coeffs[i] = 0;
             // fprintf(stderr, "Warning: Invalid coefficient %u unpacked in polyq_unpack\n", coeff);
        }
    }
}
/*************************************************
* Name:        karatsuba_mul
*
* Description: Recursive Karatsuba multiplication for integer polynomials.
* This is a helper function for poly_mul_integer.
*
* Arguments:   - int64_t *r: pointer to the output coefficient array (size 2n-1)
* - const int16_t *a: pointer to the first input polynomial's coefficients
* - const int16_t *b: pointer to the second input polynomial's coefficients
* - size_t n: the number of coefficients in the input polynomials
**************************************************/
static void karatsuba_mul(int64_t *r, const int16_t *a, const int16_t *b, size_t n) {
    
    if (n <= 32) {
        for (size_t i = 0; i < 2 * n - 1; ++i) r[i] = 0;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                r[i + j] += (int64_t)a[i] * b[j];
            }
        }
        return;
    }

    size_t n_half = n / 2;

    
    const int16_t *a0 = a, *a1 = a + n_half;
    const int16_t *b0 = b, *b1 = b + n_half;
    
    int64_t *r0 = calloc(n - 1, sizeof(int64_t));
    int64_t *r2 = calloc(n - 1, sizeof(int64_t));
    int64_t *r1_temp = calloc(n - 1, sizeof(int64_t));
    
    int16_t *a_sum = malloc(n_half * sizeof(int16_t));
    int16_t *b_sum = malloc(n_half * sizeof(int16_t));
  
    if (!r0 || !r2 || !r1_temp || !a_sum || !b_sum) {
        
        free(r0); free(r2); free(r1_temp); free(a_sum); free(b_sum);
        return;
    }

    for(size_t i = 0; i < n_half; ++i) a_sum[i] = a0[i] + a1[i];
    for(size_t i = 0; i < n_half; ++i) b_sum[i] = b0[i] + b1[i];

    karatsuba_mul(r0, a0, b0, n_half);             // z0 = a0 * b0
    karatsuba_mul(r2, a1, b1, n_half);             // z2 = a1 * b1
    karatsuba_mul(r1_temp, a_sum, b_sum, n_half);  // (a0+a1)*(b0+b1)

    // r = z2*x^n + ( (a0+a1)(b0+b1) - z0 - z2 )*x^(n/2) + z0
    for(size_t i = 0; i < n - 1; ++i) {
        r[i] += r0[i];
        r[i + n_half] += r1_temp[i] - r0[i] - r2[i];
        r[i + n] += r2[i];
    }
    
    free(r0);
    free(r2);
    free(r1_temp);
    free(a_sum);
    free(b_sum);
}

/*************************************************
* Name:        poly_mul_integer
*
* Description: Integer multiplication of two polynomials without modular reduction.
* Uses Karatsuba multiplication for efficiency.
* c = a * b mod (x^n + 1).
*
* Arguments:   - poly *c: pointer to output polynomial
* - const poly *a: pointer to first factor
* - const poly *b: pointer to second factor
**************************************************/
void poly_mul_integer(poly *c, const poly *a, const poly *b) {
    
    int64_t res[2 * N - 1] = {0};

    
    karatsuba_mul(res, a->coeffs, b->coeffs, N);
    
    
    for (int i = N; i < 2 * N - 1; ++i) {
        res[i - N] -= res[i];
    }

    
    for (int i = 0; i < N; ++i) {
        
        c->coeffs[i] = (int32_t)res[i];
    }
}

/*************************************************
* Name:        poly_caddq
*
* Description: Add q to all negative coefficients of a polynomial.
*
* Arguments:   - poly *p: pointer to input/output polynomial
**************************************************/
void poly_caddq(poly *p) {
    for (int i = 0; i < N; ++i) {
        p->coeffs[i] = caddq(p->coeffs[i]); 
    }
}

/*************************************************
* Name:        poly_double_coeffs_standard
*
* Description: Doubles all coefficients of a polynomial and reduces the result
* modulo q to the standard range [0, q-1].
*
* Arguments:   - poly *p_out: pointer to the output polynomial
* - const poly *p_in: pointer to the input polynomial
**************************************************/
void poly_double_coeffs_standard(poly *p_out, const poly *p_in) {
    for (int i = 0; i < N; ++i) {
        
        p_out->coeffs[i] = freeze(2 * p_in->coeffs[i]);
    }
}

/*************************************************
* Name:        poly_crt_lift_from_modq_lsb0
*
* Description: "Lifts" a polynomial with coefficients mod q to have coefficients
* mod 2q, with the additional constraint that all resulting
* coefficients are even.
*
* Arguments:   - poly *w: pointer to the output polynomial (coeffs mod 2q, even)
* - const poly *u_mod_q: pointer to the input polynomial (coeffs mod q)
**************************************************/
void poly_crt_lift_from_modq_lsb0(poly *w, const poly *u_mod_q) {
    for (int i = 0; i < N; ++i) {
        int16_t u_coeff = u_mod_q->coeffs[i];
        // If u_coeff is odd, add Q to make it even while preserving mod Q value.
        // (u_coeff & 1) is 1 if odd, 0 if even.
        // (Q & -(u_coeff & 1)) is Q if odd, 0 if even.
        w->coeffs[i] = u_coeff + (Q & -(u_coeff & 1));
    }
}

/*************************************************
* Name:        poly_mod_2
*
* Description: Computes polynomial coefficients modulo 2.
* Output coefficients are 0 or 1.
*
* Arguments:   - poly *r: pointer to output polynomial
* - const poly *a: pointer to input polynomial
**************************************************/
void poly_mod_2(poly *r, const poly *a) {
    for (unsigned int i = 0; i < N; ++i) {
        
        
        
        r->coeffs[i] = (a->coeffs[i] % 2 + 2) % 2;
    }
}

/*************************************************
* Name:        poly_crt_reconstruct_centered_mod_2Q
*
* Description: Reconstructs polynomial coefficients mod 2Q (centered in (-Q, Q])
* from coefficients mod 2 and mod Q using the Chinese Remainder Theorem.
*
* Arguments:   - poly *res_2Q: pointer to output polynomial (coeffs centered mod 2Q)
* - const poly *res_2: pointer to polynomial with coeffs mod 2 (0 or 1)
* - const poly *res_Q_frozen: pointer to polynomial with coeffs mod Q (in [0, Q-1])
**************************************************/
void poly_crt_reconstruct_centered_mod_2Q(poly *res_2Q, const poly *res_2, const poly *res_Q_frozen) {
    
    
    for (int k_coeff = 0; k_coeff < N; ++k_coeff) {
        
        int32_t val_Q = res_Q_frozen->coeffs[k_coeff]; 
        int32_t val_2 = res_2->coeffs[k_coeff];        
        
        int32_t val_Q_mod_2 = val_Q & 1;
        int32_t k_crt = (val_2 - val_Q_mod_2 + 2) & 1; 

        int32_t val_2Q_non_centered = val_Q + k_crt * Q;

        res_2Q->coeffs[k_coeff] = reduce32_2q(val_2Q_non_centered);
    }
}



