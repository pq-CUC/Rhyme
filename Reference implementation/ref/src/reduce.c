#include "reduce.h"
#include "params.h"
#include <stdint.h>
#include "config.h" 

// =============================================================================
// Rhyme Optimization: Fermat Reduction for Q = 257
// =============================================================================

/*************************************************
* Name:        montgomery_reduce
*
* Description: Optimized reduction for Q=257 using Fermat's Little Theorem property.
* Instead of standard Montgomery reduction which computes a*R^-1 mod q,
* this function computes a mod q directly using bitwise operations.
* Since R = 2^16 == 1 (mod 257), a*R^-1 == a (mod 257), making
* this mathematically equivalent but significantly faster.
*
* The logic relies on 256 == -1 (mod 257).
* Any integer a = Sum(b_i * 256^i) is congruent to Sum(b_i * (-1)^i).
*
* Arguments:   - int64_t a: input integer to be reduced (typically product of two 32-bit integers)
*
* Returns:     a mod q (in centered or standard range depending on sign).
**************************************************/
int32_t montgomery_reduce(int64_t a) {
    // We treat the 64-bit integer 'a' as a sequence of 8 bytes: b0, b1, ..., b7.
    // Based on 256 == -1 (mod 257), the modular reduction becomes an alternating sum:
    // a = b0 - b1 + b2 - b3 + b4 - b5 + b6 - b7 (mod 257)

    int16_t res = (int16_t)(a & 0xFF)             // + b0
                - (int16_t)((a >> 8)  & 0xFF)     // - b1
                + (int16_t)((a >> 16) & 0xFF)     // + b2
                - (int16_t)((a >> 24) & 0xFF)     // - b3
                + (int16_t)((a >> 32) & 0xFF)     // + b4
                - (int16_t)((a >> 40) & 0xFF)     // - b5
                + (int16_t)((a >> 48) & 0xFF)     // + b6
                - (int16_t)(a >> 56);             // - b7 (Keep sign for MSB)

    return (int32_t)res;
}

/*************************************************
* Name:        barrett_reduce
*
* Description: Fast reduction for 16-bit integers using Fermat property.
* Computes a mod q avoiding division or multiplication.
*
* Arguments:   - int16_t a: input integer to be reduced
*
* Returns:     a mod q.
**************************************************/
int16_t barrett_reduce(int16_t a) {
    // Decompose 16-bit 'a' into low byte and high byte.
    // a = low + high * 256
    // Since 256 == -1 (mod 257), this becomes:
    // a = low - high (mod 257)
    
    // Note: Use int16_t cast to ensure correct subtraction behavior
    int16_t t = (int16_t)(a & 0xFF) - (int16_t)(a >> 8);

    // Normalize result to be non-negative if it falls below zero
    if (t < 0) {
        t += Q;
    }

    return t;
}

/*************************************************
* Name:        caddq
*
* Description: For a given integer a, computes a + q if a is negative,
* and a otherwise. Used for constant-time modular arithmetic.
*
* Arguments:   - int32_t a: input integer
*
* Returns:     a mod q, in the range where a is congruent to the input.
**************************************************/
int32_t caddq(int32_t a) {
    a += (a >> 31) & Q;
    return a;
}

/*************************************************
* Name:        freeze
*
* Description: Full reduction of an integer modulo q to the standard range [0, q-1].
* Optimized to avoid the '%' operator (division).
*
* Arguments:   - int32_t a: input integer
*
* Returns:     a mod q, in the standard range [0, q-1].
**************************************************/
int32_t freeze(int32_t a) {
    // 1. Fast folding: Reduce 32-bit integer to a small range using byte-wise alternating sum.
    // a = b0 - b1 + b2 - b3 (mod 257)
    int32_t t = (int32_t)(a & 0xFF) 
              - (int32_t)((a >> 8) & 0xFF) 
              + (int32_t)((a >> 16) & 0xFF) 
              - (int32_t)(a >> 24); // Keep sign for MSB

    // 2. Range correction: Bring result into [0, 256].
    // The value 't' is bounded, so loops run very few times (constant-time characteristics).
    
    while (t < 0) {
        t += Q;
    }
    while (t >= Q) {
        t -= Q;
    }

    return t;
}

/*************************************************
* Name:        cryptolab_ccc_reduce32_2q_centered
*
* Description: Centered reduction of an integer modulo 2q (514).
* Used for specific signature operations requiring mod 2q arithmetic.
*
* Arguments:   - int32_t a: input integer
*
* Returns:     a mod 2q, in the centered range (-q, q].
**************************************************/
int32_t cryptolab_ccc_reduce32_2q_centered(int32_t a) {
    int32_t r = a % DQ;
    if (r < 0) r += DQ;
    if (r > Q) r -= DQ;
    return r;
}

/*************************************************
* Name:        cryptolab_ccc_freeze2q_nonnegative
*
* Description: Non-negative reduction of an integer modulo 2q.
*
* Arguments:   - int32_t a: input integer
*
* Returns:     a mod 2q, in the standard range [0, 2q-1].
**************************************************/
int32_t cryptolab_ccc_freeze2q_nonnegative(int32_t a) {
    int32_t r = a % DQ;
    if (r < 0) r += DQ;
    return r;
}