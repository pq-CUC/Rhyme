#include "reduce.h" 
#include "params.h"
#include <stdint.h>
#include <stdio.h>
#include "config.h" 

/*************************************************
* Name:        montgomery_reduce
*
* Description: Montgomery reduction; given a 32-bit integer a, computes
* a*R^{-1} mod q, where R=2^16.
*
* Arguments:   - int64_t a: input integer to be reduced
*
* Returns:     montgomery_reduce(a) mod q.
**************************************************/
int32_t montgomery_reduce(int64_t a) {
    
    int32_t t;
    int16_t v; 

    v = a * QINV; 
    t = (int64_t)v * Q; 
    t = a - t;
    t >>= 16; 
    return t;
    
}

/*************************************************
* Name:        barrett_reduce
*
* Description: Barrett reduction; given a 16-bit integer a, computes a mod q
* in a way that is faster than a naive % operator.
*
* Arguments:   - int16_t a: input integer to be reduced
*
* Returns:     a mod q.
**************************************************/
int16_t barrett_reduce(int16_t a) {
    int32_t t;
    const int32_t v = 65281; 
    t = ((int64_t)a * v) >> 24;
    t = a - (t * Q);

    return (int16_t)t;
}

/*************************************************
* Name:        caddq
*
* Description: For a given integer a, computes a + q if a is negative,
* and a otherwise.
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
* Description: Full reduction of an integer modulo q.
*
* Arguments:   - int32_t a: input integer
*
* Returns:     a mod q, in the standard range [0, q-1].
**************************************************/
int32_t freeze(int32_t a) {
    int32_t r = a % Q;
    if (r < 0) {
        r += Q;
    }
    return r;
}

/*************************************************
* Name:        cryptolab_ccc_reduce32_2q_centered
*
* Description: Centered reduction of an integer modulo 2q.
*
* Arguments:   - int32_t a: input integer
*
* Returns:     a mod 2q, in the centered range (-q, q].
**************************************************/
int32_t cryptolab_ccc_reduce32_2q_centered(int32_t a) {
    int32_t r = a % DQ;
    if (r < 0) r += DQ;
    if (r > Q) r -= DQ;
    // Now r is in [-Q+1, Q]
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
    // Now r is in [0, 2Q-1]
    return r;
}