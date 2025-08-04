#include <stdint.h>
#include <stdlib.h> 
#include <stdio.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "reduce.h" 
#include "ntt.h"

/*************************************************
* Name:        poly_double_poly_centered
*
* Description: Helper function to multiply a polynomial by 2, with coefficients
* reduced centered modulo 2q.
*
* Arguments:   - poly *p: pointer to input/output polynomial
**************************************************/
static void poly_double_poly_centered(poly *p) {
     for (int j = 0; j < N; ++j) {
         int64_t temp = (int64_t)p->coeffs[j] * 2;
         // Use centered reduction macro/function
         p->coeffs[j] = reduce32_2q((int32_t)(temp % DQ));
     }
 }

/*************************************************
* Name:        polyveck_add
*
* Description: Add two polynomial vectors of length K.
*
* Arguments:   - polyveck *w: pointer to output vector
* - const polyveck *u: pointer to first operand
* - const polyveck *v: pointer to second operand
**************************************************/
void polyveck_add(polyveck *w, const polyveck *u, const polyveck *v) {
    for (unsigned int i = 0; i < K; ++i)
        poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyveck_sub
*
* Description: Subtract one polynomial vector of length K from another.
* No modular reduction is performed.
*
* Arguments:   - polyveck *w: pointer to output vector
* - const polyveck *u: pointer to first operand
* - const polyveck *v: pointer to second operand to be subtracted
**************************************************/
void polyveck_sub(polyveck *w, const polyveck *u, const polyveck *v) {
    for (unsigned int i = 0; i < K; ++i)
        poly_sub(&w->vec[i], &u->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyveck_double
*
* Description: Multiplies all coefficients of a polynomial vector of length K
* by 2. The result is reduced modulo 2Q into the
* non-negative range [0, 2Q-1].
*
* Arguments:   - polyveck *b: pointer to the input/output vector
**************************************************/
void polyveck_double(polyveck *b) {
    for (unsigned int i = 0; i < K; ++i) {
        for (unsigned int j = 0; j < N; ++j) {
            int64_t temp = (int64_t)b->vec[i].coeffs[j] * 2;
            // Use non-negative reduction macro/function
            b->vec[i].coeffs[j] = freeze2q((int32_t)(temp % DQ));
        }
    }
}

/*************************************************
* Name:        polyveck_reduce2q
*
* Description: Apply centered reduction to all coefficients of a
* polynomial vector of length K for modulus 2q.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_reduce2q(polyveck *v) {
    for (unsigned int i = 0; i < K; ++i) {
        for (unsigned int j = 0; j < N; ++j) {
             v->vec[i].coeffs[j] = reduce32_2q(v->vec[i].coeffs[j]); // Macro points to centered
        }
    }
}

/*************************************************
* Name:        polyveck_freeze2q
*
* Description: Apply non-negative reduction to all coefficients of a
* polynomial vector of length K for modulus 2q.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_freeze2q(polyveck *v) {
    for (unsigned int i = 0; i < K; ++i) {
         for (unsigned int j = 0; j < N; ++j) {
              v->vec[i].coeffs[j] = freeze2q(v->vec[i].coeffs[j]); // Macro points to non-negative
         }
    }
}

/*************************************************
* Name:        polyveck_mod_2q
*
* Description: Applies non-negative reduction to all coefficients of a
* polynomial vector of length K for modulus 2q.
* (Alias for polyveck_freeze2q)
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_mod_2q(polyveck *v) {
    polyveck_freeze2q(v);
}

/*************************************************
* Name:        polyveck_reduce
*
* Description: Applies reduction to all coefficients of a polynomial vector of
* length K for modulus q.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_reduce(polyveck *v) {
    for (int i = 0; i < K; ++i) {
        poly_reduce(&v->vec[i]);
    }
}

/*************************************************
* Name:        polyveck_freeze
*
* Description: Applies reduction to all coefficients of a polynomial vector of
* length K for modulus q. The result is in the range [0, q-1].
* (Alias for polyveck_reduce)
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_freeze(polyveck *v) {
    for (unsigned int i = 0; i < K; ++i)
        poly_freeze(&v->vec[i]);
}

/*************************************************
* Name:        polyveck_expand
*
* Description: Generates a polynomial vector of length K with uniformly random
* coefficients in [0, q-1].
*
* Arguments:   - polyveck *v: pointer to the output vector
* - const uint8_t seed[]: pointer to the input seed (of length SEEDBYTES)
**************************************************/
void polyveck_expand(polyveck *v, const uint8_t seed[SEEDBYTES]) {
    for (unsigned int i = 0; i < K; ++i) {
        poly_uniform(&v->vec[i], seed, (uint16_t)i); // Simple nonce scheme
    }
}

/*************************************************
* Name:        polyveck_ntt
*
* Description: Apply forward NTT to all elements of a polynomial vector of length K.
*
* Arguments:   - polyveck *x: pointer to input/output vector
**************************************************/
void polyveck_ntt(polyveck *x) {
    for (unsigned int i = 0; i < K; i++) {
        poly_ntt(&x->vec[i]);
    }
}

/*************************************************
* Name:        polyveck_invntt_tomont
*
* Description: Applies the inverse NTT to all elements of a polynomial vector of length K.
*
* Arguments:   - polyveck *x: pointer to the input/output vector
**************************************************/
void polyveck_invntt_tomont(polyveck *x) {
    for (unsigned int i = 0; i < K; i++) {
        poly_invntt_tomont(&x->vec[i]);
    }
}

/*************************************************
* Name:        polyveck_pointwise_montgomery
*
* Description: Pointwise multiplies a polynomial vector of length K by a single polynomial.
*
* Arguments:   - polyveck *w: pointer to the output vector
* - const polyveck *u: pointer to the input vector
* - const poly *v: pointer to the polynomial factor
**************************************************/
void polyveck_pointwise_montgomery(polyveck *w, const polyveck *u, const poly *v) {
    for (unsigned int i = 0; i < K; i++) {
        poly_pointwise_montgomery(&w->vec[i], &u->vec[i], v);
    }
}

/*************************************************
* Name:        polyveck_negate
*
* Description: Negates all coefficients of a polynomial vector of length K.
* No modular reduction is performed.
*
* Arguments:   - polyveck *neg_v: pointer to the output vector
* - const polyveck *v: pointer to the input vector
**************************************************/
void polyveck_negate(polyveck *neg_v, const polyveck *v) {
    for(int i=0; i<K; ++i) {
         for(int j=0; j<N; ++j) {
              neg_v->vec[i].coeffs[j] = -v->vec[i].coeffs[j];
         }
    }
}

/*************************************************
* Name:        polyvecl_add
*
* Description: Add two polynomial vectors of length L+K.
*
* Arguments:   - polyvecl *w: pointer to output vector
* - const polyvecl *u: pointer to first operand
* - const polyvecl *v: pointer to second operand
**************************************************/
void polyvecl_add(polyvecl *w, const polyvecl *u, const polyvecl *v) {
    for (int i = 0; i < L + K; ++i) // **** USE L+K ****
        poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyvecl_reduce
*
* Description: Applies reduction to all coefficients of a polynomial vector
* of length L+K for modulus q. The result for each
* coefficient is in the standard range [0, q-1].
*
* Arguments:   - polyvecl *v: pointer to the input/output vector
**************************************************/
void polyvecl_reduce(polyvecl *v) {
     for (int i = 0; i < L + K; ++i) { // **** USE L+K ****
         poly_reduce(&v->vec[i]); // Reduces mod q -> [0, Q-1]
     }
 }

/*************************************************
* Name:        polyvecl_check_norm_inf
*
* Description: Check the infinity norm of a polynomial vector of length L+K.
* The norm is computed on the integer coefficients without modular reduction.
*
* Arguments:   - const polyvecl *v: pointer to the vector
* - int32_t bound: the norm bound
*
* Returns:     1 if all coefficients are within bound, 0 otherwise.
**************************************************/
int polyvecl_check_norm_inf(const polyvecl *v, int32_t bound) {
    for (int i = 0; i < L + K; ++i) { // **** USE L+K ****
        for (int j = 0; j < N; ++j) {
            int32_t coeff = v->vec[i].coeffs[j]; 
            // Assume input v is now potentially large integer coefficients

            
            // if (coeff > Q / 2) coeff -= Q;

            
            if (abs(coeff) >= bound) { // Compare absolute value with bound
                // fprintf(stderr, "Debug: Norm check failed for integer coeff: |%d| >= %d\n", coeff, bound);
                return 0; // Failed
            }
        }
    }
    return 1; // Passed
}

/*************************************************
* Name:        polyvecl_ntt
*
* Description: Apply forward NTT to all elements of a polynomial vector of length L+K.
*
* Arguments:   - polyvecl *x: pointer to input/output vector
**************************************************/
void polyvecl_ntt(polyvecl *x) {
    for (unsigned int i = 0; i < L + K; i++) { // **** USE L+K ****
        poly_ntt(&x->vec[i]);
    }
}

/*************************************************
* Name:        polyvecl_pointwise_acc_montgomery
*
* Description: Computes the dot product of two polynomial vectors of length L+K.
* Result is w = u[0]*v[0] + u[1]*v[1] + ...
*
* Arguments:   - poly *w: pointer to the output polynomial
* - const polyvecl *u: pointer to the first input vector
* - const polyvecl *v: pointer to the second input vector
**************************************************/
void polyvecl_pointwise_acc_montgomery(poly *w, const polyvecl *u, const polyvecl *v) {
    poly t;
    poly_pointwise_montgomery(w, &u->vec[0], &v->vec[0]);
    for (unsigned int i = 1; i < L + K; ++i) { // **** USE L+K ****
        poly_pointwise_montgomery(&t, &u->vec[i], &v->vec[i]);
        poly_add(w, w, &t);
        // poly_reduce(w); // Optional: reduce intermediate sum mod q if needed
    }
    // Final result w is likely in Montgomery domain, potentially needs invntt_tomont later
}

/*************************************************
* Name:        polyvecm_ntt
*
* Description: Applies the forward NTT to all elements of a polynomial vector of length M.
*
* Arguments:   - polyvecm *x: pointer to the input/output vector
**************************************************/
void polyvecm_ntt(polyvecm *x) {
    for (unsigned int i = 0; i < M; i++) { // **** USE M ****
        poly_ntt(&x->vec[i]);
    }
}

/*************************************************
* Name:        polyvecm_pointwise_acc_montgomery
*
* Description: Computes the dot product of two polynomial vectors of length M.
* Result is w = u[0]*v[0] + u[1]*v[1] + ...
*
* Arguments:   - poly *w: pointer to the output polynomial
* - const polyvecm *u: pointer to the first input vector
* - const polyvecm *v: pointer to the second input vector
**************************************************/
void polyvecm_pointwise_acc_montgomery(poly *w, const polyvecm *u, const polyvecm *v) {
    poly t;
    poly_pointwise_montgomery(w, &u->vec[0], &v->vec[0]);
    for (unsigned int i = 1; i < M; ++i) { // **** USE M ****
        poly_pointwise_montgomery(&t, &u->vec[i], &v->vec[i]);
        poly_add(w, w, &t);
        // poly_reduce(w); // Optional intermediate reduction
    }
}

/*************************************************
* Name:        calculate_Ay
*
* Description: Computes the matrix-vector product w = A*y mod 2q, where A is
* the special matrix defined in the Rhyme scheme.
* Assumes input y has coefficients reduced modulo q.
* Result w has coefficients reduced centered modulo 2q.
*
* Arguments:   - polyveck *w: pointer to the output vector
* - const polyvecm Agen[K]: pointer to the A_gen part of matrix A
* - const polyveck *Col0: pointer to the first column of matrix A
* - const polyvecl *y: pointer to the input vector y
**************************************************/
void calculate_Ay(polyveck *w, const polyvecm Agen[K], const polyveck *Col0, const polyvecl *y) {
    poly temp_prod, term_sum;
    poly temp_op1, temp_op2; // Temporaries for operands

    for (int i = 0; i < K; ++i) {
        poly_zero(&w->vec[i]); // Initialize w_i

        // Term 1: Col0_i * y_0 (mod 2q)
        // Col0 assumed to be mod 2q (non-negative). y_0 assumed mod q.
        // poly_mul handles operands and gives centered result mod 2q.
        poly_mul(&term_sum, &Col0->vec[i], &y->vec[0]);
        poly_add(&w->vec[i], &w->vec[i], &term_sum);
        poly_reduce2q(&w->vec[i]); // Keep intermediate w_i centered mod 2q

        // Term 2: sum_{p=1}^{L-1} (2*Agen[i].vec[p-1]) * y_p (mod 2q)
        // Agen coeffs are mod q. y_p coeffs are mod q.
        poly_zero(&term_sum);
        for (int p = 1; p < L; ++p) { // Corresponds to y components 1 to L-1
            temp_op1 = Agen[i].vec[p-1]; // Copy Agen element (mod q)

            // Multiply temp_op1 by 2 mod 2q (centered result)
            poly_double_poly_centered(&temp_op1); // temp_op1 = 2 * Agen[i].vec[p-1] mod 2q (centered)

            temp_op2 = y->vec[p]; // Copy y element (mod q)
            // poly_reduce2q(&temp_op2); // Center y_p? Not needed if poly_mul handles mixed inputs

            poly_mul(&temp_prod, &temp_op1, &temp_op2); // temp_prod = (2*Agen)*y_p mod 2q (centered)
            poly_add(&term_sum, &term_sum, &temp_prod);
            poly_reduce2q(&term_sum); // Reduce intermediate sum (centered)
        }
        poly_add(&w->vec[i], &w->vec[i], &term_sum);
        poly_reduce2q(&w->vec[i]); // Reduce w_i (centered)

        // Term 3: sum_{p=L}^{L+K-1} (2*I_{i, p-L}) * y_p (mod 2q)
        // Simplifies to 2 * y_{L+i} (mod 2q)
        temp_op1 = y->vec[L + i]; // Copy y_{L+i} (mod q)
        // poly_reduce2q(&temp_op1); // Center y_{L+i}? Not needed.

        poly_double_poly_centered(&temp_op1); // temp_op1 = 2 * y_{L+i} mod 2q (centered)

        poly_add(&w->vec[i], &w->vec[i], &temp_op1);
        poly_reduce2q(&w->vec[i]); // Final centered reduction for w_i
    }
    // w now contains the result A*y mod 2q, with coefficients centered in (-Q, Q]
}

/*************************************************
* Name:        calculate_Av
*
* Description: Computes the matrix-vector product w = A*v mod 2q, where A is
* the special matrix defined in the Rhyme scheme.
* Accepts input vector v with arbitrary small integer coefficients.
* Result w has coefficients reduced centered modulo 2q.
*
* Arguments:   - polyveck *w: pointer to the output vector
* - const polyvecm Agen[K]: pointer to the A_gen part of matrix A
* - const polyveck *Col0: pointer to the first column of matrix A
* - const polyvecl *v: pointer to the input vector v
**************************************************/
void calculate_Av(polyveck *w, const polyvecm Agen[K], const polyveck *Col0, const polyvecl *v) { 
    poly temp_prod, term_sum;
    poly temp_op1, temp_op2; // Temporaries for operands

    for (int i = 0; i < K; ++i) {
        poly_zero(&w->vec[i]); // Initialize w_i

        
        // Col0 assumed to be mod 2q (non-negative). v_0 has small coeffs.
        // poly_mul handles operands and gives centered result mod 2q.
        poly_mul(&term_sum, &Col0->vec[i], &v->vec[0]); 
        poly_add(&w->vec[i], &w->vec[i], &term_sum);
        poly_reduce2q(&w->vec[i]); // Keep intermediate w_i centered mod 2q

        
        // Agen coeffs are mod q. v_p coeffs are small integers.
        poly_zero(&term_sum);
        for (int p = 1; p < L; ++p) { // Corresponds to v components 1 to L-1
            temp_op1 = Agen[i].vec[p-1]; // Copy Agen element (mod q)

            // Multiply temp_op1 by 2 mod 2q (centered result)
            poly_double_poly_centered(&temp_op1); // temp_op1 = 2 * Agen[i].vec[p-1] mod 2q (centered)

            temp_op2 = v->vec[p]; 
            // poly_reduce2q(&temp_op2); // Centering input v_p is not assumed/needed by poly_mul

            poly_mul(&temp_prod, &temp_op1, &temp_op2); // temp_prod = (2*Agen)*v_p mod 2q (centered)
            poly_add(&term_sum, &term_sum, &temp_prod);
            poly_reduce2q(&term_sum); // Reduce intermediate sum (centered)
        }
        poly_add(&w->vec[i], &w->vec[i], &term_sum);
        poly_reduce2q(&w->vec[i]); // Reduce w_i (centered)

        
        // Simplifies to 2 * v_{L+i} (mod 2q)
        temp_op1 = v->vec[L + i]; 
        // poly_reduce2q(&temp_op1); // Centering input v_{L+i} not assumed/needed

        poly_double_poly_centered(&temp_op1); // temp_op1 = 2 * v_{L+i} mod 2q (centered)

        poly_add(&w->vec[i], &w->vec[i], &temp_op1);
        poly_reduce2q(&w->vec[i]); // Final centered reduction for w_i
    }
    // w now contains the result A*v mod 2q, with coefficients centered in (-Q, Q]
}

/*************************************************
* Name:        polyvecl_invntt_tomont
*
* Description: Applies the inverse NTT to all elements of a polynomial vector
* of length L+K.
*
* Arguments:   - polyvecl *v: pointer to the input/output vector
**************************************************/
void polyvecl_invntt_tomont(polyvecl *v) {
    for (unsigned int i = 0; i < L + K; ++i) { 
        
        poly_invntt_tomont(&v->vec[i]);
    }
}


/*************************************************
* Name:        calculate_A_vec_prod_crt
*
* Description: Computes the matrix-vector product w = A*v using a mix of
* NTT-based multiplication for modulus q and direct computation
* for modulus 2, then combines the results using CRT.
*
* Arguments:   - polyveck *w_out: pointer to the output vector
* - const polyvecm Agen_ntt[K]: pointer to the NTT form of matrix A_gen
* - const polyveck *Col0_ntt: pointer to the NTT form of the first column of A
* - const polyvecl *z_input: pointer to the input vector
**************************************************/
void calculate_A_vec_prod_crt(polyveck *w, const polyvecm Agen_ntt[K], const polyveck *Col0_ntt, const polyvecl *v)
{
    polyvecl v_ntt;
    polyveck w_mod_q, w_mod_2;
    poly temp_prod_ntt;

    
    for(int i=0; i < L+K; ++i) {
        v_ntt.vec[i] = v->vec[i];
        //poly_reduce(&v_ntt.vec[i]);
        poly_ntt(&v_ntt.vec[i]);
    }

    for (int i = 0; i < K; ++i) {
        poly_pointwise_montgomery(&w_mod_q.vec[i], &Col0_ntt->vec[i], &v_ntt.vec[0]);
        
        for (int j = 0; j < M; ++j) {
            poly Agen_doubled_ntt;
            for(int l=0; l<N; ++l) {
                Agen_doubled_ntt.coeffs[l] = fqmul(Agen_ntt[i].vec[j].coeffs[l], mont_two);
            }
            poly_pointwise_montgomery(&temp_prod_ntt, &Agen_doubled_ntt, &v_ntt.vec[j + 1]);
            poly_add(&w_mod_q.vec[i], &w_mod_q.vec[i], &temp_prod_ntt);
        }

        poly v_doubled_ntt;
        for(int l=0; l<N; ++l) {
            v_doubled_ntt.coeffs[l] = fqmul(v_ntt.vec[L+i].coeffs[l], mont_two);
        }
        poly_add(&w_mod_q.vec[i], &w_mod_q.vec[i], &v_doubled_ntt);

        poly_invntt_tomont(&w_mod_q.vec[i]);
    }
    
    
    poly v0_mod_2;
    poly_mod_2(&v0_mod_2, &v->vec[0]);
    w_mod_2.vec[0] = v0_mod_2;
    for (int i = 1; i < K; ++i) {
        poly_zero(&w_mod_2.vec[i]);
    }

    
    for (int i = 0; i < K; ++i) {
        poly_crt_reconstruct_centered_mod_2Q(&w->vec[i], &w_mod_2.vec[i], &w_mod_q.vec[i]);
    }
}