#ifndef CCC_POLYVEC_H // Changed guard name
#define CCC_POLYVEC_H

#include "params.h"
#include "poly.h"
#include <stdint.h>

/* Vectors of polynomials of length K */
typedef struct {
    poly vec[K];
} polyveck;

#define polyveck_add ccc_NAMESPACE(polyveck_add)
void polyveck_add(polyveck *w, const polyveck *u, const polyveck *v);
#define polyveck_sub ccc_NAMESPACE(polyveck_sub)
void polyveck_sub(polyveck *w, const polyveck *u, const polyveck *v);
#define polyveck_double ccc_NAMESPACE(polyveck_double)
void polyveck_double(polyveck *b); // Multiply by 2 mod 2Q

#define polyveck_reduce2q ccc_NAMESPACE(polyveck_reduce2q)
void polyveck_reduce2q(polyveck *v); // Apply centered reduce mod 2q
#define polyveck_freeze2q ccc_NAMESPACE(polyveck_freeze2q)
void polyveck_freeze2q(polyveck *v); // Apply non-negative reduce mod 2q
#define polyveck_freeze ccc_NAMESPACE(polyveck_freeze)
void polyveck_freeze(polyveck *v); // Apply reduce mod q -> [0, Q-1]

#define polyveck_expand ccc_NAMESPACE(polyveck_expand)
void polyveck_expand(polyveck *v, const uint8_t seed[SEEDBYTES]);

#define polyveck_ntt ccc_NAMESPACE(polyveck_ntt)
void polyveck_ntt(polyveck *x);
#define polyveck_invntt_tomont ccc_NAMESPACE(polyveck_invntt_tomont)
void polyveck_invntt_tomont(polyveck *x);

#define polyveck_pointwise_montgomery ccc_NAMESPACE(polyveck_pointwise_montgomery) // Renamed from polyveck_poly_pointwise_montgomery for consistency
void polyveck_pointwise_montgomery(polyveck *w, const polyveck *u, const poly *v);

#define polyveck_reduce ccc_NAMESPACE(polyveck_reduce)
void polyveck_reduce(polyveck *v); // Apply reduce mod q

#define polyveck_mod_2q ccc_NAMESPACE(polyveck_mod_2q)
void polyveck_mod_2q(polyveck *v); // Apply non-negative reduce mod 2q (same as freeze2q)

#define polyveck_negate ccc_NAMESPACE(polyveck_negate) // Added declaration
void polyveck_negate(polyveck *neg_v, const polyveck *v); // Negate coefficients

/* Vectors of polynomials of length L+K */ // **** UPDATED dimension ****
typedef struct {
    poly vec[L+K];
} polyvecl;

#define polyvecl_add ccc_NAMESPACE(polyvecl_add)
void polyvecl_add(polyvecl *w, const polyvecl *u, const polyvecl *v);

#define polyvecl_reduce ccc_NAMESPACE(polyvecl_reduce)
void polyvecl_reduce(polyvecl *v); // Reduce mod q

#define polyvecl_ntt ccc_NAMESPACE(polyvecl_ntt)
void polyvecl_ntt(polyvecl *x);

#define polyvecl_check_norm_inf ccc_NAMESPACE(polyvecl_check_norm_inf)
int polyvecl_check_norm_inf(const polyvecl *v, int32_t bound);

#define polyvecl_pointwise_acc_montgomery ccc_NAMESPACE(polyvecl_pointwise_acc_montgomery)
void polyvecl_pointwise_acc_montgomery(poly *w, const polyvecl *u, const polyvecl *v);


typedef struct {
  poly vec[M];
} polyvecm;

#define polyvecm_ntt ccc_NAMESPACE(polyvecm_ntt)
void polyvecm_ntt(polyvecm *x);

#define polyvecm_pointwise_acc_montgomery ccc_NAMESPACE(polyvecm_pointwise_acc_montgomery)
void polyvecm_pointwise_acc_montgomery(poly *w, const polyvecm *u, const polyvecm *v);

#define polyvecl_zero ccc_NAMESPACE(polyvecl_zero)
void polyvecl_zero(polyvecl *v);
#define polyvecl_invntt_tomont ccc_NAMESPACE(polyvecl_invntt_tomont)
void polyvecl_invntt_tomont(polyvecl *v);
// Declaration for the function to compute Ay mod 2q
void calculate_Ay(polyveck *w, const polyvecm Agen[K], const polyveck *Col0, const polyvecl *y);
void calculate_Av(polyveck *w, const polyvecm Agen[K], const polyveck *Col0, const polyvecl *v); // Renamed y to v
void calculate_A_vec_prod_crt(polyveck *w_out,
  const polyvecm Agen[K],   // K x M matrix
  const polyveck *Col0_orig, // K-dim vector, coeffs [0, 2Q-1]
  const polyvecl *z_input);  // (L+K)-dim vector, coeffs [-B_infty, B_infty]
#endif // CCC_POLYVEC_H