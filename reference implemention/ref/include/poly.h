#ifndef CCC_POLY_H
#define CCC_POLY_H

#include "params.h"
#include "reduce.h"
#include <stddef.h>
#include <stdint.h>
#include "config.h" // For ccc_NAMESPACE

typedef struct {
    int16_t coeffs[N];
} poly;

#define poly_add ccc_NAMESPACE(poly_add)
void poly_add(poly *c, const poly *a, const poly *b);

#define poly_sub ccc_NAMESPACE(poly_sub)
void poly_sub(poly *c, const poly *a, const poly *b);

#define poly_pointwise_montgomery ccc_NAMESPACE(poly_pointwise_montgomery)
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b);

#define poly_reduce2q ccc_NAMESPACE(poly_reduce2q)
void poly_reduce2q(poly *a);

#define poly_freeze2q ccc_NAMESPACE(poly_freeze2q)
void poly_freeze2q(poly *a);

#define poly_reduce ccc_NAMESPACE(poly_reduce)
void poly_reduce(poly *p);

#define poly_freeze ccc_NAMESPACE(poly_freeze)
void poly_freeze(poly *a);

#define poly_uniform ccc_NAMESPACE(poly_uniform)
void poly_uniform(poly *a, const uint8_t seed[SEEDBYTES], uint16_t nonce);

#define poly_ntt ccc_NAMESPACE(poly_ntt)
void poly_ntt(poly *a);

#define poly_invntt_tomont ccc_NAMESPACE(poly_invntt_tomont)
void poly_invntt_tomont(poly *a);

#define poly_zero ccc_NAMESPACE(poly_zero)
void poly_zero(poly *p);

#define poly_mul_const ccc_NAMESPACE(poly_mul_const)
void poly_mul_const(poly *r, const poly *a, int32_t factor);

#define poly_compare ccc_NAMESPACE(poly_compare)
int poly_compare(const poly *a, const poly *b);

#define poly_mul ccc_NAMESPACE(poly_mul)
void poly_mul(poly *c, const poly *a, const poly *b);

void ccc_NAMESPACE(polyq_pack)(uint8_t *r, const poly *a);
void ccc_NAMESPACE(polyq_unpack)(poly *r, const uint8_t *a);

#define polyq_pack ccc_NAMESPACE(polyq_pack)
#define polyq_unpack ccc_NAMESPACE(polyq_unpack)

#define poly_mul_integer ccc_NAMESPACE(poly_mul_integer)
void poly_mul_integer(poly *c, const poly *a, const poly *b);

#define poly_tomont ccc_NAMESPACE(poly_tomont)
void poly_tomont(poly *p); 

#define poly_caddq ccc_NAMESPACE(poly_caddq)
void poly_caddq(poly *p);

#define poly_double_coeffs_standard ccc_NAMESPACE(poly_double_coeffs_standard)
void poly_double_coeffs_standard(poly *p_out, const poly *p_in);

#define poly_crt_lift_from_modq_lsb0 ccc_NAMESPACE(poly_crt_lift_from_modq_lsb0)
void poly_crt_lift_from_modq_lsb0(poly *w, const poly *u_mod_q);

#define poly_mod_2 ccc_NAMESPACE(poly_mod_2)
void poly_mod_2(poly *r, const poly *a);

#define poly_crt_reconstruct_centered_mod_2Q ccc_NAMESPACE(poly_crt_reconstruct_centered_mod_2Q)
void poly_crt_reconstruct_centered_mod_2Q(poly *res_2Q, const poly *res_2, const poly *res_Q_frozen);


#endif // CCC_POLY_H