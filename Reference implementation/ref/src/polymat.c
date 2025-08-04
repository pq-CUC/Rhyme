#include "polymat.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "sampler.h"
#include <stdint.h>

/*************************************************
 * Name:        polymat_expand
 *
 * Description: Implementation of ExpandA. Generates matrix A with uniformly
 *              random coefficients a_{i,j} by performing rejection
 *              sampling on the output stream of SHAKE128(rho|j|i)
 *              or AES256CTR(rho,j|i).
 *
 * Arguments:   - polyvecm mat[K]: output matrix k \times m
 *              - const uint8_t rho[]: byte array containing seed rho
 **************************************************/
void polymatkl_expand(polyvecl mat[K], const uint8_t rho[SEEDBYTES]) {
    unsigned int i, j;

    for (i = 0; i < K; ++i)
        for (j = 0; j < M; ++j)
            poly_uniform(&mat[i].vec[j + 1], rho, (i << 8) + j);
}

/*************************************************
 * Name:        polymat_expand
 *
 * Description: Implementation of ExpandA. Generates matrix A with uniformly
 *              random coefficients a_{i,j} by performing rejection
 *              sampling on the output stream of SHAKE128(rho|j|i)
 *              or AES256CTR(rho,j|i).
 *
 * Arguments:   - polyvecm mat[K]: output matrix k \times m
 *              - const uint8_t rho[]: byte array containing seed rho
 **************************************************/
void polymatkm_expand(polyvecm mat[K], const uint8_t rho[SEEDBYTES]) {
    unsigned int i, j;

    for (i = 0; i < K; ++i)
        for (j = 0; j < M; ++j)
            poly_uniform(&mat[i].vec[j], rho, (i << 8) + j);
}

/*************************************************
* Name:        polymatkl_double
*
* Description: Doubles the coefficients of the k x m sub-matrix within a k x l matrix.
* No modular reduction is performed.
*
* Arguments:   - polyvecl mat[K]: pointer to the input/output k x l matrix
**************************************************/
void polymatkl_double(polyvecl mat[K]) {
    unsigned int i, j, k;
    for (i = 0; i < K; ++i) {
        for (j = 1; j < L; ++j) {
            for (k = 0; k < N; ++k) {
                mat[i].vec[j].coeffs[k] *= 2;
            }
        }
    }
}

/*************************************************
* Name:        polymatkl_pointwise_montgomery
*
* Description: Pointwise multiplies a k x l matrix with a vector of length l.
* The result is a vector of length k.
*
* Arguments:   - polyveck *t: pointer to the output vector
* - const polyvecl mat[K]: pointer to the input k x l matrix
* - const polyvecl *v: pointer to the input vector of length l
**************************************************/
void polymatkl_pointwise_montgomery(polyveck *t, const polyvecl mat[K],
                                    const polyvecl *v) {
    unsigned int i;

    for (i = 0; i < K; ++i) {
        polyvecl_pointwise_acc_montgomery(&t->vec[i], &mat[i], v);
    }
}

/*************************************************
* Name:        polymatkm_pointwise_montgomery
*
* Description: Pointwise multiplies a k x m matrix with a vector of length m.
* The result is a vector of length k.
*
* Arguments:   - polyveck *t: pointer to the output vector
* - const polyvecm mat[K]: pointer to the input k x m matrix
* - const polyvecm *v: pointer to the input vector of length m
**************************************************/
void polymatkm_pointwise_montgomery(polyveck *t, const polyvecm mat[K],
                                    const polyvecm *v) {
    unsigned int i;

    for (i = 0; i < K; ++i) {
        polyvecm_pointwise_acc_montgomery(&t->vec[i], &mat[i], v);
    }
}

/*************************************************
* Name:        polymatkm_expand_ntt
*
* Description: Generates a k x m matrix A_gen with uniformly random coefficients
* and applies the NTT to each polynomial.
*
* Arguments:   - polyvecm mat[K]: pointer to the output k x m matrix in NTT domain
* - const uint8_t rho[]: pointer to the seed for matrix generation
**************************************************/
void polymatkm_expand_ntt(polyvecm mat[K], const uint8_t rho[SEEDBYTES]) {
    for (unsigned int i = 0; i < K; ++i) {
        for (unsigned int j = 0; j < M; ++j) {
            poly_uniform_ntt(&mat[i].vec[j], rho, (i << 8) | j);
        }
    }
}