#include "sign.h"
#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polymat.h"
#include "polyvec.h"
#include "randombytes.h"
#include "symmetric.h"
#include "ntt.h" 
#include "sampler.h"
#include "reduce.h"
#include <stdlib.h> 
#include <inttypes.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>  
#include <limits.h> 
//unsigned long long global_sign_cycles_count = 0;

/*************************************************
* Name:        center_coeff
*
* Description: Centers a coefficient modulo q to the range [-q/2, q/2].
*
* Arguments:   - int32_t coeff: input coefficient
*
* Returns:     Centered coefficient.
**************************************************/
static inline int32_t center_coeff(int32_t coeff) {
    // Assumes coeff is already reduced mod q, e.g., [0, Q-1]
    int32_t c = coeff;
    // Use freeze first to ensure it's in [0, Q-1] range robustly
    c = freeze(coeff); // Ensure input is [0, Q-1] before centering
    if (c > Q / 2) {
        c -= Q;
    }
    // Now c is in [-Q/2, Q/2] (or close, depending on Q odd/even)
    return c;
}


/*************************************************
* Name:        crypto_sign_keypair
*
* Description: Generates a public and private key pair for the Rhyme signature scheme.
* This function implements Algorithm 8 from the Rhyme paper.
*
* Arguments:   - uint8_t *pk: pointer to output public key byte array
* - size_t *pklen: pointer to output length of public key
* - uint8_t *sk: pointer to output private key byte array
* - size_t *sklen: pointer to output length of private key
*
* Returns:     0 on success.
**************************************************/
int crypto_sign_keypair(uint8_t *pk, size_t *pklen, uint8_t *sk, size_t *sklen) {
    uint8_t seedbuf[SEEDBYTES + CRHBYTES + CRHBYTES + CRHBYTES + SEEDBYTES]; 
    const uint8_t *seedA, *seedsk, *seed_s_prime, *seed_s_bar_prime, *key; 

    polyvecm Agen_ntt[K];
    polyveck b_ntt, ntt_neg_two_b;
    polyvecm sgen, sgen_ntt;
    polyveck egen, egen_ntt;
    polyvecl s_prime;
    poly s_bar_prime_0; 
    
    randombytes(seedbuf, SEEDBYTES);
    shake256(seedbuf, sizeof(seedbuf), seedbuf, SEEDBYTES);
    seedA = seedbuf;
    seedsk = seedA + SEEDBYTES;
    seed_s_prime = seedsk + CRHBYTES;
    seed_s_bar_prime = seed_s_prime + CRHBYTES;
    key = seed_s_bar_prime + CRHBYTES;

    polymatkm_expand_ntt(Agen_ntt, seedA);

    for(int i=0; i < M; ++i) SampleSigma1(&sgen.vec[i], seedsk, (uint16_t)i);
    for(int i=0; i < K; ++i) SampleSigma1(&egen.vec[i], seedsk, (uint16_t)(M+i));
    for(int i = 0; i < L + K; ++i) {
        SampleSigma1(&s_prime.vec[i], seed_s_prime, (uint16_t)i);
    }
 
    SampleSigma1(&s_bar_prime_0, seed_s_bar_prime, 0);

    sgen_ntt = sgen;
    egen_ntt = egen;
    polyvecm_ntt(&sgen_ntt);
    polyveck_ntt(&egen_ntt);

    polymatkm_pointwise_montgomery(&b_ntt, Agen_ntt, &sgen_ntt);
    polyveck_add(&b_ntt, &b_ntt, &egen_ntt);

    polyveck_negate(&ntt_neg_two_b, &b_ntt);
    polyveck_freeze(&ntt_neg_two_b);
    *pklen = pack_pk(pk, seedA, &ntt_neg_two_b);
    *sklen = pack_sk(sk, pk, *pklen, &sgen, &egen, &s_prime, &s_bar_prime_0, key);

    return 0;
}


/*************************************************
* Name:        crypto_sign_signature
*
* Description: Computes a signature for a given message.
* This function implements Algorithm 9 from the Rhyme paper,
* using an iterative process to find a signature that
* satisfies the norm bound.
*
* Arguments:   - uint8_t *sig: pointer to output signature byte array
* - size_t *siglen: pointer to output length of signature
* - const uint8_t *m: pointer to the message to be signed
* - size_t mlen: length of the message
* - const uint8_t *sk: pointer to the private key
*
* Returns:     0 if a signature was successfully generated, -1 otherwise.
**************************************************/
int crypto_sign_signature(uint8_t *sig, size_t *siglen, const uint8_t *m,
                          size_t mlen, const uint8_t *sk)
{

    uint8_t key[SEEDBYTES];
    uint8_t mu[CRHBYTES];
    uint8_t c_seed[CRHBYTES];
    uint8_t seedA[SEEDBYTES];
    uint8_t pk_temp[CRYPTO_PUBLICKEYBYTES];
    size_t pklen_temp;
    xof256_state state;
    
    uint8_t seedysc_base[CRHBYTES];
    uint8_t seedx_base[CRHBYTES];
    uint8_t iter_seeds_buf[5 * CRHBYTES];
    uint8_t current_seedx[CRHBYTES];
    uint16_t count = 0;
    uint8_t count_bytes[2];

    poly y0;
    polyvecd_rest y_rest;
    poly z0;
    polyvecd_rest z_rest;
    stream256_state rng_state;
    int loop_done;

    polyvecm sgen;
    polyveck egen;
    polyvecl s;
    polyvecl s_prime;
    poly c;
    poly x_poly;
    poly c_plus_x;
    polyveck w, w_for_hash;
    uint8_t packed_w[K * POLY_PACKEDBYTES_2Q];
    poly s_bar_prime_0; 
    
    polyvecm Agen_ntt[K];
    polyveck ntt_neg_two_b, ntt_Col0;
    
    unpack_sk(pk_temp, &pklen_temp, &sgen, &egen, &s_prime, &s_bar_prime_0, key, sk, CRYPTO_SECRETKEYBYTES);    
    unpack_pk(seedA, &ntt_neg_two_b, pk_temp, pklen_temp);

    polymatkm_expand_ntt(Agen_ntt, seedA);

    ntt_Col0 = ntt_neg_two_b; 


    poly_zero(&s.vec[0]); s.vec[0].coeffs[0] = 1;
    for(int i=0; i < M; ++i) s.vec[1+i] = sgen.vec[i];
    for(int i=0; i < K; ++i) s.vec[M+1+i] = egen.vec[i];


#if ccc_MODE != 2
    polyvecl s_ntt;
    s_ntt = s;
    polyvecl_ntt(&s_ntt); 
#endif

    shake256_init(&state);
    shake256_absorb(&state, pk_temp, pklen_temp);
    shake256_absorb(&state, m, mlen);
    shake256_finalize(&state);
    shake256_squeeze(mu, CRHBYTES, &state);

    uint8_t base_seeds_buf[CRHBYTES * 2];
    shake256_init(&state);
    shake256_absorb(&state, key, SEEDBYTES);
    shake256_absorb(&state, mu, CRHBYTES);
    shake256_finalize(&state);
    shake256_squeeze(base_seeds_buf, sizeof(base_seeds_buf), &state);
    memcpy(seedysc_base, base_seeds_buf, CRHBYTES);
    memcpy(seedx_base, base_seeds_buf + CRHBYTES, CRHBYTES);

    shake256_init(&state);
    shake256_absorb(&state, seedx_base, CRHBYTES);
    count_bytes[0] = 0;
    count_bytes[1] = 0;
    shake256_absorb(&state, count_bytes, 2);
    shake256_finalize(&state);
    shake256_squeeze(current_seedx, CRHBYTES, &state);

    loop_done = 0;
    while(!loop_done) {
        count++;

        //global_sign_cycles_count++;
        
        //if (count > MAX_SIGN_ITERATIONS) return -1;
        
        count_bytes[0] = count & 0xFF;
        count_bytes[1] = (count >> 8) & 0xFF;

        shake256_init(&state);
        shake256_absorb(&state, seedysc_base, CRHBYTES);
        shake256_absorb(&state, count_bytes, 2);
        shake256_finalize(&state);
        shake256_squeeze(iter_seeds_buf, sizeof(iter_seeds_buf), &state);

        /* Sample y */
        SampleY0(&y0, iter_seeds_buf, count); 
        SampleY_rest(&y_rest, iter_seeds_buf + CRHBYTES, count); 

        /* Compute w = Ay */

        polyvecl y_combined;
        y_combined.vec[0] = y0; 
        for(int k=0; k<D_REST; k++) {
            y_combined.vec[k+1] = y_rest.vec[k]; 
        }

        calculate_A_vec_prod_crt(&w, Agen_ntt, &ntt_Col0, &y_combined);
        
        polyveck_freeze2q(&w);

        w_for_hash = w; 
        pack_polyveck_2q(packed_w, &w_for_hash);

        shake256_init(&state);
        shake256_absorb(&state, packed_w, K * POLY_PACKEDBYTES_2Q);
        shake256_absorb(&state, mu, CRHBYTES);
        shake256_finalize(&state);
        shake256_squeeze(c_seed, CRHBYTES, &state);
        
        SampleChallenge(&c, c_seed);
        SampleX_poly(&x_poly, &c, current_seedx, count);

        stream256_init(&rng_state, seedysc_base, count);
        

        if (!Sample_Z(&z0, &x_poly, &c, &y0, &rng_state)) {
            continue; 
        }

        poly_add(&c_plus_x, &c, &x_poly);
        int z_rest_fail = 0;


#if ccc_MODE == 2
        /* Mode 2: Sparse multiplication */
        for(int i=0; i < D_REST; ++i) {
            poly temp;
            poly_mul_sparse_ternary(&temp, &s.vec[i+1], &c, &x_poly);
            poly_add(&z_rest.vec[i], &y_rest.vec[i], &temp);
            if(poly_check_norm_inf(&z_rest.vec[i], B1)) { z_rest_fail = 1; break; }
        }
#else
        /* Mode 3/5: NTT-based multiplication */
        poly c_plus_x_ntt = c_plus_x;
        poly_ntt(&c_plus_x_ntt);

        for(int i=0; i < D_REST; ++i) {
            poly temp;
            
            poly_pointwise_montgomery(&temp, &s_ntt.vec[i+1], &c_plus_x_ntt);
            poly_invntt_tomont(&temp);
            
            /* Add and center coefficients */
            for(int k=0; k<N; k++) {
                int16_t val = temp.coeffs[k];
                /* Fast reduction mod 257 */
                int16_t t = (val & 0xFF) - (val >> 8);
                
                if (t > 128) t -= 257;
                else if (t < -128) t += 257;
                
                z_rest.vec[i].coeffs[k] = y_rest.vec[i].coeffs[k] + t;
            }
            
            if(poly_check_norm_inf(&z_rest.vec[i], B1)) { 
                z_rest_fail = 1; 
                break; 
            }
        }
#endif
        if (z_rest_fail) {
            continue;
        }

        loop_done = 1;
    }

    *siglen = pack_sig(sig, &z0, &z_rest, &c);
    if (*siglen == 0) return -1; 

    return 0;
}


/*************************************************
* Name:        crypto_sign_verify
*
* Description: Verifies a signature against a public key and message.
* This function implements Algorithm 10 from the Rhyme paper.
*
* Arguments:   - const uint8_t *sig: pointer to the signature to be verified
* - size_t siglen: length of the signature
* - const uint8_t *m: pointer to the message
* - size_t mlen: length of the message
* - const uint8_t *pk: pointer to the public key
* - size_t pklen: length of the public key
*
* Returns:     0 if the signature is valid, -1 otherwise.
**************************************************/


int crypto_sign_verify(const uint8_t *sig, size_t siglen, const uint8_t *m, size_t mlen,
                       const uint8_t *pk, size_t pklen)
{
    uint8_t seedA[SEEDBYTES];
    uint8_t mu[CRHBYTES];
    uint8_t cprime_seed[CRHBYTES];
    
    poly z0;
    polyvecd_rest z_rest; 
    poly c, cprime;
    
    polyvecm Agen_ntt[K];
    polyveck ntt_neg_two_b, ntt_Col0;
    polyveck Az_mod_2q;
    polyveck w_recomputed, w_hash_input;
    polyveck qcj; 
    uint8_t packed_w_recomputed[K * POLY_PACKEDBYTES_2Q];
    
    xof256_state state;

    if (unpack_pk(seedA, &ntt_neg_two_b, pk, pklen) != 0) {
         fprintf(stderr, "Error: unpack_pk\n");
        return -1;
    }

    if (unpack_sig(&z0, &z_rest, &c, sig, siglen) != 0) {
        fprintf(stderr, "Error: unpack_sig\n");
        return -1;
    }

    // check z0 <= B0
    if (poly_check_norm_inf(&z0, B0)) {
        // fprintf(stderr, "Error: z0 norm > B0\n");
        return -1;
    }
    // check z_rest <= B1
    for(int i=0; i<D_REST; i++) {
        if (poly_check_norm_inf(&z_rest.vec[i], B1)) {
            // fprintf(stderr, "Error: z_rest[%d] norm > B1\n", i);
            return -1;
        }
    }

    polymatkm_expand_ntt(Agen_ntt, seedA);

    ntt_Col0 = ntt_neg_two_b; 

    polyvecl z_combined;
    z_combined.vec[0] = z0;
    for(int k=0; k<D_REST; k++) {
        z_combined.vec[k+1] = z_rest.vec[k];
    }


    calculate_A_vec_prod_crt(&Az_mod_2q, Agen_ntt, &ntt_Col0, &z_combined);

    /* Compute w = Az - q*c*j (mod 2q) */

    poly_mul_const(&qcj.vec[0], &c, Q); // qcj[0] = c * Q
    for(int i = 1; i < K; ++i) {
        poly_zero(&qcj.vec[i]);     // qcj[1..K-1] = 0
    }
    
    polyveck_sub(&w_recomputed, &Az_mod_2q, &qcj);
    polyveck_reduce2q(&w_recomputed); 

    // Hash verification
    shake256_init(&state);
    shake256_absorb(&state, pk, pklen);
    shake256_absorb(&state, m, mlen);
    shake256_finalize(&state);
    shake256_squeeze(mu, CRHBYTES, &state);

    w_hash_input = w_recomputed;
    polyveck_freeze2q(&w_hash_input); 
    pack_polyveck_2q(packed_w_recomputed, &w_hash_input);

    shake256_init(&state);
    shake256_absorb(&state, packed_w_recomputed, K * POLY_PACKEDBYTES_2Q);
    shake256_absorb(&state, mu, CRHBYTES);
    shake256_finalize(&state);
    shake256_squeeze(cprime_seed, CRHBYTES, &state);

    SampleChallenge(&cprime, cprime_seed);

    if (poly_compare(&c, &cprime) != 0) {
        return -1; 
    }

    return 0; 
}



/*************************************************
* Name:        crypto_sign_sign
*
* Description: Standard API wrapper for creating a signed message.
* The output format is | 2-byte siglen | message | signature |.
*
* Arguments:   - uint8_t *sm: pointer to output signed message
* - size_t *smlen: pointer to output length of signed message
* - const uint8_t *m: pointer to the message to be signed
* - size_t mlen: length of the message
* - const uint8_t *sk: pointer to the private key
*
* Returns:     0 on success, -1 on failure.
**************************************************/
int crypto_sign_sign(uint8_t *sm, size_t *smlen, const uint8_t *m, size_t mlen,
                     const uint8_t *sk)
{
    size_t siglen; 
    int ret;
    const size_t siglen_field_bytes = sizeof(uint16_t); 

    if (!smlen) return -1;
    *smlen = 0; 

    uint8_t *sig_ptr = sm + siglen_field_bytes + mlen; 
    ret = crypto_sign_signature(sig_ptr, &siglen, m, mlen, sk);

    if (ret != 0) {
        
        return -1; 
    }
    if (siglen > UINT16_MAX) {
        fprintf(stderr, "Error: Actual signature length %zu exceeds the maximum value of uint16_t\n", siglen);
        return -1; 
    }

    size_t final_smlen = siglen_field_bytes + mlen + siglen;
    *smlen = final_smlen;

    sm[0] = (uint8_t)(siglen & 0xFF);
    sm[1] = (uint8_t)((siglen >> 8) & 0xFF);

    memmove(sm + siglen_field_bytes, m, mlen);

    return 0; 
}


/*************************************************
* Name:        crypto_sign_open
*
* Description: Standard API wrapper for verifying a signed message and
* recovering the original message.
*
* Arguments:   - uint8_t *m: pointer to output message
* - size_t *mlen: pointer to output length of message
* - const uint8_t *sm: pointer to the signed message
* - size_t smlen: length of the signed message
* - const uint8_t *pk: pointer to the public key
*
* Returns:     0 on success, -1 on failure.
**************************************************/
 int crypto_sign_open(uint8_t *m, size_t *mlen, const uint8_t *sm, size_t smlen,
                     const uint8_t *pk, const size_t pklen)
{
    const size_t siglen_field_bytes = sizeof(uint16_t); 
    size_t actual_siglen;
    size_t calculated_mlen;

    *mlen = 0;

    if (smlen < siglen_field_bytes) {
        fprintf(stderr, "Error (crypto_sign_open): smlen (%zu) is too small to read signature length field (%zu)\n", smlen, siglen_field_bytes);
        return -1;
    }

    actual_siglen = (uint16_t)sm[0] | ((uint16_t)sm[1] << 8);

    if (smlen < siglen_field_bytes + actual_siglen) {
        fprintf(stderr, "Error (crypto_sign_open): smlen (%zu) is less than the signature length field (%zu) + the read signature length (%zu)\n", smlen, siglen_field_bytes, actual_siglen);
        return -1; 
    }

    calculated_mlen = smlen - siglen_field_bytes - actual_siglen;

    const uint8_t *message_ptr = sm + siglen_field_bytes;
    const uint8_t *sig_ptr = sm + siglen_field_bytes + calculated_mlen; 

    if (crypto_sign_verify(sig_ptr, actual_siglen, message_ptr, calculated_mlen, pk, pklen) != 0) {
        return -1; 
    }

    *mlen = calculated_mlen;
    memmove(m, message_ptr, calculated_mlen);

    return 0; 
}