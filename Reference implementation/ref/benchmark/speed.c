
#include <stdio.h>
#include <stdlib.h> 
#include <string.h> 
#include <time.h>   


#include "../include/params.h"
#include "../include/poly.h"
#include "../include/polyvec.h"
#include "../include/packing.h"
#include "../include/sign.h" 
#include "../include/sampler.h" 
#include "../include/ntt.h"     
#include "../include/polymat.h" 
#include "../include/symmetric.h" 


#include "../include/randombytes.h"
#include "cpucycles.h"
#include "speed_print.h"

#define NTESTS 1000 

uint64_t t[NTESTS]; 


static void randomize_poly(poly *p) {
    randombytes((uint8_t *)p->coeffs, N * sizeof(int32_t)); 
    
}

static void randomize_polyveck(polyveck *v) {
    for(int i=0; i<K; ++i) randomize_poly(&v->vec[i]);
}
static void randomize_polyvecl(polyvecl *v) {
    for(int i=0; i<L+K; ++i) randomize_poly(&v->vec[i]);
}
static void randomize_polyvecm(polyvecm *v) {
    for(int i=0; i<M; ++i) randomize_poly(&v->vec[i]);
}

int main() {
    
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES];
    uint8_t msg[SEEDBYTES * 2]; 
    size_t pklen=0, sklen=0, siglen=0;
    int i = 0;
    uint16_t nonce = 0;
    clock_t srt, ed; 
    long double time_avg;

    uint8_t seed[SEEDBYTES];
    uint8_t crh_seed[CRHBYTES]; 
    uint8_t key_seed[SEEDBYTES]; 
    uint8_t seed_y1[CRHBYTES];
    uint8_t seed_c_prime[CRHBYTES];
    uint8_t seed_x_prime[CRHBYTES];
    uint8_t seed_s_bar_prime_0[CRHBYTES]; 
    uint8_t seed_c_prime2[CRHBYTES]; 
    uint8_t seed_x_prime2[CRHBYTES]; 
    uint8_t seed_s_prime[CRHBYTES]; 

    polyveck b; 
    polyvecm sgen; 
    polyveck egen; 
    polyvecl s_prime; 
    poly s_bar_prime_0; 
    polyvecl z; 
    poly c; 
    poly x_poly; 
    polyvecl y; 
    polyveck temp_veck; 
    polyvecl temp_vecl; 
    poly test_poly_a, test_poly_b, test_poly_res; 
    polyvecm Agen[K]; 
    polyvecl s; 
    polyvecl s_ntt; 
    poly c_plus_x, c_plus_x_ntt; 
    polyvecl S_c_plus_x; 


    printf("Benchmarking CCC implementation (based on ccc_MODE=%d)\n", ccc_MODE);
    printf("Parameters: N=%d, Q=%d, K=%d, L=%d (L+K=%d, M=%d), TAU=%d, B_INFTY=%d\n", N, Q, K, L, L+K, M, TAU, B_INFTY);
    printf("Test iterations: %d\n", NTESTS);
    printf("Timing unit: CPU cycles (median/average)\n");
    printf("--------------------------------------------------\n");

    randombytes(seed, SEEDBYTES);
    randombytes(crh_seed, CRHBYTES);
    randombytes(key_seed, SEEDBYTES); 
    randombytes(msg, sizeof(msg)); 
    
    randombytes(seed_y1, CRHBYTES);
    randombytes(seed_c_prime, CRHBYTES);
    randombytes(seed_x_prime, CRHBYTES);
    randombytes(seed_s_bar_prime_0, CRHBYTES); 
    randombytes(seed_s_prime, CRHBYTES); 

    randombytes(seed_s_bar_prime_0, CRHBYTES); 
    randombytes(seed_c_prime2, CRHBYTES);
    randombytes(seed_x_prime2, CRHBYTES);
    

    // Keypair Generation
    srt = clock();
    for (i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        crypto_sign_keypair(pk, &pklen, sk, &sklen);
    }
    ed = clock();
    print_results("crypto_sign_keypair", t, NTESTS);
    time_avg = (long double)(ed - srt) * 1000.0 / CLOCKS_PER_SEC / NTESTS;
    printf("  Approx time: %.4Lf ms\n\n", time_avg);


    // Signature Generation
    crypto_sign_keypair(pk, &pklen, sk, &sklen);// Need a valid keypair first
    srt = clock();
    for (i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        // Assume crypto_sign_signature uses the modified SampleY internally correctly
        crypto_sign_signature(sig, &siglen, msg, sizeof(msg), sk);
    }
    ed = clock();
    print_results("crypto_sign_signature", t, NTESTS);
    time_avg = (long double)(ed - srt) * 1000.0 / CLOCKS_PER_SEC / NTESTS;
    printf("  Approx time: %.4Lf ms\n\n", time_avg);


    // Signature Verification
    // Ensure a valid signature exists for verification benchmark
    if (crypto_sign_signature(sig, &siglen, msg, sizeof(msg), sk) != 0) {
         fprintf(stderr, "Error: Failed to generate signature for verification benchmark setup.\n");
         // Handle error appropriately, maybe return 1
    }
    srt = clock();
    for (i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        crypto_sign_verify(sig, siglen, msg, sizeof(msg), pk, pklen);
    }
    ed = clock();
    print_results("crypto_sign_verify", t, NTESTS);
    time_avg = (long double)(ed - srt) * 1000.0 / CLOCKS_PER_SEC / NTESTS;
    printf("  Approx time: %.4Lf ms\n\n", time_avg);

    printf("--- Component Benchmarks ---\n");

    // Expand A_gen matrix (K x M)
    for (i = 0; i < NTESTS; ++i) {
        seed[0] = (uint8_t)i; // Vary seed slightly
        t[i] = cpucycles();
        polymatkm_expand(Agen, seed); // Uses poly_uniform
    }
    print_results("polymatkm_expand (A_gen)", t, NTESTS);

    // Sample secret component (s_gen or e_gen element)
    for (i = 0; i < NTESTS; ++i) {
        crh_seed[0] = (uint8_t)i;
        t[i] = cpucycles();
        SampleSigma1(&test_poly_a, crh_seed, (uint16_t)i);
    }
    print_results("SampleSigma1 (s_gen/e_gen element)", t, NTESTS);

    // NTT domain computation of b = Agen * sgen + egen
    // Prepare NTT inputs (outside timing loop)
    polyvecm sgen_ntt;
    polyveck egen_ntt; // This might be unused if keypair was modified - check sign.c
    polyvecm Agen_ntt[K];
    polyveck b_ntt;
    randomize_polyvecm(&sgen); // Random sgen
    randomize_polyveck(&egen); // Random egen
    for(int k_idx=0; k_idx<K; ++k_idx) randomize_polyvecm(&Agen[k_idx]); // Random Agen
    // Convert to NTT domain
    for(int k_idx=0; k_idx<K; ++k_idx) {
        Agen_ntt[k_idx] = Agen[k_idx]; // Copy
        polyvecm_ntt(&Agen_ntt[k_idx]);
    }
    sgen_ntt = sgen; // Copy
    polyvecm_ntt(&sgen_ntt);
    // If egen_ntt is truly unused in keypair, don't transform it.
    // Let's assume standard keypair uses NTT for egen add after INTT(A*s)
    // egen_ntt = egen; // Copy
    // polyveck_ntt(&egen_ntt); // Transform egen if needed by keypair logic

    for (i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        polymatkm_pointwise_montgomery(&b_ntt, Agen_ntt, &sgen_ntt); // b_ntt = NTT(A)*NTT(s)
        // Assuming standard keygen: INTT then add egen
        polyveck_invntt_tomont(&b_ntt); // b = INTT(NTT(A)*NTT(s))
        polyveck_add(&b_ntt, &b_ntt, &egen); // b = b + egen (non-NTT addition)
        polyveck_reduce(&b_ntt); // Reduce final b
    }
    print_results("Keypair Core NTT Calc (A*s + e)", t, NTESTS);

    // Sample challenge c
    for (i = 0; i < NTESTS; ++i) {
         crh_seed[0] = (uint8_t)i;
        t[i] = cpucycles();
        SampleChallenge(&c, crh_seed);
    }
    print_results("SampleChallenge", t, NTESTS);

    // Sample auxiliary x
    SampleChallenge(&c, crh_seed); // Need a challenge c first
    for (i = 0; i < NTESTS; ++i) {
        nonce = (uint16_t)i;
        crh_seed[1] = (uint8_t)i;
        t[i] = cpucycles();
        SampleX_poly(&x_poly, &c, crh_seed, nonce);
    }
    print_results("SampleX_poly", t, NTESTS);

    // **** Benchmark SampleY (with corrected call) ****
    srt = clock();
    for (i = 0; i < NTESTS; ++i) {
        nonce = (uint16_t)i; // Use i as nonce variation
        t[i] = cpucycles();
        // **** Corrected SampleY call with 9 arguments ****
        SampleY(&y, &s_prime, &s_bar_prime_0,
            seed_y1,
            seed_c_prime,
            seed_x_prime,
            seed_c_prime2,   
            seed_x_prime2,   
            nonce);
    }
    ed = clock();
    print_results("SampleY (3C Sampling - Corrected Call)", t, NTESTS);
    time_avg = (long double)(ed - srt) * 1000.0 / CLOCKS_PER_SEC / NTESTS;
    printf("  Approx time: %.4Lf ms\n\n", time_avg);
    // --- End SampleY Benchmark ---

    // NTT domain computation of s*(c+x)
    // Prepare inputs (outside timing loop)
    
    crypto_sign_keypair(pk, &pklen, sk, &sklen);
    
    uint8_t pk_tmp[CRYPTO_PUBLICKEYBYTES];
    size_t pklen_tmp;
    unpack_sk(pk_tmp, &pklen_tmp, &sgen, &egen, &s_prime, &s_bar_prime_0, key_seed, sk, sklen);

    // Reconstruct s
    poly_zero(&s.vec[0]); s.vec[0].coeffs[0] = 1;
    for(int m_idx=0; m_idx < M; ++m_idx) s.vec[1+m_idx] = sgen.vec[m_idx];
    for(int k_idx=0; k_idx < K; ++k_idx) s.vec[M+1+k_idx] = egen.vec[k_idx];
    s_ntt = s;
    polyvecl_ntt(&s_ntt); // s in NTT domain
    SampleChallenge(&c, crh_seed);
    SampleX_poly(&x_poly, &c, crh_seed, 0);
    poly_add(&c_plus_x, &c, &x_poly);
    poly_reduce(&c_plus_x); // Reduce c+x mod q
    c_plus_x_ntt = c_plus_x;
    poly_ntt(&c_plus_x_ntt); // c+x in NTT domain

    // Benchmark the pointwise multiplication loop only
     for (i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        // Assuming S_c_plus_x is used to store result
        for (int j = 0; j < L + K; ++j) { // Dimension L+K
            poly_pointwise_montgomery(&S_c_plus_x.vec[j], &s_ntt.vec[j], &c_plus_x_ntt);
            // Don't do INTT here, just benchmark the multiplication part
        }
     }
    print_results("poly_pointwise_montgomery loop (NTT(s)*NTT(c+x))", t, NTESTS);

    // Forward NTT
    randomize_poly(&test_poly_a);
    poly_reduce(&test_poly_a); // Ensure input is somewhat reasonable
    for (i = 0; i < NTESTS; ++i) {
        test_poly_b = test_poly_a; // Use a copy
        t[i] = cpucycles();
        poly_ntt(&test_poly_b); // Test NTT on copy
    }
    print_results("poly_ntt", t, NTESTS);

    // Inverse NTT
    poly_ntt(&test_poly_a); // Ensure input is in NTT form
    for (i = 0; i < NTESTS; ++i) {
        test_poly_b = test_poly_a; // Use a copy
        t[i] = cpucycles();
        poly_invntt_tomont(&test_poly_b); // Test INTT on copy
    }
    print_results("poly_invntt_tomont", t, NTESTS);

    // Pointwise Polynomial Multiplication (NTT domain)
    randomize_poly(&test_poly_a); poly_reduce(&test_poly_a);
    randomize_poly(&test_poly_b); poly_reduce(&test_poly_b);
    poly_ntt(&test_poly_a); // Prepare NTT operands
    poly_ntt(&test_poly_b);
    for (i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        poly_pointwise_montgomery(&test_poly_res, &test_poly_a, &test_poly_b);
    }
    print_results("poly_pointwise_montgomery (poly * poly)", t, NTESTS);

    // Polynomial Addition (Normal domain)
    randomize_poly(&test_poly_a); poly_reduce(&test_poly_a);
    randomize_poly(&test_poly_b); poly_reduce(&test_poly_b);
    for (i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        poly_add(&test_poly_res, &test_poly_a, &test_poly_b);
        // poly_reduce(&test_poly_res); // Reduce if needed, depends on context
    }
    print_results("poly_add", t, NTESTS);

    
    // Need valid pk, sk, sig first
    crypto_sign_keypair(pk, &pklen, sk, &sklen);
    if(crypto_sign_signature(sig, &siglen, msg, sizeof(msg), sk) != 0) {
         fprintf(stderr, "Error: Failed to generate sig for packing benchmark setup.\n");
         return 1;
    }
    // Need valid z (centered) and c for pack_sig test
    polyvecl z_unpacked;
    poly c_unpacked;
    if (unpack_sig(&z_unpacked, &c_unpacked, sig, siglen) != 0) {
        fprintf(stderr, "Error: unpack_sig failed during benchmark setup.\n");
        return 1;
    }

    // Pack Signature
    for (i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        pack_sig(sig, &z_unpacked, &c_unpacked);
    }
    print_results("pack_sig", t, NTESTS);

    // Unpack Signature
    for (i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        unpack_sig(&z, &c, sig, siglen); // Unpack into different vars if needed
    }
    print_results("unpack_sig", t, NTESTS);

    // Pack Public Key
    unpack_pk(seed, &b, pk, pklen);// Get components for packing
    for (i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        pack_pk(pk, seed, &b);
    }
    print_results("pack_pk", t, NTESTS);
    
    polyveck Col0;
    { 
        polyveck qj_vec;
        polyveck two_b = b; 
        poly_zero(&qj_vec.vec[0]); 
        qj_vec.vec[0].coeffs[0] = Q;
        for(int i = 1; i < K; ++i) { 
            poly_zero(&qj_vec.vec[i]); 
        }
        polyveck_mod_2q(&qj_vec);
        polyveck_double(&two_b);
        polyveck_sub(&Col0, &qj_vec, &two_b);
        polyveck_mod_2q(&Col0);
    }
    
    // Unpack Public Key
     for (i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        unpack_pk(seed, &b, pk, pklen);
    }
    print_results("unpack_pk", t, NTESTS);

    // Pack Secret Key
    unpack_sk(pk_tmp, &pklen_tmp, &sgen, &egen, &s_prime, &s_bar_prime_0, key_seed, sk, sklen); 
    for (i = 0; i < NTESTS; ++i) {
       t[i] = cpucycles();
        pack_sk(sk, pk, pklen, &sgen, &egen, &s_prime, &s_bar_prime_0, key_seed);
     }
    print_results("pack_sk", t, NTESTS);

    // Unpack Secret Key
     for (i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
         unpack_sk(pk_tmp, &pklen_tmp, &sgen, &egen, &s_prime, &s_bar_prime_0, key_seed, sk, sklen);
    }
    print_results("unpack_sk", t, NTESTS);
    
    polyveck w_recomputed; 

    if (unpack_sig(&z_unpacked, &c_unpacked, sig, siglen) != 0) {
        fprintf(stderr, "Error: unpack_sig failed while preparing inputs for calculate_A_vec_prod_crt benchmark.\n");
        return 1;
    }

    for (i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        
        calculate_A_vec_prod_crt(&w_recomputed, Agen, &Col0, &z_unpacked);
    }
    print_results("calculate_A_vec_prod_crt", t, NTESTS);
    printf("\n--------------------------------------------------\n");
    printf("Benchmarking complete.\n");

    return 0;
}