#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "reduce.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "encoding.h" 


/*************************************************
* Name:        center_for_packing
*
* Description: Centers a value modulo Q to the range [-(Q-1)/2, (Q-1)/2].
* Used before packing to ensure correct sign representation.
*
* Arguments:   - int32_t v: input value
*
* Returns:     Centered value.
**************************************************/
static int32_t center_for_packing(int32_t v) {
    v = v % Q;
    if (v < 0) v += Q;
    if (v > Q / 2) { 
        v -= Q;      
    }
    return v;
}
/*************************************************
* Name:        pack_pk
*
* Description: Packs the public key into a byte array without a length field.
* The format is | seedA | fixed-size main data | overflow data |.
*
* Arguments:   - uint8_t *pk: pointer to output byte array
* - const uint8_t seedA[]: pointer to seed for matrix A
* - const polyveck *b: pointer to polynomial vector b
*
* Returns:     The total number of bytes written to pk.
**************************************************/
size_t pack_pk(uint8_t *pk, const uint8_t seedA[SEEDBYTES], const polyveck *b) {
    uint8_t *pk_ptr = pk;
    
    // Copy seedA
    memcpy(pk_ptr, seedA, SEEDBYTES);
    pk_ptr += SEEDBYTES;

    //  Prepare pointers for fixed-size main data and temporary overflow buffer
    uint8_t *main_buf_ptr = pk_ptr;
    uint8_t overflow_buf[(K * N + 7) / 8] = {0};
    int total_overflow_bits = 0;

    //  Process all polynomials into a fixed-size main area and a variable overflow area
    for (int i = 0; i < K; ++i) {
        const poly *p = &b->vec[i];
        for (int j = 0; j < N; ++j) {
            int32_t coeff = p->coeffs[j];
            if (coeff < 255) {
                main_buf_ptr[i * N + j] = (uint8_t)coeff;
            } else {
                main_buf_ptr[i * N + j] = 0xFF; // Mark as overflow
                if (coeff == 256) {
                    overflow_buf[total_overflow_bits / 8] |= (1 << (total_overflow_bits % 8));
                }
                total_overflow_bits++;
            }
        }
    }

    //  Advance pointer past the fixed-size main data area
    pk_ptr += K * N;

    //  Append the overflow data
    size_t overflow_byte_len = (total_overflow_bits + 7) / 8;
    memcpy(pk_ptr, overflow_buf, overflow_byte_len);
    pk_ptr += overflow_byte_len;

    //  Return the total actual size
    return (size_t)(pk_ptr - pk);
}

/*************************************************
* Name:        unpack_pk
*
* Description: Unpacks a public key from a byte array that has no length field.
*
* Arguments:   - uint8_t seedA[]: pointer to output seed for matrix A
* - polyveck *b: pointer to output polynomial vector b
* - const uint8_t *pk: pointer to input byte array
* - size_t pklen: length of the input byte array
*
* Returns:     0 on success, non-zero on failure.
**************************************************/
int unpack_pk(uint8_t seedA[SEEDBYTES], polyveck *b, const uint8_t *pk, size_t pklen) {
    if (pklen < SEEDBYTES + (size_t)K * N) return -1; // Basic check for minimal possible size

    const uint8_t *pk_ptr = pk;

    // 1. Unpack seedA
    memcpy(seedA, pk_ptr, SEEDBYTES);
    pk_ptr += SEEDBYTES;

    // 2. Scan the fixed-size main data area to determine overflow length
    const uint8_t *main_buf_ptr = pk_ptr;
    int total_overflow_bits = 0;
    for (size_t i = 0; i < (size_t)K * N; ++i) {
        if (main_buf_ptr[i] == 0xFF) {
            total_overflow_bits++;
        }
    }
    
    size_t overflow_byte_len = (total_overflow_bits + 7) / 8;
    
    // 3. Verify consistency of the provided pklen
    if (pklen != SEEDBYTES + (size_t)K * N + overflow_byte_len) {
        return -2; // Length mismatch
    }

    // 4. Set pointer to the start of the overflow data
    const uint8_t *overflow_buf_ptr = main_buf_ptr + K * N;
    int overflow_bit_pos = 0;

    // 5. Reconstruct the polynomial vector b
    for (int i = 0; i < K; ++i) {
        poly *p = &b->vec[i];
        for (int j = 0; j < N; ++j) {
            uint8_t byte = main_buf_ptr[i * N + j];
            if (byte == 0xFF) {
                if ((size_t)(overflow_bit_pos / 8) >= overflow_byte_len) return -3; // Malformed data
                
                int bit = (overflow_buf_ptr[overflow_bit_pos / 8] >> (overflow_bit_pos % 8)) & 1;
                overflow_bit_pos++;
                p->coeffs[j] = 255 + bit;
            } else {
                p->coeffs[j] = byte;
            }
        }
    }

    return 0; 
}

/*************************************************
* Name:        pack_sk
*
* Description: Packs the secret key.
**************************************************/
size_t pack_sk(uint8_t *sk, const uint8_t *pk, size_t pklen, const polyvecm *s_gen, const polyveck *e_gen, const polyvecl *s_prime, const poly *s_bar_prime_0, const uint8_t key[SEEDBYTES]) {
    uint8_t *sk_ptr = sk;
    const size_t sgen_bytes = M * SK_POLY_PACKEDBYTES;
    const size_t egen_bytes = K * SK_POLY_PACKEDBYTES;
    const size_t sprime_bytes = S_PRIME_PACKEDBYTES;
    const size_t s_bar_prime_bytes = S_BAR_PRIME_0_PACKEDBYTES; 

    memcpy(sk_ptr, pk, pklen);
    sk_ptr += pklen;
    pack_polyvecm_sigma1(sk_ptr, s_gen);
    sk_ptr += sgen_bytes;
    pack_polyveck_sigma1(sk_ptr, e_gen);
    sk_ptr += egen_bytes;


    pack_polyvecl_s_prime(sk_ptr, s_prime); 
    sk_ptr += sprime_bytes;
    {
        size_t idx = 0;
        uint8_t byte = 0;
        int bit = 0;
        
        for (int j = 0; j < N; ++j) {
            
            byte |= ((s_bar_prime_0->coeffs[j] + 1) >> 1) << bit;
            if (++bit == 8) {
                if (idx >= s_bar_prime_bytes) break; 
                sk_ptr[idx++] = byte;
                byte = 0;
                bit = 0;
            }
        }
            if (bit > 0 && idx < s_bar_prime_bytes) {
                sk_ptr[idx] = byte;
            }
    }

    sk_ptr += s_bar_prime_bytes;
    memcpy(sk_ptr, key, SEEDBYTES);
    sk_ptr += SEEDBYTES;
    return (size_t)(sk_ptr - sk);
}

/*************************************************
* Name:        unpack_sk
*
* Description: Unpacks a secret key.
**************************************************/
int unpack_sk(uint8_t *pk, size_t *pklen, polyvecm *s_gen, polyveck *e_gen, polyvecl *s_prime, poly *s_bar_prime_0, uint8_t key[SEEDBYTES], const uint8_t *sk, size_t sklen) {
    const uint8_t *sk_ptr = sk;
    const size_t sgen_bytes = M * POLYSIGMA1_PACKEDBYTES;
    const size_t egen_bytes = K * POLYSIGMA1_PACKEDBYTES;
    const size_t sprime_bytes = S_PRIME_PACKEDBYTES; 
    const size_t s_bar_prime_bytes = S_BAR_PRIME_0_PACKEDBYTES;

    // Deduce the length of the embedded public key
    if (sklen < SEEDBYTES + (size_t)K * N) return -1; // Not enough data for even a minimal PK
    
    const uint8_t *main_buf_ptr = sk_ptr + SEEDBYTES;
    int total_overflow_bits = 0;
    for (size_t i = 0; i < (size_t)K * N; ++i) {
        if (main_buf_ptr[i] == 0xFF) {
            total_overflow_bits++;
        }
    }
    size_t overflow_byte_len = (total_overflow_bits + 7) / 8;
    size_t embedded_pklen = SEEDBYTES + (size_t)K * N + overflow_byte_len;

    if (sklen < embedded_pklen) return -1; // sklen too short for the deduced pklen
    *pklen = embedded_pklen;
    memcpy(pk, sk_ptr, *pklen);
    sk_ptr += *pklen;

    if ((size_t)(sk_ptr - sk) + sgen_bytes + egen_bytes + sprime_bytes + SEEDBYTES > sklen) return -1;
    unpack_polyvecm_sigma1(s_gen, sk_ptr);
    sk_ptr += sgen_bytes;
    unpack_polyveck_sigma1(e_gen, sk_ptr);
    sk_ptr += egen_bytes;


    unpack_polyvecl_s_prime(s_prime, sk_ptr); 
    sk_ptr += sprime_bytes;
    if ((size_t)(sk_ptr - sk) + s_bar_prime_bytes > sklen) return -1;
    {
        size_t idx = 0;
        uint8_t byte = sk_ptr[0];
        int bit = 0;
        for (int j = 0; j < N; ++j) {
            if (bit == 8) {
                if (++idx >= s_bar_prime_bytes && j < N - 1) break; 
                byte = idx < s_bar_prime_bytes ? sk_ptr[idx] : 0;
                bit = 0;
            }
            
            s_bar_prime_0->coeffs[j] = ((byte >> bit) & 1) * 2 - 1;
            bit++;
        }
    }
    sk_ptr += s_bar_prime_bytes;
    
    if ((size_t)(sk_ptr - sk) + SEEDBYTES > sklen) return -1;
    memcpy(key, sk_ptr, SEEDBYTES);

    return 0;
}

/*************************************************
* Name:        pack_sig
*
* Description: Packs the signature (vector z and polynomial c) into a byte array.
* Uses rANS for compressing vector z.
* The format is | 2-byte z_len | rANS encoded z | packed c |.
*
* Arguments:   - uint8_t sig[]: pointer to output signature byte array
* - const polyvecl *z: pointer to input vector z
* - const poly *c: pointer to input challenge polynomial
*
* Returns:     The total number of bytes written to sig, or 0 on failure.
**************************************************/
size_t pack_sig(uint8_t sig[CRYPTO_BYTES], const poly *z0, const polyvecd_rest *z_rest, const poly *c) {
    uint8_t *sig_ptr = sig;
    uint16_t size_encoded_z;

    /* Combine z0 and z_rest into a single vector for rANS encoding */
    polyvecl z_combined;
    
    /* Center coefficients for z0 */
    for(int j=0; j<N; ++j) {
        z_combined.vec[0].coeffs[j] = center_for_packing(z0->coeffs[j]);
    }
    
    /* Center coefficients for z_rest */
    for(int i = 0; i < D_REST; ++i) {
        for(int j=0; j<N; ++j) {
            z_combined.vec[i+1].coeffs[j] = center_for_packing(z_rest->vec[i].coeffs[j]);
        }
    }

    /* Perform rANS encoding */
    size_encoded_z = encode_z(sig_ptr + CCC_RANS_SIZE_BYTES, &z_combined);
    if (size_encoded_z == 0) {
        fprintf(stderr, "[MODE %d] Error: encode_z failed in pack_sig.\n", ccc_MODE);
        return 0; 
    }

    /* Compute size and write */
    size_t total_needed = CCC_RANS_SIZE_BYTES + size_encoded_z + POLYC_PACKEDBYTES;
    if (total_needed > CRYPTO_BYTES) {
        return 0; 
    }

    sig_ptr[0] = (uint8_t)(size_encoded_z & 0xFF);
    sig_ptr[1] = (uint8_t)((size_encoded_z >> 8) & 0xFF);
    sig_ptr += CCC_RANS_SIZE_BYTES;
    sig_ptr += size_encoded_z;

    /*  Pack challenge c */
    pack_poly_challenge(sig_ptr, c);
    
    return (size_t)(sig_ptr + POLYC_PACKEDBYTES - sig);
}



/*************************************************
* Name:        unpack_sig
*
* Description: Unpacks a signature from a byte array into vector z and polynomial c.
*
* Arguments:   - polyvecl *z: pointer to output vector z
* - poly *c: pointer to output challenge polynomial
* - const uint8_t *sig: pointer to input signature byte array
* - size_t siglen: length of the input signature
*
* Returns:     0 on success, non-zero on failure.
**************************************************/
int unpack_sig(poly *z0, polyvecd_rest *z_rest, poly *c, const uint8_t *sig, size_t siglen) {
    const uint8_t *sig_ptr = sig;
    uint16_t size_encoded_z;
    size_t min_required_len = CCC_RANS_SIZE_BYTES + POLYC_PACKEDBYTES;

    if (siglen < min_required_len) return -1;

    /*  Read compressed data length */
    size_encoded_z = (uint16_t)sig_ptr[0] | ((uint16_t)sig_ptr[1] << 8);
    sig_ptr += CCC_RANS_SIZE_BYTES;

    /* Verify length */
    size_t expected_total_len = CCC_RANS_SIZE_BYTES + size_encoded_z + POLYC_PACKEDBYTES;
    if (siglen < expected_total_len) {
    fprintf(stderr, "Buffer too small: siglen %zu < expected %zu\n", siglen, expected_total_len);
    return -2;
}
    
    /* Decode into a temporary combined vector */
    polyvecl z_combined;
    int decode_status = decode_z(&z_combined, sig_ptr, size_encoded_z);
    if (decode_status != 0) return -3; 
    
    sig_ptr += size_encoded_z;

    /* Distribute decoded data to z0 and z_rest, correcting centered values */
    /* Process z0 */
    for(int j=0; j<N; ++j) {
        int32_t val = z_combined.vec[0].coeffs[j];
        if (val > Q/2) {
            val -= Q;
        }
        z0->coeffs[j] = val;
    }
    
    /* Process z_rest */
    for(int i = 0; i < D_REST; ++i) {
        for(int j=0; j<N; ++j) {
            int32_t val = z_combined.vec[i+1].coeffs[j];
            if (val > Q/2) {
                val -= Q;
            }
            z_rest->vec[i].coeffs[j] = val;
        }
    }
    unpack_poly_challenge(c, sig_ptr);
    
    return 0;
}

/*************************************************
* Name:        pack_polyvecl_sigma1
*
* Description: Packs a polynomial vector of length L+K with coefficients in {-1, 1}.
* Each coefficient is represented by a single bit.
*
* Arguments:   - uint8_t *buf: pointer to the output byte array
* - const polyvecl *s: pointer to the input polynomial vector
**************************************************/
void pack_polyvecl_sigma1(uint8_t *buf, const polyvecl *s) {
    size_t idx = 0; uint8_t byte = 0; int bit = 0;
    const size_t total_bytes = (K + L) * POLYSIGMA1_PACKEDBYTES;
    for (int i = 0; i < K + L; ++i) for (int j = 0; j < N; ++j) {
        byte |= ((s->vec[i].coeffs[j] + 1) >> 1) << bit;
        if (++bit == 8) { if(idx>=total_bytes) return; buf[idx++] = byte; byte = 0; bit = 0; }
    }
    if (bit > 0 && idx < total_bytes) buf[idx++] = byte;
}

/*************************************************
* Name:        unpack_polyvecl_sigma1
*
* Description: Unpacks a polynomial vector of length L+K with coefficients in {-1, 1}.
*
* Arguments:   - polyvecl *s: pointer to the output polynomial vector
* - const uint8_t *buf: pointer to the input byte array
**************************************************/
void unpack_polyvecl_sigma1(polyvecl *s, const uint8_t *buf) {
     size_t idx = 0; uint8_t byte = buf[0]; int bit = 0;
     const size_t total_bytes = (K + L) * POLYSIGMA1_PACKEDBYTES;
     for (int i = 0; i < K + L; ++i) for (int j = 0; j < N; ++j) {
          if (bit == 8) { if(++idx >= total_bytes && (i<K+L-1 || j<N-1)) return; byte = idx < total_bytes ? buf[idx] : 0; bit = 0; }
          s->vec[i].coeffs[j] = ((byte >> bit) & 1) * 2 - 1; bit++;
     }
}

/*************************************************
* Name:        pack_polyvecm_sigma1
*
* Description: Packs a polynomial vector of length M with coefficients in {-1, 1}.
* Each coefficient is represented by a single bit.
*
* Arguments:   - uint8_t *buf: pointer to the output byte array
* - const polyvecm *s_m: pointer to the input polynomial vector
**************************************************/
void pack_polyvecm_sigma1(uint8_t *buf, const polyvecm *s_m) {
    size_t idx = 0; uint8_t byte = 0; int bit = 0;
    const size_t total_bytes = M * POLYSIGMA1_PACKEDBYTES;
    for (int i = 0; i < M; ++i) for (int j = 0; j < N; ++j) {
        byte |= ((s_m->vec[i].coeffs[j] + 1) >> 1) << bit;
        if (++bit == 8) { if(idx>=total_bytes) return; buf[idx++] = byte; byte = 0; bit = 0; }
    }
    if (bit > 0 && idx < total_bytes) buf[idx++] = byte;
}

/*************************************************
* Name:        unpack_polyvecm_sigma1
*
* Description: Unpacks a polynomial vector of length M with coefficients in {-1, 1}.
*
* Arguments:   - polyvecm *s_m: pointer to the output polynomial vector
* - const uint8_t *buf: pointer to the input byte array
**************************************************/
void unpack_polyvecm_sigma1(polyvecm *s_m, const uint8_t *buf) {
     size_t idx = 0; uint8_t byte = buf[0]; int bit = 0;
     const size_t total_bytes = M * POLYSIGMA1_PACKEDBYTES;
     for (int i = 0; i < M; ++i) for (int j = 0; j < N; ++j) {
          if (bit == 8) { if(++idx >= total_bytes && (i<M-1 || j<N-1)) return; byte = idx < total_bytes ? buf[idx] : 0; bit = 0; }
          s_m->vec[i].coeffs[j] = ((byte >> bit) & 1) * 2 - 1; bit++;
     }
}

/*************************************************
* Name:        pack_polyveck_sigma1
*
* Description: Packs a polynomial vector of length K with coefficients in {-1, 1}.
* Each coefficient is represented by a single bit.
*
* Arguments:   - uint8_t *buf: pointer to the output byte array
* - const polyveck *s_k: pointer to the input polynomial vector
**************************************************/
void pack_polyveck_sigma1(uint8_t *buf, const polyveck *s_k) {
    size_t idx = 0; uint8_t byte = 0; int bit = 0;
    const size_t total_bytes = K * POLYSIGMA1_PACKEDBYTES;
    for (int i = 0; i < K; ++i) for (int j = 0; j < N; ++j) {
        byte |= ((s_k->vec[i].coeffs[j] + 1) >> 1) << bit;
        if (++bit == 8) { if(idx>=total_bytes) return; buf[idx++] = byte; byte = 0; bit = 0; }
    }
    if (bit > 0 && idx < total_bytes) buf[idx++] = byte;
}

/*************************************************
* Name:        unpack_polyveck_sigma1
*
* Description: Unpacks a polynomial vector of length K with coefficients in {-1, 1}.
*
* Arguments:   - polyveck *s_k: pointer to the output polynomial vector
* - const uint8_t *buf: pointer to the input byte array
**************************************************/
void unpack_polyveck_sigma1(polyveck *s_k, const uint8_t *buf) {
     size_t idx = 0; uint8_t byte = buf[0]; int bit = 0;
     const size_t total_bytes = K * POLYSIGMA1_PACKEDBYTES;
     for (int i = 0; i < K; ++i) for (int j = 0; j < N; ++j) {
          if (bit == 8) { if(++idx >= total_bytes && (i<K-1 || j<N-1)) return; byte = idx < total_bytes ? buf[idx] : 0; bit = 0; }
          s_k->vec[i].coeffs[j] = ((byte >> bit) & 1) * 2 - 1; bit++;
     }
}

/*************************************************
* Name:        pack_poly_challenge
*
* Description: Packs a challenge polynomial c, which is sparse and has TAU
* coefficients equal to 1. Only the indices of the
* non-zero coefficients are stored.
*
* Arguments:   - uint8_t *buf: pointer to the output byte array (of length TAU)
* - const poly *c: pointer to the input challenge polynomial
**************************************************/
void pack_poly_challenge(uint8_t *buf, const poly *c) {
    uint32_t count = 0;
    for(int i=0; i<N && count < TAU; ++i) if(c->coeffs[i] == 1) buf[count++] = (uint8_t)i;
    if (count < TAU) memset(buf + count, 0, TAU - count);
}

/*************************************************
* Name:        unpack_poly_challenge
*
* Description: Unpacks a challenge polynomial c from a byte array containing
* the indices of its non-zero coefficients.
*
* Arguments:   - poly *c: pointer to the output challenge polynomial
* - const uint8_t *buf: pointer to the input byte array (of length TAU)
**************************************************/
void unpack_poly_challenge(poly *c, const uint8_t *buf) {
    poly_zero(c); uint8_t seen[N] = {0};
    for(int i=0; i<TAU; ++i) if (!seen[buf[i]]) { c->coeffs[buf[i]] = 1; seen[buf[i]] = 1; }
}

/*************************************************
* Name:        pack_polyveck_2q
*
* Description: Packs a polynomial vector of length K with coefficients in [0, 2q-1].
*
* Arguments:   - uint8_t *buf: pointer to the output byte array
* - const polyveck *v: pointer to the input polynomial vector
**************************************************/
void pack_polyveck_2q(uint8_t *buf, const polyveck *v) {
    const int bits = POLY2Q_BITS; size_t offset = 0;
    const size_t bytes_poly = POLY_PACKEDBYTES_2Q;
    for(int i=0; i<K; ++i){
        uint64_t B=0; int b=0; size_t idx=0;
        for(int j=0; j<N; ++j){
            uint32_t C = (uint32_t)v->vec[i].coeffs[j]; if(C>=DQ) C=DQ-1;
            B |= (uint64_t)C << b; b += bits;
            while(b>=8){ if(idx>=bytes_poly) return; buf[offset+idx++] = B&0xFF; B>>=8; b-=8; }
        }
        if(b>0 && idx<bytes_poly) buf[offset+idx++] = B&0xFF;
        while(idx<bytes_poly) buf[offset+idx++] = 0;
        offset += bytes_poly;
    }
}

/*************************************************
* Name:        pack_polyvecl_s_prime
*
* Description: Packs the polynomial vector s', where all coefficients
* are from {-1, 1}. Each coefficient is packed into a single bit.
*
* Arguments:   - uint8_t *buf: pointer to the output byte array
* - const polyvecl *s_prime: pointer to the input vector s'
**************************************************/
void pack_polyvecl_s_prime(uint8_t *buf, const polyvecl *s_prime) {
    size_t idx = 0;
    uint8_t byte = 0;
    int bit = 0;
    const size_t total_bytes = S_PRIME_PACKEDBYTES;

    for (int i = 0; i < L + K; ++i) {
        for (int j = 0; j < N; ++j) {
            // Map {-1, 1} to {0, 1}
            byte |= ((s_prime->vec[i].coeffs[j] + 1) >> 1) << bit;
            if (++bit == 8) {
                if (idx >= total_bytes) return; // Prevent buffer overflow
                buf[idx++] = byte;
                byte = 0;
                bit = 0;
            }
        }
    }
    if (bit > 0 && idx < total_bytes) {
        buf[idx] = byte;
    }
}

/*************************************************
* Name:        unpack_polyvecl_s_prime
*
* Description: Unpacks the polynomial vector s' into coefficients in {-1, 1}.
*
* Arguments:   - polyvecl *s_prime: pointer to the output vector s'
* - const uint8_t *buf: pointer to the input byte array
**************************************************/
void unpack_polyvecl_s_prime(polyvecl *s_prime, const uint8_t *buf) {
    size_t idx = 0;
    uint8_t byte = buf[0];
    int bit = 0;
    const size_t total_bytes = S_PRIME_PACKEDBYTES;

    for (int i = 0; i < L + K; ++i) {
        for (int j = 0; j < N; ++j) {
            if (bit == 8) {
                if (++idx >= total_bytes && (i < L + K -1 || j < N - 1)) return; // Bounds check
                byte = idx < total_bytes ? buf[idx] : 0;
                bit = 0;
            }
            // Map {0, 1} back to {-1, 1}
            s_prime->vec[i].coeffs[j] = ((byte >> bit) & 1) * 2 - 1;
            bit++;
        }
    }
}