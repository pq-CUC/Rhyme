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
* Name:        pack_pk
*
* Description: Packs the public key into a byte array without a length field.
* The format is | seedA | fixed-size main data | overflow data |.
* This version is modified to be more space-efficient, like the Lore scheme.
*
* Arguments:   - uint8_t *pk: pointer to output byte array
* - const uint8_t seedA[]: pointer to seed for matrix A
* - const polyveck *b: pointer to polynomial vector b
*
* Returns:     The total number of bytes written to pk.
**************************************************/
size_t pack_pk(uint8_t *pk, const uint8_t seedA[SEEDBYTES], const polyveck *b) {
    uint8_t *pk_ptr = pk;
    
    // 1. Copy seedA
    memcpy(pk_ptr, seedA, SEEDBYTES);
    pk_ptr += SEEDBYTES;

    // 2. Prepare pointers for fixed-size main data and temporary overflow buffer
    uint8_t *main_buf_ptr = pk_ptr;
    uint8_t overflow_buf[(K * N + 7) / 8] = {0};
    int total_overflow_bits = 0;

    // 3. Process all polynomials into a fixed-size main area and a variable overflow area
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

    // 4. Advance pointer past the fixed-size main data area
    pk_ptr += K * N;

    // 5. Append the overflow data
    size_t overflow_byte_len = (total_overflow_bits + 7) / 8;
    memcpy(pk_ptr, overflow_buf, overflow_byte_len);
    pk_ptr += overflow_byte_len;

    // 6. Return the total actual size
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
    if (pklen < SEEDBYTES + K * N) return -1; // Basic check for minimal possible size

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
* Description: Packs the secret key. This function remains largely unchanged as
* it just copies the already packed public key.
**************************************************/
size_t pack_sk(uint8_t *sk, const uint8_t *pk, size_t pklen, const polyvecm *s_gen, const polyveck *e_gen, const polyvecl *s_prime, const poly *s_bar_prime_0, const uint8_t key[SEEDBYTES]) {
    // This function's logic does not need to change.
    // It correctly uses the pklen returned by the new pack_pk.
    uint8_t *sk_ptr = sk;
    const size_t sgen_bytes = M * POLYSIGMA1_PACKEDBYTES;
    const size_t egen_bytes = K * POLYSIGMA1_PACKEDBYTES;
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
* Description: Unpacks a secret key. This function is modified to deduce the
* public key length instead of reading it from the buffer.
**************************************************/
int unpack_sk(uint8_t *pk, size_t *pklen, polyvecm *s_gen, polyveck *e_gen, polyvecl *s_prime, poly *s_bar_prime_0, uint8_t key[SEEDBYTES], const uint8_t *sk, size_t sklen) {
    const uint8_t *sk_ptr = sk;
    const size_t sgen_bytes = M * POLYSIGMA1_PACKEDBYTES;
    const size_t egen_bytes = K * POLYSIGMA1_PACKEDBYTES;
    const size_t sprime_bytes = S_PRIME_PACKEDBYTES; 
    const size_t s_bar_prime_bytes = S_BAR_PRIME_0_PACKEDBYTES;

    // --- Start of modification ---
    // Deduce the length of the embedded public key
    if (sklen < SEEDBYTES + K*N) return -1; // Not enough data for even a minimal PK
    
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
    // --- End of modification ---

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
size_t pack_sig(uint8_t sig[CRYPTO_BYTES], const polyvecl *z, const poly *c) {
    uint8_t *sig_ptr = sig;
    uint16_t size_encoded_z;

    // Encode z first into the buffer starting after the size field
    size_encoded_z = encode_z(sig_ptr + CCC_RANS_SIZE_BYTES, z);
    if (size_encoded_z == 0) {
        fprintf(stderr, "[MODE %d] Error: encode_z failed in pack_sig.\n", ccc_MODE);
        return 0; // Return 0 to indicate error
    }

    // Calculate total space needed
    size_t total_needed = CCC_RANS_SIZE_BYTES + size_encoded_z + POLYC_PACKEDBYTES;

    // Check if the required space exceeds the allocated buffer size
    if (total_needed > CRYPTO_BYTES) {
        fprintf(stderr, "[MODE %d] Error: Packed sig size needed (%zu) > CRYPTO_BYTES (%d). z_size=%u\n",
                ccc_MODE, total_needed, CRYPTO_BYTES, size_encoded_z);
        return 0; // Return 0 to indicate error
    }

    // Write the size of encoded z (2 bytes)
    sig_ptr[0] = (uint8_t)(size_encoded_z & 0xFF);
    sig_ptr[1] = (uint8_t)((size_encoded_z >> 8) & 0xFF);
    sig_ptr += CCC_RANS_SIZE_BYTES;

    // Advance pointer past the encoded z data (already written by encode_z)
    sig_ptr += size_encoded_z;

    // Pack the challenge c
    pack_poly_challenge(sig_ptr, c);
    sig_ptr += POLYC_PACKEDBYTES;
    // Calculate the final packed length
    size_t packed_len = (size_t)(sig_ptr - sig);

    return packed_len; // Return the actual packed length
}
// **** END MODIFIED pack_sig ****


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
int unpack_sig(polyvecl *z, poly *c, const uint8_t *sig, size_t siglen) {
    const uint8_t *sig_ptr = sig;
    uint16_t size_encoded_z;
    size_t min_required_len = CCC_RANS_SIZE_BYTES + POLYC_PACKEDBYTES;

    // Check if siglen is at least minimum possible length
    if (siglen < min_required_len) return -1; // Signature too short

    // Read the size of the encoded z component
    size_encoded_z = (uint16_t)sig_ptr[0] | ((uint16_t)sig_ptr[1] << 8);
    sig_ptr += CCC_RANS_SIZE_BYTES;

    // Check if the claimed total length makes sense
    size_t expected_total_len = CCC_RANS_SIZE_BYTES + size_encoded_z + POLYC_PACKEDBYTES;
    if (siglen < expected_total_len) return -2; // siglen is shorter than implied by encoded_z size
    // Optional: Check if expected_total_len exceeds max buffer size if needed elsewhere,
    // but here we only care about decoding what's given in siglen.
    // Check if there's enough data left for the claimed size_encoded_z
    if ((size_t)(sig_ptr - sig) + size_encoded_z > siglen) return -2; // Not enough data for z

    // Decode z
    int decode_status = decode_z(z, sig_ptr, size_encoded_z);
    if (decode_status != 0) return -3; // rANS decoding failed
    sig_ptr += size_encoded_z;

    // Check if there's enough data left for the challenge c
    if ((size_t)(sig_ptr - sig) + POLYC_PACKEDBYTES > siglen) return -2; // Not enough data for c

    // Unpack the challenge c
    unpack_poly_challenge(c, sig_ptr);
    sig_ptr += POLYC_PACKEDBYTES;

    // Check if we consumed exactly siglen bytes (no trailing non-zero bytes expected if sender used correct length)
    // Note: Trailing zero bytes might exist if sender padded, but shouldn't affect verification.
    // We've already checked siglen >= expected_total_len. If siglen > expected_total_len,
    // it implies trailing padding, which we can ignore or optionally check for zeros.
    if ((size_t)(sig_ptr - sig) != expected_total_len) {
       // This case might occur if siglen provided was larger than necessary.
       // We successfully decoded based on the internal length field.
       // Let's not treat this as an error, but maybe issue a warning if needed.
       // fprintf(stderr, "Warning: unpack_sig consumed %zu bytes, but siglen was %zu\n", (size_t)(sig_ptr - sig), siglen);
    }
     // Optional check for non-zero trailing bytes if strict conformance is required:
     for (size_t i = (size_t)(sig_ptr - sig); i < siglen; ++i) {
         if (sig[i] != 0) { return -4; /* Trailing non-zero bytes */ }
     }

    return 0; // Success
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