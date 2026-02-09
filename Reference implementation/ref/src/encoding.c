#include "encoding.h"
#include "params.h"        // For K, L, N, B_INFTY
#include "polyvec.h"       // For polyvecl
#include "rans_byte.h"     // rANS definitions
#include "config.h"        // For ccc_NAMESPACE
#include <string.h>        // For memcpy
#include <stdint.h>
#include <stdio.h>         // For fprintf
#include <stdlib.h>        // For malloc/free
#include <limits.h>        // For UINT16_MAX
#include <inttypes.h>      // For PRId32

// --- rANS Parameters ---
#define SCALE_BITS_Z 16  
#define SCALE_Z (1U << SCALE_BITS_Z) //  65536 



#if ccc_MODE == 2 && TAU == 30
    #include "rans_tables_tau30_split.h" 
    #define DSYMS_Z0      dsyms_z0_tau30
    #define ESYMS_Z0      esyms_z0_tau30
    #define SYMBOL_Z0     symbol_z0_tau30
    #define DSYMS_Z_REST  dsyms_z_rest_tau30
    #define ESYMS_Z_REST  esyms_z_rest_tau30
    #define SYMBOL_Z_REST symbol_z_rest_tau30

#elif ccc_MODE == 3 && TAU == 60 
    #include "rans_tables_tau60_split.h"
    #define DSYMS_Z0      dsyms_z0_tau60
    #define ESYMS_Z0      esyms_z0_tau60
    #define SYMBOL_Z0     symbol_z0_tau60
    #define DSYMS_Z_REST  dsyms_z_rest_tau60
    #define ESYMS_Z_REST  esyms_z_rest_tau60
    #define SYMBOL_Z_REST symbol_z_rest_tau60

#elif ccc_MODE == 5 && TAU == 128
    #include "rans_tables_tau128_split.h"
    #define DSYMS_Z0      dsyms_z0_tau128
    #define ESYMS_Z0      esyms_z0_tau128
    #define SYMBOL_Z0     symbol_z0_tau128
    #define DSYMS_Z_REST  dsyms_z_rest_tau128
    #define ESYMS_Z_REST  esyms_z_rest_tau128
    #define SYMBOL_Z_REST symbol_z_rest_tau128
#else
    #error "Unsupported ccc_MODE or TAU value for rANS tables in encoding.c."
#endif


#define MIN_VAL (-(B_INFTY - 1))
#define MAX_VAL (B_INFTY - 1)
#define NUM_SYMBOLS_Z (MAX_VAL - MIN_VAL + 1)
#define OFFSET_Z (-MIN_VAL)

/*************************************************
* Name:        encode_z
*
* Description: Compresses a polynomial vector z using range Asymmetric Numeral
* Systems (rANS) entropy coding.
*
* Arguments:   - uint8_t *buf: pointer to the output byte array for the compressed data
* - const polyvecl *z: pointer to the input polynomial vector z
*
* Returns:     The length of the compressed byte array, or 0 on error.
**************************************************/
uint16_t encode_z(uint8_t *buf, const polyvecl *z) {
    size_t total_coeffs = (size_t)(K + L) * N;
    size_t temp_buf_size = total_coeffs * 3 + 256;
    uint8_t *encoding_buffer = malloc(temp_buf_size);
    if (!encoding_buffer) { return 0; }
    uint8_t *ptr = encoding_buffer + temp_buf_size;
    RansState rans;
    uint16_t encoded_size = 0;

    RansEncInit(&rans);

    /* Iterate over polynomial vector (reverse order) */
    for (size_t vec_idx_rev = 0; vec_idx_rev < K + L; ++vec_idx_rev) {
        size_t vec_idx = K + L - 1 - vec_idx_rev; 
        
        /* Select appropriate encoding and decoding tables */
        const RansEncSymbol *current_esyms;
        const RansDecSymbol *current_dsyms; 

        if (vec_idx == 0) {
            current_esyms = ESYMS_Z0;
            current_dsyms = DSYMS_Z0;       /* Select z0 table */
        } else {
            current_esyms = ESYMS_Z_REST;
            current_dsyms = DSYMS_Z_REST;  /* Select z_rest table */
        }

        for (size_t coeff_idx_rev = 0; coeff_idx_rev < N; ++coeff_idx_rev) {
            size_t coeff_idx = N - 1 - coeff_idx_rev;
            int32_t val = z->vec[vec_idx].coeffs[coeff_idx];
            
            if (val < MIN_VAL || val > MAX_VAL) {
                 fprintf(stderr, "[MODE %d] Error: Coeff %" PRId32 " out of encodable range [%d, %d] in encode_z!\n",
                         ccc_MODE, val, MIN_VAL, MAX_VAL);
                 free(encoding_buffer); return 0;
            }
            
            uint16_t symbol = (uint16_t)(val + OFFSET_Z);

            /* Use current_dsyms for frequency check */
            if (symbol >= NUM_SYMBOLS_Z || current_dsyms[symbol].freq == 0) {
                 fprintf(stderr, "[MODE %d] Error: Attempting to encode symbol %u (val %d) with freq 0!\n",
                         ccc_MODE, symbol, val);
                 free(encoding_buffer); return 0;
            }

            /* Use current_esyms for encoding */
            RansEncPutSymbol(&rans, &ptr, &current_esyms[symbol]);

            if (ptr < encoding_buffer) {
                 fprintf(stderr, "[MODE %d] Error: rANS temp buffer underflow!\n", ccc_MODE);
                 free(encoding_buffer); return 0;
            }
        }
    }
    RansEncFlush(&rans, &ptr);
    size_t size_encoded_temp = (encoding_buffer + temp_buf_size) - ptr;
    if (size_encoded_temp > UINT16_MAX) {
         fprintf(stderr, "[MODE %d] Error: Encoded size %zu exceeds uint16_t limit!\n", ccc_MODE, size_encoded_temp);
         free(encoding_buffer); return 0;
    }
    encoded_size = (uint16_t)size_encoded_temp;
    memcpy(buf, ptr, encoded_size);
    free(encoding_buffer);
    return encoded_size;
}

/*************************************************
* Name:        decode_z
*
* Description: Decompresses a polynomial vector z from a byte array using
* range Asymmetric Numeral Systems (rANS) entropy decoding.
*
* Arguments:   - polyvecl *z: pointer to the output (decompressed) polynomial vector z
* - const uint8_t *buf: pointer to the input byte array of compressed data
* - uint16_t size_in: the length of the compressed byte array
*
* Returns:     0 on success, non-zero on failure.
**************************************************/
int decode_z(polyvecl *z, const uint8_t *buf, uint16_t size_in) {
    RansState rans;
    const uint8_t *ptr = buf;
    const uint8_t *buf_end = buf + size_in;
    size_t total_coeffs = (size_t)(K + L) * N;
    size_t coeffs_decoded = 0;
    
    if (size_in < 4) { return 1; }
    if (RansDecInit(&rans, (uint8_t **)&ptr)) { return 2; }

    for (size_t vec_idx = 0; vec_idx < K + L; ++vec_idx) {
        const RansDecSymbol *current_dsyms;
        const uint16_t *current_lut;
        
        if (vec_idx == 0) {
            current_dsyms = DSYMS_Z0;
            current_lut   = SYMBOL_Z0;
        } else {
            current_dsyms = DSYMS_Z_REST;
            current_lut   = SYMBOL_Z_REST;
        }
        for (size_t coeff_idx = 0; coeff_idx < N; ++coeff_idx) {
            if (ptr > buf_end) { return 9; }
            uint16_t slot = RansDecGet(&rans, SCALE_BITS_Z);
            uint16_t symbol = current_lut[slot];
            
            if (symbol >= NUM_SYMBOLS_Z) { return 4; }
            
            int32_t val = (int32_t)symbol - OFFSET_Z;
            z->vec[vec_idx].coeffs[coeff_idx] = val;
            coeffs_decoded++;
            
            RansDecAdvanceSymbol(&rans, (uint8_t **)&ptr, buf_end, &current_dsyms[symbol], SCALE_BITS_Z);
            
            if (ptr > buf_end && coeffs_decoded < total_coeffs) { return 5; }
        }
    }

    if (RansDecVerify(&rans)) { return 6; }
    if (ptr != buf_end) { return 7; }
    if (coeffs_decoded != total_coeffs) { return 8; }
    
    return 0;
}
