#include "aes256ctr.h"
#include <string.h>
#include <stdlib.h>
#include <openssl/evp.h> 

/*************************************************
* Name:        make_iv
*
* Description: Helper function to construct the 16-byte Initialization Vector (IV)
* from a 64-bit nonce. The nonce is copied into the first 8 bytes,
* and the remaining bytes are zeroed.
*
* Arguments:   - uint8_t iv[16]: output IV buffer
* - uint64_t nonce: 64-bit nonce value
**************************************************/
static void make_iv(uint8_t iv[16], uint64_t nonce) {
    memset(iv, 0, 16);
    // Copy the 64-bit nonce to the first 8 bytes of the IV
    memcpy(iv, &nonce, sizeof(nonce)); 
}

/*************************************************
* Name:        aes256ctr_set_nonce
*
* Description: Resets the nonce (IV) in the AES state without re-initializing the key.
* This optimization avoids the overhead of key expansion when only
* the nonce changes.
*
* Arguments:   - aes256ctr_ctx *state: pointer to the AES-CTR state
* - uint64_t nonce: new 64-bit nonce
**************************************************/
void aes256ctr_set_nonce(aes256ctr_ctx *state, uint64_t nonce) {
    uint8_t iv[16];
    make_iv(iv, nonce);
    
    EVP_CIPHER_CTX *ctx = (EVP_CIPHER_CTX *)state->ctx;
    if (!ctx) return;

    // OpenSSL allows resetting the IV while keeping the Key by passing NULL for key.
    // This is more efficient than full initialization.
    EVP_EncryptInit_ex(ctx, NULL, NULL, NULL, iv);
}

/*************************************************
* Name:        aes256ctr_init
*
* Description: Initializes the AES-256-CTR state with a key and a nonce.
* Allocates the OpenSSL cipher context.
*
* Arguments:   - aes256ctr_ctx *state: pointer to the state to initialize
* - const uint8_t key[32]: pointer to the 32-byte secret key
* - uint64_t nonce: 64-bit nonce
**************************************************/
void aes256ctr_init(aes256ctr_ctx *state, const uint8_t key[32], uint64_t nonce) {
    uint8_t iv[16];
    make_iv(iv, nonce);

    // Allocate context
    EVP_CIPHER_CTX *ctx = EVP_CIPHER_CTX_new();
    state->ctx = (void *)ctx; 

    if (ctx) {
        // Initialize as AES-256-CTR
        EVP_EncryptInit_ex(ctx, EVP_aes_256_ctr(), NULL, key, iv);
    }
}

/*************************************************
* Name:        aes256ctr_squeezeblocks
*
* Description: Generates pseudo-random bytes by encrypting a zero-buffer.
* The output is written to the 'out' buffer.
*
* Arguments:   - uint8_t *out: pointer to the output buffer
* - size_t nblocks: number of 16-byte blocks to generate
* - aes256ctr_ctx *state: pointer to the AES-CTR state
**************************************************/
void aes256ctr_squeezeblocks(uint8_t *out, size_t nblocks, aes256ctr_ctx *state) {
    int outlen;
    size_t len = nblocks * AES256CTR_BLOCKBYTES;
    
    // Retrieve context
    EVP_CIPHER_CTX *ctx = (EVP_CIPHER_CTX *)state->ctx;
    if (!ctx) return;

    // Encrypt zero-bytes to generate the keystream (pseudo-random output)
    uint8_t *in_zeros = calloc(1, len); 
    if (!in_zeros) return;

    EVP_EncryptUpdate(ctx, out, &outlen, in_zeros, len);
    
    free(in_zeros);
}

/*************************************************
* Name:        aes256ctr_prf
*
* Description: AES-256-CTR based Pseudo-Random Function.
* Generates 'outlen' bytes of output given a seed and nonce.
* This is a self-contained function that initializes and frees the state.
*
* Arguments:   - uint8_t *out: pointer to output buffer
* - size_t outlen: number of bytes to generate
* - const uint8_t seed[32]: key/seed
* - uint64_t nonce: nonce value
**************************************************/
void aes256ctr_prf(uint8_t *out, size_t outlen, const uint8_t seed[32], uint64_t nonce) {
    aes256ctr_ctx state;
    aes256ctr_init(&state, seed, nonce);

    int len;
    uint8_t *in_zeros = calloc(1, outlen);
    
    EVP_CIPHER_CTX *ctx = (EVP_CIPHER_CTX *)state.ctx;
    if (ctx && in_zeros) {
        EVP_EncryptUpdate(ctx, out, &len, in_zeros, outlen);
    }
    
    if (in_zeros) free(in_zeros);
    aes256ctr_free(&state);
}

/*************************************************
* Name:        aes256ctr_free
*
* Description: Frees the OpenSSL cipher context associated with the state.
* Must be called to avoid memory leaks.
*
* Arguments:   - aes256ctr_ctx *state: pointer to the state to free
**************************************************/
void aes256ctr_free(aes256ctr_ctx *state) {
    if (state && state->ctx) {
        EVP_CIPHER_CTX_free((EVP_CIPHER_CTX *)state->ctx);
        state->ctx = NULL;
    }
}