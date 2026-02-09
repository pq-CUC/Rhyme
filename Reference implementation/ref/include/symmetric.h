#ifndef SYMMETRIC_H
#define SYMMETRIC_H

#include "params.h"
#include <stdint.h>
#include "fips202.h"
#include "aes256ctr.h" 

// Cryptographic XOF function: shake256
typedef keccak_state xof256_state;

#define ccc_shake256_absorb_twice                                           \
    ccc_NAMESPACE(ccc_shake256_absorb_twice)
void ccc_shake256_absorb_twice(keccak_state *state, const uint8_t *in1,
                                size_t in1len, const uint8_t *in2, size_t in2len);

#define XOF256_BLOCKBYTES SHAKE256_RATE

#define xof256_absorbe_once(STATE, IN, IN_LEN)                                  \
    shake256_absorb_once(STATE, IN, IN_LEN)
#define xof256_absorbe_twice(STATE, IN, IN_LEN, IN2, IN2_LEN)                   \
    ccc_shake256_absorb_twice(STATE, IN, IN_LEN, IN2, IN2_LEN) 
#define xof256_squeeze(OUT, OUT_LEN, STATE)                                     \
    shake256_squeeze(OUT, OUT_LEN, STATE)
#define xof256_squeezeblocks(OUT, OUTBLOCKS, STATE)                             \
    shake256_squeezeblocks(OUT, OUTBLOCKS, STATE)


// Stream function: aes256 or shake128|256
#ifdef ccc_USE_AES 

typedef aes256ctr_ctx stream128_state;
typedef aes256ctr_ctx stream256_state;

#define STREAM128_BLOCKBYTES AES256CTR_BLOCKBYTES
#define STREAM256_BLOCKBYTES AES256CTR_BLOCKBYTES

#define stream128_init(STATE, SEED, NONCE) \
        aes256ctr_init(STATE, SEED, NONCE)
#define stream128_squeezeblocks(OUT, OUTBLOCKS, STATE) \
        aes256ctr_squeezeblocks(OUT, OUTBLOCKS, STATE)
#define stream256_init(STATE, SEED, NONCE) \
        aes256ctr_init(STATE, SEED, NONCE)
#define stream256_squeezeblocks(OUT, OUTBLOCKS, STATE) \
        aes256ctr_squeezeblocks(OUT, OUTBLOCKS, STATE)

#else 

/* ---------------------------------------------------------
   Hybrid Mode: AES for Matrix Generation, SHAKE for others
   --------------------------------------------------------- */

/* 1. stream128 (For Agen) -> Forced to AES */
typedef aes256ctr_ctx stream128_state;
#define STREAM128_BLOCKBYTES AES256CTR_BLOCKBYTES

#define stream128_init(STATE, SEED, NONCE) \
        aes256ctr_init(STATE, SEED, NONCE)

#define stream128_free(STATE) \
        aes256ctr_free(STATE)

#define stream128_squeezeblocks(OUT, OUTBLOCKS, STATE) \
        aes256ctr_squeezeblocks(OUT, OUTBLOCKS, STATE)

/* 2. stream256 (For others) -> SHAKE256 */
typedef keccak_state stream256_state;
#define STREAM256_BLOCKBYTES SHAKE256_RATE

#define ccc_shake256_stream_init                                            \
    ccc_NAMESPACE(ccc_shake256_stream_init)
void ccc_shake256_stream_init(keccak_state *state,
                                 const uint8_t seed[CRHBYTES], uint16_t nonce);

#define stream256_init(STATE, SEED, NONCE)                                     \
    ccc_shake256_stream_init(STATE, SEED, NONCE)
#define stream256_squeezeblocks(OUT, OUTBLOCKS, STATE)                         \
    shake256_squeezeblocks(OUT, OUTBLOCKS, STATE)


#endif // stream

#endif //SYMMETRIC_H