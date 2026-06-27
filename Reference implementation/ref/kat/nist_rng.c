/* NIST AES256-CTR-DRBG for KAT generation, with a built-in table-based
 * AES-256 (encrypt direction only).  Bit-compatible with the NIST rng.c
 * shipped with PQC submissions, but without the OpenSSL dependency. */
#include <string.h>
#include "nist_rng.h"

/* ------------------------------------------------------------- AES-256 */
static const uint8_t SBOX[256] = {
0x63,0x7c,0x77,0x7b,0xf2,0x6b,0x6f,0xc5,0x30,0x01,0x67,0x2b,0xfe,0xd7,0xab,0x76,
0xca,0x82,0xc9,0x7d,0xfa,0x59,0x47,0xf0,0xad,0xd4,0xa2,0xaf,0x9c,0xa4,0x72,0xc0,
0xb7,0xfd,0x93,0x26,0x36,0x3f,0xf7,0xcc,0x34,0xa5,0xe5,0xf1,0x71,0xd8,0x31,0x15,
0x04,0xc7,0x23,0xc3,0x18,0x96,0x05,0x9a,0x07,0x12,0x80,0xe2,0xeb,0x27,0xb2,0x75,
0x09,0x83,0x2c,0x1a,0x1b,0x6e,0x5a,0xa0,0x52,0x3b,0xd6,0xb3,0x29,0xe3,0x2f,0x84,
0x53,0xd1,0x00,0xed,0x20,0xfc,0xb1,0x5b,0x6a,0xcb,0xbe,0x39,0x4a,0x4c,0x58,0xcf,
0xd0,0xef,0xaa,0xfb,0x43,0x4d,0x33,0x85,0x45,0xf9,0x02,0x7f,0x50,0x3c,0x9f,0xa8,
0x51,0xa3,0x40,0x8f,0x92,0x9d,0x38,0xf5,0xbc,0xb6,0xda,0x21,0x10,0xff,0xf3,0xd2,
0xcd,0x0c,0x13,0xec,0x5f,0x97,0x44,0x17,0xc4,0xa7,0x7e,0x3d,0x64,0x5d,0x19,0x73,
0x60,0x81,0x4f,0xdc,0x22,0x2a,0x90,0x88,0x46,0xee,0xb8,0x14,0xde,0x5e,0x0b,0xdb,
0xe0,0x32,0x3a,0x0a,0x49,0x06,0x24,0x5c,0xc2,0xd3,0xac,0x62,0x91,0x95,0xe4,0x79,
0xe7,0xc8,0x37,0x6d,0x8d,0xd5,0x4e,0xa9,0x6c,0x56,0xf4,0xea,0x65,0x7a,0xae,0x08,
0xba,0x78,0x25,0x2e,0x1c,0xa6,0xb4,0xc6,0xe8,0xdd,0x74,0x1f,0x4b,0xbd,0x8b,0x8a,
0x70,0x3e,0xb5,0x66,0x48,0x03,0xf6,0x0e,0x61,0x35,0x57,0xb9,0x86,0xc1,0x1d,0x9e,
0xe1,0xf8,0x98,0x11,0x69,0xd9,0x8e,0x94,0x9b,0x1e,0x87,0xe9,0xce,0x55,0x28,0xdf,
0x8c,0xa1,0x89,0x0d,0xbf,0xe6,0x42,0x68,0x41,0x99,0x2d,0x0f,0xb0,0x54,0xbb,0x16
};
static const uint8_t RCON[15] = {0x00,0x01,0x02,0x04,0x08,0x10,0x20,0x40,
                                 0x80,0x1b,0x36,0x6c,0xd8,0xab,0x4d};

typedef struct { uint8_t rk[15][16]; } aes256_ctx;

static void aes256_key_expand(aes256_ctx *ctx, const uint8_t key[32]) {
    uint8_t w[60][4];
    memcpy(w, key, 32);
    for (int i = 8; i < 60; i++) {
        uint8_t t[4];
        memcpy(t, w[i - 1], 4);
        if (i % 8 == 0) {
            uint8_t tmp = t[0];
            t[0] = (uint8_t)(SBOX[t[1]] ^ RCON[i / 8]);
            t[1] = SBOX[t[2]];
            t[2] = SBOX[t[3]];
            t[3] = SBOX[tmp];
        } else if (i % 8 == 4) {
            for (int j = 0; j < 4; j++) t[j] = SBOX[t[j]];
        }
        for (int j = 0; j < 4; j++) w[i][j] = (uint8_t)(w[i - 8][j] ^ t[j]);
    }
    for (int r = 0; r < 15; r++) memcpy(ctx->rk[r], w[4 * r], 16);
}

static uint8_t xt(uint8_t x) { return (uint8_t)((x << 1) ^ ((x >> 7) * 0x1b)); }

static void aes256_encrypt_block(const aes256_ctx *ctx, const uint8_t in[16], uint8_t out[16]) {
    uint8_t s[16];
    for (int i = 0; i < 16; i++) s[i] = (uint8_t)(in[i] ^ ctx->rk[0][i]);
    for (int round = 1; round <= 14; round++) {
        uint8_t t[16];
        /* SubBytes + ShiftRows */
        for (int c = 0; c < 4; c++)
            for (int r = 0; r < 4; r++)
                t[4 * c + r] = SBOX[s[4 * ((c + r) & 3) + r]];
        if (round < 14) {
            /* MixColumns */
            for (int c = 0; c < 4; c++) {
                uint8_t a0 = t[4*c], a1 = t[4*c+1], a2 = t[4*c+2], a3 = t[4*c+3];
                s[4*c]   = (uint8_t)(xt(a0) ^ xt(a1) ^ a1 ^ a2 ^ a3);
                s[4*c+1] = (uint8_t)(a0 ^ xt(a1) ^ xt(a2) ^ a2 ^ a3);
                s[4*c+2] = (uint8_t)(a0 ^ a1 ^ xt(a2) ^ xt(a3) ^ a3);
                s[4*c+3] = (uint8_t)(xt(a0) ^ a0 ^ a1 ^ a2 ^ xt(a3));
            }
        } else {
            memcpy(s, t, 16);
        }
        for (int i = 0; i < 16; i++) s[i] = (uint8_t)(s[i] ^ ctx->rk[round][i]);
    }
    memcpy(out, s, 16);
}

/* expose for self-test */
void nist_aes256_ecb(const uint8_t key[32], const uint8_t in[16], uint8_t out[16]) {
    aes256_ctx ctx;
    aes256_key_expand(&ctx, key);
    aes256_encrypt_block(&ctx, in, out);
}

/* ------------------------------------------------------------- NIST DRBG */
typedef struct {
    uint8_t key[32];
    uint8_t v[16];
    int reseed_counter;
} drbg_ctx;
static drbg_ctx DRBG;

static void drbg_update(const uint8_t *provided, uint8_t key[32], uint8_t v[16]) {
    uint8_t temp[48];
    aes256_ctx ctx;
    aes256_key_expand(&ctx, key);
    for (int i = 0; i < 3; i++) {
        for (int j = 15; j >= 0; j--) {
            if (v[j] == 0xff) v[j] = 0x00;
            else { v[j]++; break; }
        }
        aes256_encrypt_block(&ctx, v, temp + 16 * i);
    }
    if (provided)
        for (int i = 0; i < 48; i++) temp[i] ^= provided[i];
    memcpy(key, temp, 32);
    memcpy(v, temp + 32, 16);
}

void randombytes_init(unsigned char *entropy_input,
                      unsigned char *personalization_string,
                      int security_strength) {
    uint8_t seed_material[48];
    (void)security_strength;
    memcpy(seed_material, entropy_input, 48);
    if (personalization_string)
        for (int i = 0; i < 48; i++) seed_material[i] ^= personalization_string[i];
    memset(DRBG.key, 0, 32);
    memset(DRBG.v, 0, 16);
    drbg_update(seed_material, DRBG.key, DRBG.v);
    DRBG.reseed_counter = 1;
}

int randombytes(uint8_t *x, size_t xlen) {
    uint8_t block[16];
    aes256_ctx ctx;
    aes256_key_expand(&ctx, DRBG.key);
    size_t i = 0;
    while (xlen > 0) {
        for (int j = 15; j >= 0; j--) {
            if (DRBG.v[j] == 0xff) DRBG.v[j] = 0x00;
            else { DRBG.v[j]++; break; }
        }
        aes256_encrypt_block(&ctx, DRBG.v, block);
        size_t take = xlen > 16 ? 16 : xlen;
        memcpy(x + i, block, take);
        i += take;
        xlen -= take;
    }
    drbg_update(NULL, DRBG.key, DRBG.v);
    DRBG.reseed_counter++;
    return 0;
}
