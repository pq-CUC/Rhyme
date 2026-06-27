#ifndef RHYME_ENCODING_H
#define RHYME_ENCODING_H
#include <stdint.h>
#include <stddef.h>
#include "params.h"
#include "poly.h"

#define encode_z RHYME_NAMESPACE(encode_z)
#define decode_z RHYME_NAMESPACE(decode_z)

/* Compress (z1, z_rest) into buf (capacity cap).
 * Layout: [u16 rans_len][rans stream][raw low-bit stream].
 * Returns total bytes written, or 0 on overflow/error. */
size_t encode_z(uint8_t *buf, size_t cap, const poly *z1, const poly z_rest[D_REST]);

/* Decompress; len must be exact. Returns 0 ok, nonzero error. */
int decode_z(poly *z1, poly z_rest[D_REST], const uint8_t *buf, size_t len);

#endif
