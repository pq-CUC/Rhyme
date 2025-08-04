// include/encoding.h
#ifndef CCC_ENCODING_H
#define CCC_ENCODING_H

#include "params.h"
#include "polyvec.h"
#include <stdint.h>
#include <stddef.h>
#include "config.h"

#define encode_z ccc_NAMESPACE(encode_z)
uint16_t encode_z(uint8_t *buf, const polyvecl *z);

#define decode_z ccc_NAMESPACE(decode_z)
int decode_z(polyvecl *z, const uint8_t *buf, uint16_t size_in);

#endif // CCC_ENCODING_H