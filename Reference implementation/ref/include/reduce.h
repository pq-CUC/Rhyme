
#ifndef CCC_REDUCE_H
#define CCC_REDUCE_H

#include "params.h"
#include <stdint.h>

#define MONT 1 
#define QINV -255

#include "config.h" 
#define montgomery_reduce ccc_NAMESPACE(montgomery_reduce)
int32_t montgomery_reduce(int64_t a);

#define barrett_reduce ccc_NAMESPACE(barrett_reduce)
int16_t barrett_reduce(int16_t a);

#define caddq ccc_NAMESPACE(caddq)
int32_t caddq(int32_t a);

#define freeze ccc_NAMESPACE(freeze)
int32_t freeze(int32_t a);

int32_t cryptolab_ccc_reduce32_2q_centered(int32_t a);
int32_t cryptolab_ccc_freeze2q_nonnegative(int32_t a);

#define reduce32_2q cryptolab_ccc_reduce32_2q_centered
#define freeze2q cryptolab_ccc_freeze2q_nonnegative

#endif // CCC_REDUCE_H