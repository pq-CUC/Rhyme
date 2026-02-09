#ifndef CONFIG_H
#define CONFIG_H
// #define ccc_USE_AES
#if ccc_MODE == 2
#define CRYPTO_ALGNAME "ccc2"
#define ccc_NAMESPACETOP ccc2
#define ccc_NAMESPACE(s) cryptolab_ccc2_##s
#elif ccc_MODE == 3
#define CRYPTO_ALGNAME "ccc3"
#define ccc_NAMESPACETOP ccc3
#define ccc_NAMESPACE(s) cryptolab_ccc3_##s
#elif ccc_MODE == 5
#define CRYPTO_ALGNAME "ccc5"
#define ccc_NAMESPACETOP ccc5
#define ccc_NAMESPACE(s) cryptolab_ccc5_##s
#endif
#endif