#ifndef CCC_PARAMS_H // Changed guard name slightly
#define CCC_PARAMS_H

#include "config.h" // Keep for ccc_MODE and namespace

#define SEEDBYTES 32 // Random seed bytes
#define CRHBYTES 64  // Hash output bytes (for deriving internal seeds)
#define N 256        // Polynomial degree
#define CCC_RANS_Z 1 // rANS for z is ALWAYS enabled
#define CCC_RANS_SIZE_BYTES 2 // Size field for rANS data is always present

// --- CCC Core Parameters ---
#define Q 257      
#define DQ (Q << 1) 
#define mont_two 2
#if ccc_MODE == 2
    #define K 2         
    #define L 3         
    #define TAU 30      
    #define B_INFTY 32 
    
#elif ccc_MODE == 3 
    #define K 3         
    #define L 4         
    #define TAU 60      
    #define B_INFTY 48

#elif ccc_MODE == 5 
    #define K 4         
    #define L 5         
    #define TAU 128      
    #define B_INFTY 60
#else
    #error "Unsupported ccc_MODE"
#endif

#define M (L-1) 

#define POLYQ_BITS 9
#define POLYQ_PACKEDBYTES ((N * POLYQ_BITS + 7) / 8) 

#define POLYSIGMA1_BITS 1
#define POLYSIGMA1_PACKEDBYTES ((N * POLYSIGMA1_BITS + 7) / 8) 


#define S_PRIME_PACKEDBYTES ((((K + L) * N * 1) + 7) / 8)

#define S_BAR_PRIME_0_PACKEDBYTES POLYSIGMA1_PACKEDBYTES

#define POLYC_PACKEDBYTES TAU 

#define POLY2Q_BITS 10
#define POLY_PACKEDBYTES_2Q ((N * POLY2Q_BITS + 7) / 8) 

/* --- Crypto Byte Sizes (Now common for all modes) --- */
#define CRYPTO_PUBLICKEYBYTES (SEEDBYTES + K*N + (K*N + 7)/8 + 8)
#define CRYPTO_SECRETKEYBYTES (CRYPTO_PUBLICKEYBYTES + M * POLYSIGMA1_PACKEDBYTES + K * POLYSIGMA1_PACKEDBYTES + S_PRIME_PACKEDBYTES + S_BAR_PRIME_0_PACKEDBYTES + SEEDBYTES)

#define CRYPTO_BYTES 3800

#define MAX_SIGN_ITERATIONS 10 // Signing loop limit


#endif // CCC_PARAMS_H