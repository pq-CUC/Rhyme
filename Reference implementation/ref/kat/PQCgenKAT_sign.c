#include "api.h"
#include "randombytes.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h> 

#define MAX_MARKER_LEN 50

#define KAT_SUCCESS 0
#define KAT_FILE_OPEN_ERROR -1
#define KAT_DATA_ERROR -3
#define KAT_CRYPTO_FAILURE -4

int FindMarker(FILE *infile, const char *marker);
int ReadHex(FILE *infile, unsigned char *A, int Length, char *str);
void fprintBstr(FILE *fp, char *S, unsigned char *A, unsigned long long LEN); 

int main(int argc, const char **argv) { 
    char fn_req[64]; 
    char fn_rsp[64];
    FILE *fp_req = NULL; 
    FILE *fp_rsp = NULL; 
    unsigned char seed[48];
    unsigned char msg_buf[3300]; 
    unsigned char entropy_input[48];
    unsigned char *m = NULL, *sm = NULL, *m1 = NULL; 
    size_t mlen = 0, smlen = 0, mlen1 = 0, pklen = 0, sklen = 0; 
    int count = 0;
    int done = 0;
    unsigned char pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    int ret_val = 0;
    
    sprintf(fn_req, "PQCsignKAT_%d.req", CRYPTO_SECRETKEYBYTES);
    fp_req = fopen(fn_req, "w"); 
    if (fp_req == NULL) {
        printf("Error: Unable to open request file <%s> for writing\n", fn_req);
        return KAT_FILE_OPEN_ERROR;
    }

    for (int i = 0; i < 48; i++)
        entropy_input[i] = i;
    randombytes_init(entropy_input, NULL, 256);

    
    for (int i = 0; i < 100; i++) { 
        fprintf(fp_req, "count = %d\n", i);
        randombytes(seed, 48);
        fprintBstr(fp_req, "seed = ", seed, 48);
        mlen = 33 * (i + 1); 
        fprintf(fp_req, "mlen = %zu\n", mlen); 
        if (mlen > sizeof(msg_buf)) { 
            printf("Error: Message length %zu is too large for buffer of size %zu\n", mlen, sizeof(msg_buf));
            fclose(fp_req);
            return KAT_DATA_ERROR;
        }
        randombytes(msg_buf, mlen); 
        fprintBstr(fp_req, "msg = ", msg_buf, mlen);
        
        fprintf(fp_req, "pk =\n");
        fprintf(fp_req, "sk =\n");
        fprintf(fp_req, "smlen =\n");
        fprintf(fp_req, "sm =\n\n");
    }
    fclose(fp_req); 
    fp_req = NULL; 

    printf("Request file generated: %s\n", fn_req);
    sprintf(fn_rsp, "PQCsignKAT_%d.rsp", CRYPTO_SECRETKEYBYTES);

    
    fp_req = fopen(fn_req, "r");
    if (fp_req == NULL) {
        printf("Error: Unable to open request file <%s> for reading\n", fn_req);
        return KAT_FILE_OPEN_ERROR;
    }

    
    fp_rsp = fopen(fn_rsp, "w");
    if (fp_rsp == NULL) {
        printf("Error: Unable to open response file <%s> for writing\n", fn_rsp);
        fclose(fp_req); 
        return KAT_FILE_OPEN_ERROR;
    }

    fprintf(fp_rsp, "# %s\n\n", CRYPTO_ALGNAME); 

    done = 0;
    do {
        
        if (!FindMarker(fp_req, "count = ")) { 
             done = 1; 
             break;
        }
        fscanf(fp_req, "%d", &count);
        fprintf(fp_rsp, "count = %d\n", count);

        
        if (!ReadHex(fp_req, seed, 48, "seed = ")) {
            printf("Error: Could not read 'seed' from <%s> (count=%d)\n", fn_req, count);
            
            fclose(fp_req); fclose(fp_rsp); free(m); free(m1); free(sm);
            return KAT_DATA_ERROR;
        }
        fprintBstr(fp_rsp, "seed = ", seed, 48);

        
        randombytes_init(seed, NULL, 256); 

        
        if (!FindMarker(fp_req, "mlen = ")) {
             printf("Error: Could not read 'mlen' from <%s> (count=%d)\n", fn_req, count);
             fclose(fp_req); fclose(fp_rsp); free(m); free(m1); free(sm); return KAT_DATA_ERROR;
        }

        
        if (fscanf(fp_req, "%zu", &mlen) != 1) { 
             printf("Error: Format error while reading 'mlen' (count=%d)\n", count);
             fclose(fp_req); fclose(fp_rsp); free(m); free(m1); free(sm); return KAT_DATA_ERROR;
        }
        fprintf(fp_rsp, "mlen = %zu\n", mlen);
 
        free(m); free(m1); free(sm);
        m = NULL; m1 = NULL; sm = NULL; 

        m = calloc(mlen ? mlen : 1, sizeof(unsigned char)); 

        size_t sm_alloc_size = (mlen > SIZE_MAX - CRYPTO_BYTES) ? SIZE_MAX : mlen + CRYPTO_BYTES; 
        size_t m1_alloc_size = mlen ? mlen : 1; 

        sm = calloc(sm_alloc_size, sizeof(unsigned char));
        m1 = calloc(m1_alloc_size, sizeof(unsigned char));

        if (!m || !sm || !m1) {
            printf("Error: Memory allocation failed (count=%d)\n", count);
            fclose(fp_req); fclose(fp_rsp); free(m); free(m1); free(sm); return KAT_CRYPTO_FAILURE;
        }

        if (mlen > INT_MAX) {
             printf("Error: Message length %zu is too large for ReadHex (count=%d)\n", mlen, count);
             fclose(fp_req); fclose(fp_rsp); free(m); free(m1); free(sm); return KAT_DATA_ERROR;
        }
        if (!ReadHex(fp_req, m, (int)mlen, "msg = ")) {
            printf("Error: Could not read 'msg' from <%s> (count=%d)\n", fn_req, count);
            fclose(fp_req); fclose(fp_rsp); free(m); free(m1); free(sm); return KAT_DATA_ERROR;
        }
        fprintBstr(fp_rsp, "msg = ", m, mlen);

        
        if (!FindMarker(fp_req, "pk =")) { break; }
        if (!FindMarker(fp_req, "sk =")) {  break; }
        if (!FindMarker(fp_req, "smlen =")) {  break; }
        if (!FindMarker(fp_req, "sm =")) {  break; }

        
        if ((ret_val = kat_crypto_sign_keypair(pk, &pklen, sk, &sklen)) != 0) {
            printf("kat_crypto_sign_keypair returned <%d> (count=%d)\n", ret_val, count);
            fclose(fp_req); fclose(fp_rsp); free(m); free(m1); free(sm); return KAT_CRYPTO_FAILURE;
        }
        fprintBstr(fp_rsp, "pk = ", pk, pklen); 
        fprintBstr(fp_rsp, "sk = ", sk, sklen); 

        
        if ((ret_val = kat_crypto_sign(sm, &smlen, m, mlen, sk)) != 0) {
            printf("kat_crypto_sign_open returned <%d> (count=%d)\n", ret_val, count);
             fclose(fp_req); fclose(fp_rsp); free(m); free(m1); free(sm); return KAT_CRYPTO_FAILURE;
        }
        fprintf(fp_rsp, "smlen = %zu\n", smlen); 
        fprintBstr(fp_rsp, "sm = ", sm, smlen);
        fprintf(fp_rsp, "\n");

        
        
        if ((ret_val = kat_crypto_sign_open(m1, &mlen1, sm, smlen, pk, pklen)) != 0) {
            printf("kat_crypto_sign_open returned <%d> (count=%d)\n", ret_val, count);
             fclose(fp_req); fclose(fp_rsp); free(m); free(m1); free(sm); return KAT_CRYPTO_FAILURE;
        }

        if (mlen != mlen1) {
            printf("kat_crypto_sign_open returned bad 'mlen': got <%zu>, expected <%zu> (count=%d)\n", mlen1, mlen, count);
             fclose(fp_req); fclose(fp_rsp); free(m); free(m1); free(sm); return KAT_CRYPTO_FAILURE;
        }

        if (memcmp(m, m1, mlen)) {
            printf("kat_crypto_sign_open returned bad 'm' value (count=%d)\n", count);
             fclose(fp_req); fclose(fp_rsp); free(m); free(m1); free(sm); return KAT_CRYPTO_FAILURE;
        }
        
        free(m); m = NULL;
        free(m1); m1 = NULL;
        free(sm); sm = NULL;

    } while (!done);

    
    if (fp_req) fclose(fp_req);
    if (fp_rsp) fclose(fp_rsp);

    
    free(m);
    free(m1);
    free(sm);

    printf("Response file generated: %s\n", fn_rsp);
    return KAT_SUCCESS;
}

//
// ALLOW TO READ HEXADECIMAL ENTRY (KEYS, DATA, TEXT, etc.)
//
int FindMarker(FILE *infile, const char *marker) {
    char line[MAX_MARKER_LEN];
    int i, len;
    int curr_char; 

    if (!infile || !marker) return 0; 

    len = (int)strlen(marker);
    if (len == 0 || len > MAX_MARKER_LEN - 1) return 0; 

    
    for (i = 0; i < len; i++) {
        curr_char = fgetc(infile);
        if (curr_char == EOF) return 0; 
        line[i] = (char)curr_char;
    }
    line[len] = '\0';

    
    while (1) {
        if (strncmp(line, marker, len) == 0) return 1; 
        
        for (i = 0; i < len - 1; i++)
            line[i] = line[i + 1];
        
        curr_char = fgetc(infile);
        if (curr_char == EOF) return 0; 
        line[len - 1] = (char)curr_char;
        
    }

    
    return 0;
}

//
// ALLOW TO READ HEXADECIMAL ENTRY (KEYS, DATA, TEXT, etc.)
//
int ReadHex(FILE *infile, unsigned char *A, int Length, char *str) {
    int i, ch, started;
    unsigned char ich; 

    if (Length == 0) {
        
        return 1; 
    }
    if (!A) return 0; 
    memset(A, 0x00, Length); 
    started = 0; 

    if (FindMarker(infile, str)) {
        int nybble_count = 0; 
        unsigned char current_byte = 0; 

        while ((ch = fgetc(infile)) != EOF) {
            if (!isxdigit(ch)) { 
                if (!started) { 
                    if (ch == '\n') break; 
                    else continue; 
                } else break; 
            }
            started = 1; 

            if ((ch >= '0') && (ch <= '9')) ich = ch - '0';
            else if ((ch >= 'A') && (ch <= 'F')) ich = ch - 'A' + 10;
            else if ((ch >= 'a') && (ch <= 'f')) ich = ch - 'a' + 10;
            else ich = 0; 

            current_byte = (current_byte << 4) | ich;
            nybble_count++;

            if (nybble_count % 2 == 0) {
                int byte_index = (nybble_count / 2) - 1;
                if (byte_index < Length) {
                    A[byte_index] = current_byte;
                } else {
                     
                     return 0; 
                }
                current_byte = 0; 
            }
        }
        
        if (nybble_count % 2 != 0) {    
             
        }
         
         if (nybble_count / 2 != Length && started) { 
             
         }

    } else {
        return 0; 
    }

    return 1; 
}

void fprintBstr(FILE *fp, char *S, unsigned char *A, unsigned long long LEN) {
    unsigned long long i;

    fprintf(fp, "%s", S);

    for (i = 0; i < LEN; i++)
        fprintf(fp, "%02X", A[i]);

    if (LEN == 0) 
        fprintf(fp, "00");

    fprintf(fp, "\n");
}