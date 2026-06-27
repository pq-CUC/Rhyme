/* Coprimality-test benchmark / correctness harness.
 *
 * Compares three deciders on the SAME sampled short-column sets:
 *   (A) GOLD   : exact-resultant + multi-xgcd  (current kg_solve_check_only path)
 *   (D) METHOD-D: bounded multi-prime resultant test, no CRT rebuild
 *
 * For each random sample of d*(d-1) short columns (n coeffs each), we record
 * whether each decider says "coprime/solvable" and how long it takes.
 *
 * This file is standalone (its own main); link against the keygen objects.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "kg_ntt.h"
#include "kg_zint.h"
#include "kg_solver.h"

/* det of an m x m matrix over F_p, row-major, entries in [0,p). */
static uint32_t detmod(uint32_t *M, unsigned m, uint32_t p) {
    uint32_t det = 1;
    for (unsigned col = 0; col < m; col++) {
        /* find pivot */
        unsigned piv = col;
        while (piv < m && M[piv * m + col] == 0) piv++;
        if (piv == m) return 0;
        if (piv != col) {
            for (unsigned j = 0; j < m; j++) {
                uint32_t t = M[piv * m + j];
                M[piv * m + j] = M[col * m + j];
                M[col * m + j] = t;
            }
            det = (uint32_t)(p - det); /* row swap flips sign */
            if (det == p) det = 0;
        }
        uint32_t pivval = M[col * m + col];
        /* inverse of pivot mod p */
        uint32_t inv = kg_powmod(pivval, p - 2, p);
        det = (uint32_t)(((uint64_t)det * pivval) % p);
        for (unsigned r = col + 1; r < m; r++) {
            uint32_t factor = (uint32_t)(((uint64_t)M[r * m + col] * inv) % p);
            if (!factor) continue;
            for (unsigned j = col; j < m; j++) {
                uint32_t sub = (uint32_t)(((uint64_t)factor * M[col * m + j]) % p);
                uint32_t cur = M[r * m + j];
                cur = (cur >= sub) ? (cur - sub) : (cur + p - sub);
                M[r * m + j] = cur;
            }
        }
    }
    return det;
}

/* gcd of two uint32 */
static uint32_t gcd_u32(uint32_t a, uint32_t b) {
    while (b) { uint32_t t = a % b; a = b; b = t; }
    return a;
}

/* METHOD-D core: bounded multi-prime coprimality test.
 *
 * The d minors M_0..M_{d-1} (each a degree-<n poly over Z) are coprime in
 * Z[x]/(x^n+1) iff gcd over i of Res(M_i, x^n+1) is 1.  Instead of rebuilding
 * each integer resultant via CRT, we probe with K distinct NTT-friendly
 * primes p_t: for each prime we compute Res(M_i) mod p_t (= product over the n
 * NTT points of the minor value at that point), giving residues R_{i,t}.
 *
 * Decision (deterministic up to the prime budget K):
 *   - If for SOME prime p_t, gcd_i(R_{i,t} mod p_t) == 1  -> coprime (accept).
 *   - If across ALL K primes every prime divides the running gcd -> not coprime.
 *
 * We track, per sample, the smallest #primes needed to decide accept.
 *
 * Returns 1 = coprime, 0 = not coprime, and writes #primes used to *kused.
 */
static int method_D(const int32_t *cols, unsigned n, unsigned d,
                    unsigned Kbudget, unsigned *kused) {
    unsigned dm1 = d - 1;
    uint32_t *H  = malloc((size_t)d * dm1 * n * 4);
    uint32_t *Mt = malloc((size_t)dm1 * dm1 * 4);
    /* per-minor accumulated integer-gcd surrogate: we keep, across primes,
     * the gcd of the *lifted* residues is not available without CRT.  Instead
     * we use the standard multi-prime trick: the integer Res(M_i) is coprime
     * across i iff there EXISTS a modulus under which the residue-vector gcd
     * is 1.  A single prime p where gcd_i(Res(M_i) mod p) == 1 already proves
     * integer coprimality (since any common integer factor would divide that
     * gcd mod p too).  So accept on first such prime. */
    int decided = 0, result = 0;
    unsigned t;
    for (t = 0; t < Kbudget; t++) {
        unsigned pi = t;
        if (pi >= kg_nprimes()) break;
        uint32_t p = kg_prime(pi);
        /* NTT each short column mod p */
        for (unsigned i = 0; i < d; i++)
            for (unsigned j = 0; j < dm1; j++) {
                uint32_t *h = H + ((size_t)i * dm1 + j) * n;
                const int32_t *src = cols + ((size_t)j * d + i) * n;
                for (unsigned u = 0; u < n; u++) {
                    int64_t v = src[u] % (int64_t)p;
                    if (v < 0) v += p;
                    h[u] = (uint32_t)v;
                }
                kg_ntt_fwd(h, n, pi);
            }
        /* Res(M_i) mod p = product over NTT points s of det of (d-1)x(d-1) minor */
        uint32_t gcd_acc = 0;
        for (unsigned di = 0; di < d; di++) {
            uint32_t res_i = 1;
            for (unsigned s = 0; s < n; s++) {
                unsigned rr = 0;
                for (unsigned i = 0; i < d; i++) {
                    if (i == di) continue;
                    for (unsigned j = 0; j < dm1; j++)
                        Mt[rr * dm1 + j] = H[((size_t)i * dm1 + j) * n + s];
                    rr++;
                }
                uint32_t dval = detmod(Mt, dm1, p);
                res_i = (uint32_t)(((uint64_t)res_i * dval) % p);
            }
            gcd_acc = gcd_u32(gcd_acc, res_i);
            if (gcd_acc == 1) break; /* already coprime mod p */
        }
        if (gcd_acc == 1) { decided = 1; result = 1; t++; break; }
        /* gcd_acc != 1 (==0 means p | every Res_i, or a shared factor mod p).
         * Inconclusive for "coprime"; try next prime. If exhausted, we cannot
         * prove coprime -> conservatively report not-coprime. */
    }
    if (!decided) { result = 0; t = (t > Kbudget ? Kbudget : t); }
    *kused = t;
    free(H); free(Mt);
    return result;
}

/* sample d*(d-1) short columns: CBD-like small ints in [-eta,eta] */
static void sample_cols(int32_t *cols, unsigned n, unsigned d, unsigned eta,
                        uint64_t *st) {
    unsigned dm1 = d - 1;
    for (size_t k = 0; k < (size_t)d * dm1 * n; k++) {
        /* simple CBD: sum of eta bits minus sum of eta bits */
        int v = 0;
        for (unsigned b = 0; b < eta; b++) {
            *st = *st * 6364136223846793005ULL + 1442695040888963407ULL;
            v += (int)((*st >> 33) & 1);
            *st = *st * 6364136223846793005ULL + 1442695040888963407ULL;
            v -= (int)((*st >> 33) & 1);
        }
        cols[k] = v;
    }
}

int main(int argc, char **argv) {
    unsigned n = 256, d = 4, eta = 2, nsamp = 200, Kbudget = 8;
    if (argc > 1) n = (unsigned)atoi(argv[1]);
    if (argc > 2) d = (unsigned)atoi(argv[2]);
    if (argc > 3) nsamp = (unsigned)atoi(argv[3]);
    if (argc > 4) eta = (unsigned)atoi(argv[4]);

    printf("# n=%u d=%u eta=%u nsamp=%u Kbudget=%u nprimes=%u\n",
           n, d, eta, nsamp, Kbudget, kg_nprimes());

    int32_t *cols = malloc((size_t)d * (d - 1) * n * sizeof(int32_t));
    int32_t *Fout = malloc((size_t)d * n * sizeof(int32_t));
    uint64_t st = 0x12345678ABCDEFULL;

    unsigned agree = 0, disagree = 0, gold_cop = 0, d_cop = 0;
    unsigned long ksum = 0; unsigned kmax = 0;
    double t_gold = 0, t_d = 0;
    struct timespec a, b;

    kg_solve_check_only = 1;
    for (unsigned s = 0; s < nsamp; s++) {
        sample_cols(cols, n, d, eta, &st);

        clock_gettime(CLOCK_MONOTONIC, &a);
        int rg = kg_solve_last_column(Fout, cols, n, d);
        clock_gettime(CLOCK_MONOTONIC, &b);
        t_gold += (b.tv_sec - a.tv_sec) + (b.tv_nsec - a.tv_nsec) * 1e-9;
        int gold = (rg == 0) ? 1 : 0;

        unsigned kused = 0;
        clock_gettime(CLOCK_MONOTONIC, &a);
        int dd = (kg_solve_coprime_descent(cols, n, d) == 0) ? 1 : 0;
        clock_gettime(CLOCK_MONOTONIC, &b);
        (void)kused; (void)Kbudget; (void)method_D;
        t_d += (b.tv_sec - a.tv_sec) + (b.tv_nsec - a.tv_nsec) * 1e-9;

        gold_cop += gold; d_cop += dd;
        ksum += kused; if (kused > kmax) kmax = kused;
        if (gold == dd) agree++; else {
            disagree++;
            if (disagree <= 8)
                printf("  DISAGREE sample %u: gold=%d D=%d (k=%u)\n",
                       s, gold, dd, kused);
        }
    }

    printf("\n== results ==\n");
    printf("samples           : %u\n", nsamp);
    printf("gold coprime      : %u (%.1f%%)\n", gold_cop, 100.0*gold_cop/nsamp);
    printf("descent coprime   : %u (%.1f%%)\n", d_cop, 100.0*d_cop/nsamp);
    printf("agree / disagree  : %u / %u\n", agree, disagree);
    printf("time gold total   : %.4f s  (%.3f ms/sample)\n", t_gold, 1e3*t_gold/nsamp);
    printf("time descent total: %.4f s  (%.3f ms/sample)\n", t_d, 1e3*t_d/nsamp);
    printf("speedup (gold/desc): %.2fx\n", t_d > 0 ? t_gold/t_d : 0);

    free(cols); free(Fout);
    return 0;
}
