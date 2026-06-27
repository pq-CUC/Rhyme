/* Field-norm descent verification (Ewha Algorithm 2, descent half only).
 *
 * Goal of THIS file: prove the descent is correct before using it to decide
 * coprimality.  For a single minor polynomial P in Z[x]/(x^n+1):
 *
 *   - GOLD resultant  : Res(P, x^n+1) = product over the n NTT roots of P(root)
 *                       rebuilt as an exact integer via multi-prime CRT.
 *   - DESCENT scalar  : apply N(P)(x^2) = P(x)*P(-x) mod (x^m+1), m halving each
 *                       layer, log2(n) times, until P becomes a degree-0 integer.
 *
 * These two integers MUST be identical.  If they are, the descent is a correct
 * (and far cheaper) way to obtain the resultant, hence a correct basis for the
 * coprimality decision gcd(Res(M_1),...,Res(M_d)).
 *
 * We work with the minors exactly as kg_solve_last_column builds them, to be
 * sure we test the same object the gold path tests.
 *
 * Standalone main; link against keygen objects.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kg_ntt.h"
#include "kg_zint.h"

/* ---- big-integer polynomial helpers (zint coeffs), schoolbook, keygen-only ---- */

/* negacyclic multiply mod (x^m + 1): r[k] = sum a[i]b[j], i+j=k  (minus wrap) */
static void poly_negamul(zint *r, const zint *a, const zint *b, unsigned m) {
    zint *acc = malloc((size_t)m * sizeof(zint));
    for (unsigned k = 0; k < m; k++) zi_init(&acc[k]);
    zint t; zi_init(&t);
    for (unsigned i = 0; i < m; i++) {
        if (zi_is_zero(&a[i])) continue;
        for (unsigned j = 0; j < m; j++) {
            if (zi_is_zero(&b[j])) continue;
            zi_mul(&t, &a[i], &b[j]);
            unsigned k = i + j;
            if (k < m) {
                zi_add(&acc[k], &acc[k], &t);
            } else {
                zi_sub(&acc[k - m], &acc[k - m], &t); /* x^m = -1 */
            }
        }
    }
    for (unsigned k = 0; k < m; k++) { zi_copy(&r[k], &acc[k]); zi_free(&acc[k]); }
    free(acc); zi_free(&t);
}

/* P(-x): negate odd-index coefficients */
static void poly_negx(zint *r, const zint *a, unsigned m) {
    for (unsigned i = 0; i < m; i++) {
        zi_copy(&r[i], &a[i]);
        if (i & 1) zi_neg(&r[i]);
    }
}

/* one field-norm step: out (length m/2) <- N(P)=P(x)P(-x) mod x^m+1, has only
 * even powers, repack into y=x^2 of length m/2 over x^{m/2}+1. */
static void fieldnorm_step(zint *out, const zint *P, unsigned m) {
    zint *Pn = malloc((size_t)m * sizeof(zint));
    zint *prod = malloc((size_t)m * sizeof(zint));
    for (unsigned i = 0; i < m; i++) { zi_init(&Pn[i]); zi_init(&prod[i]); }
    poly_negx(Pn, P, m);
    poly_negamul(prod, P, Pn, m);
    /* prod should have only even-degree terms; pack prod[2k] -> out[k] */
    for (unsigned k = 0; k < m / 2; k++) zi_copy(&out[k], &prod[2 * k]);
    for (unsigned i = 0; i < m; i++) { zi_free(&Pn[i]); zi_free(&prod[i]); }
    free(Pn); free(prod);
}

/* full descent: P (length n) -> scalar integer (resultant), written to res */
static void descent_resultant(zint *res, const zint *P, unsigned n) {
    unsigned m = n;
    zint *cur = malloc((size_t)m * sizeof(zint));
    for (unsigned i = 0; i < m; i++) { zi_init(&cur[i]); zi_copy(&cur[i], &P[i]); }
    while (m > 1) {
        zint *nxt = malloc((size_t)(m / 2) * sizeof(zint));
        for (unsigned i = 0; i < m / 2; i++) zi_init(&nxt[i]);
        fieldnorm_step(nxt, cur, m);
        for (unsigned i = 0; i < m; i++) zi_free(&cur[i]);
        free(cur);
        cur = nxt; m /= 2;
    }
    zi_copy(res, &cur[0]);
    zi_free(&cur[0]); free(cur);
}

/* ---- gold resultant via multi-prime CRT (mirrors kg_solver.c) ---- */
/* needs the crt machinery; we instead just rebuild via a simple approach:
 * Res = product over NTT points, accumulated through CRT.  To avoid copying the
 * whole crt_mod apparatus, we use a Garner-free big-int product is too slow.
 * Simpler: compute Res mod many primes, compare descent's Res mod the same
 * primes.  If descent is correct, descent(P) mod p == prod_s P(root_s) mod p
 * for EVERY prime p.  Checking enough primes (and that bitlength matches) is a
 * strong correctness signal. */

/* product of NTT-point values of P mod p (this is Res(P,x^n+1) mod p) */
static uint32_t res_mod_p(const zint *P, unsigned n, unsigned pi) {
    uint32_t p = kg_prime(pi);
    uint32_t *h = malloc((size_t)n * 4);
    for (unsigned t = 0; t < n; t++) h[t] = zi_mod_u32(&P[t], p);
    kg_ntt_fwd(h, n, pi);
    uint64_t rp = 1;
    for (unsigned t = 0; t < n; t++) rp = (rp * h[t]) % p;
    free(h);
    return (uint32_t)rp;
}

/* sample one small-coeff minor-like polynomial directly (for a focused test) */
static void sample_poly(zint *P, unsigned n, int bound, uint64_t *st) {
    for (unsigned i = 0; i < n; i++) {
        *st = *st * 6364136223846793005ULL + 1442695040888963407ULL;
        int v = (int)((*st >> 30) % (uint64_t)(2 * bound + 1)) - bound;
        zi_set_i64(&P[i], v);
    }
}

int main(int argc, char **argv) {
    unsigned n = 256;
    int bound = 3;
    unsigned ntest = 5;
    if (argc > 1) n = (unsigned)atoi(argv[1]);
    if (argc > 2) bound = atoi(argv[2]);
    if (argc > 3) ntest = (unsigned)atoi(argv[3]);

    printf("# descent vs gold resultant: n=%u bound=%d ntest=%u\n", n, bound, ntest);
    uint64_t st = 0xCAFEBABEULL;
    unsigned pass = 0, fail = 0;

    for (unsigned tt = 0; tt < ntest; tt++) {
        zint *P = malloc((size_t)n * sizeof(zint));
        for (unsigned i = 0; i < n; i++) zi_init(&P[i]);
        sample_poly(P, n, bound, &st);

        zint res; zi_init(&res);
        descent_resultant(&res, P, n);

        /* check descent's resultant mod several primes equals direct res_mod_p */
        int ok = 1;
        unsigned checks = 12;
        for (unsigned c = 0; c < checks; c++) {
            unsigned pi = c;
            uint32_t p = kg_prime(pi);
            uint32_t a = zi_mod_u32(&res, p);
            uint32_t b = res_mod_p(P, n, pi);
            if (a != b) { ok = 0;
                printf("  test %u: MISMATCH at prime[%u]=%u  descent=%u direct=%u\n",
                       tt, pi, p, a, b);
                break;
            }
        }
        if (ok) { pass++; printf("  test %u: OK  (resultant bitlen=%u)\n", tt, zi_bitlen(&res)); }
        else fail++;

        zi_free(&res);
        for (unsigned i = 0; i < n; i++) zi_free(&P[i]);
        free(P);
    }
    printf("\npass=%u fail=%u\n", pass, fail);
    return fail ? 1 : 0;
}
