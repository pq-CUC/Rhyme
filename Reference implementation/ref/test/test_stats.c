#include <stdio.h>
#include <math.h>
#include <string.h>
#include "params.h"
#include "poly.h"
#include "sampler.h"
#include "symmetric.h"

int main(void) {
    uint8_t seed[CRHBYTES] = {9};
    poly y, z, c, v;
    double s1 = 0, s2 = 0; long cnt = 0;
    /* y1 CDT std check */
    for (int t = 0; t < 200; t++) {
        SampleY1(&y, seed, (uint16_t)t);
        for (int i = 0; i < N; i++) { double x = y.coeffs[i]; s1 += x; s2 += x*x; cnt++; }
    }
    printf("y1: mean %.3f std %.3f (expect ~0, ~%d)\n", s1/cnt, sqrt(s2/cnt), RHYME_SIGMAY);
    /* sigma_base CDT */
    s1 = s2 = 0; cnt = 0;
    for (int t = 0; t < 200; t++) {
        SampleGauss(&y, seed, (uint16_t)(1000 + t));
        for (int i = 0; i < N; i++) { double x = y.coeffs[i]; s1 += x; s2 += x*x; cnt++; }
    }
    printf("gauss: mean %.4f std %.4f (expect ~0, ~1.30)\n", s1/cnt, sqrt(s2/cnt));
    /* Algorithm 5: acceptance, z1 std, parity */
    int acc = 0, tot = 600;
    s1 = s2 = 0; cnt = 0;
    long parity_bad = 0; int maxz = 0;
    uint8_t cseed[CRHBYTES] = {7};
    for (int t = 0; t < tot; t++) {
        SampleY1(&y, seed, (uint16_t)(2000 + t));
        cseed[1] = (uint8_t)t; cseed[2] = (uint8_t)(t >> 8);
        SampleChallenge(&c, cseed);
        stream256_state st;
        uint8_t rseed[CRHBYTES] = {3, (uint8_t)t, (uint8_t)(t>>8)};
        stream256_init(&st, rseed, 0);
        if (Sample_Z(&z, &v, &c, &y, &st)) {
            acc++;
            for (int i = 0; i < N; i++) {
                double x = z.coeffs[i]; s1 += x; s2 += x*x; cnt++;
                if (abs(z.coeffs[i]) > maxz) maxz = abs(z.coeffs[i]);
                if (((v.coeffs[i] - c.coeffs[i]) & 1) != 0) parity_bad++;
                if (abs(v.coeffs[i]) > RHYME_R) parity_bad++;
            }
        }
    }
    printf("Alg5: acc %.3f (design ~>0.5incl-L2), z1 std %.2f, max|z1| %d (B0=%d), parity_bad %ld\n",
           (double)acc/tot, sqrt(s2/cnt), maxz, B0, parity_bad);
    return 0;
}
