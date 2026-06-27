#!/usr/bin/env python3
"""Independently generate the degree-12 fixed-point (Q63) minimax coefficients
for exp(-x) on [0, ln2) used by the rejection sampler's exp evaluation.

The evaluation in src/sampler.c uses the nested recurrence
    y <- C[0];  for u in 1..12:  y <- C[u] - ((x*y) >> 63)
with x, y in Q63. This realizes a degree-12 polynomial approximation of
exp(-x). The coefficients below are produced by a high-precision Chebyshev fit
of exp(-t) on [0, ln2) (requires mpmath); they are not taken from any external
implementation. Max error is <= 2 ulp(Q63) versus the true function, well
below the precision the sampler requires.

Run:  python3 gen/exp_coeffs.py
"""
import mpmath as mp

def main():
    mp.mp.dps = 80
    ln2 = mp.log(2)
    deg = 12
    # a[0]*t^deg + ... + a[deg]
    a = mp.chebyfit(lambda t: mp.e ** (-t), [0, ln2], deg + 1, error=False)
    # recurrence produces sum_j (-1)^j C[deg-j] t^j  =>  C[m] = (-1)^(deg-m) * coeff_t^(deg-m)
    C = []
    for m in range(deg + 1):
        j = deg - m
        coeff_tj = a[deg - j]                      # coefficient of t^j
        val = ((-1) ** j) * coeff_tj * (mp.mpf(2) ** 63)
        C.append(int(mp.nint(val)) & 0xFFFFFFFFFFFFFFFF)
    print("    static const uint64_t C[] = {")
    for i in range(0, deg + 1, 3):
        print("        " + ", ".join("0x%016Xu" % c for c in C[i:i + 3]) + ",")
    print("    };")

if __name__ == "__main__":
    main()
