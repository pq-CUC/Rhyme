#!/usr/bin/env python3
"""Generate rANS frequency tables for Rhyme signatures (cut-F, split coding).

Under Construction 4 ("cut-F") every retained z_rest row is a SHORT-row output
with the SAME marginal statistics, so there is no separate long F-row.  Two
coded streams remain:
  z1     : |v| <= B0, Gaussian sigma_y.  v' = v + B0 in [0, 2*B0]; hi raw-split.
  z_rest : Gaussian, marginal std s_rest = sqrt(v_marg),
           v_marg = sigma_base^2 * (1 + 4*D*N*(eta/2))  [D = widened matrix col count]
           v' = v + CENTER_R; hi-split with an ESC symbol for the tail.

The "big" table (zb) is emitted identical to the small one (zs) so the existing
encoding format/machinery stays intact while all z_rest rows code uniformly.
"""
import math
from decimal import Decimal, getcontext
getcontext().prec = 60

SCALE_BITS = 16
SCALE = 1 << SCALE_BITS

# Table 7 parameters (mirrors gen/params_tables.py)
MODES = {
    128: dict(n=256,  q=3329,  k=2, l=3, eta=2,  B0=388, sigma_y=122.0, sigma_base=1.295870, DMAT=10),
    256: dict(n=512,  q=9473,  k=2, l=3, eta=4,  B0=598, sigma_y=180.0, sigma_base=1.316036, DMAT=10),
    384: dict(n=512,  q=11777, k=3, l=4, eta=5,  B0=605, sigma_y=182.0, sigma_base=1.326374, DMAT=14),
    512: dict(n=1024, q=18433, k=2, l=3, eta=6,  B0=916, sigma_y=266.0, sigma_base=1.335898, DMAT=10),
}

def row_sigma_rest(m):
    n, eta, sb = m["n"], m["eta"], m["sigma_base"]
    D = m["DMAT"]
    # e_bottom is sampled at width 2*sigma_base, so the identity-block variance
    # coefficient is 4 (not 1):  v_marg = sb^2 * (4 + 4*D*n*(eta/2)).
    v_marg = sb * sb * (4.0 + 4.0 * D * n * (eta / 2.0))
    return math.sqrt(v_marg)

def real_bits_per_coeff(nbuckets, base, L, sigma, reserve=0):
    """Expected code length in bits/coeff actually spent by the encoder:
    L verbatim low bits, plus the high part (hi = v' >> L) coded with the SAME
    quantized frequencies the emitted table uses (via bucket_probs + normalize).
    This is the real cost, including integer-frequency rounding, not the ideal
    Shannon entropy."""
    ps = bucket_probs(nbuckets, base, L, sigma)
    f = normalize(ps, reserve=reserve)
    bits = 0.0
    for p, fi in zip(ps, f):
        pf = float(p)
        if pf > 0.0:
            bits += pf * (L - math.log2(fi / SCALE))
    return bits


def pick_L(base, sigma, reserve=0, Lmax=16):
    """Select the split point L that minimizes the real signature size.

    The low L bits of each coefficient are stored verbatim, while the high part
    is rANS-coded. We sweep every feasible L and return the one with the
    smallest real bits/coeff, so L is chosen to minimize the signature rather
    than to match any fixed bucket count. 'base' is B0 for z1 or the center
    offset for z_rest; 'reserve' matches the ESC reservation of the emitted
    table."""
    best_L, best_bits = 0, None
    for L in range(0, Lmax + 1):
        nbuckets = ((2 * base) >> L) + 1
        if nbuckets < 2:
            break
        b = real_bits_per_coeff(nbuckets, base, L, sigma, reserve)
        if best_bits is None or b < best_bits - 1e-9:
            best_L, best_bits = L, b
    return best_L

def gauss_cdf(x, sigma):
    return Decimal(0.5) * (1 + Decimal(math.erf(float(x) / (sigma * math.sqrt(2)))))

def bucket_probs(nbuckets, base, L, sigma):
    ps = []
    for h in range(nbuckets):
        lo = h * (1 << L) - base - 0.5
        hi = (h + 1) * (1 << L) - base - 0.5
        ps.append(gauss_cdf(hi, sigma) - gauss_cdf(lo, sigma))
    return ps

def normalize(ps, reserve=0):
    total = SCALE - reserve
    raw = [max(1, int(p * total)) for p in ps]
    s = sum(raw)
    while s != total:
        idx = max(range(len(raw)), key=lambda i: raw[i])
        if s > total:
            if raw[idx] > 1: raw[idx] -= 1; s -= 1
            else: raise RuntimeError("cannot normalize")
        else: raw[idx] += 1; s += 1
    return raw

def emit_mode(level, out):
    m = MODES[level]
    B0, sy = m["B0"], m["sigma_y"]
    s_rest = row_sigma_rest(m)
    cen_s = int(8.5 * s_rest)
    L1 = pick_L(B0, sy)
    LS = pick_L(cen_s, s_rest, reserve=1)
    LB = LS                       # big == small under cut-F
    nz1 = ((2 * B0) >> L1) + 1
    cen_b = cen_s
    nzs = ((2 * cen_s) >> LS) + 1
    nzb = nzs
    f1 = normalize(bucket_probs(nz1, B0, L1, sy))
    fs = normalize(bucket_probs(nzs, cen_s, LS, s_rest), reserve=1); fs.append(1)
    fb = list(fs)                 # identical big table
    out.append(f"#if RHYME_MODE == {level}")
    out.append(f"#define RANS_L1 {L1}")
    out.append(f"#define RANS_LS {LS}")
    out.append(f"#define RANS_LB {LB}")
    out.append(f"#define RANS_NZ1 {nz1}")
    out.append(f"#define RANS_NZS {len(fs)}  /* z_rest rows, ESC at {len(fs)-1} */")
    out.append(f"#define RANS_NZB {len(fb)}  /* == NZS under cut-F */")
    out.append(f"#define RANS_CENTER_S {cen_s}")
    out.append(f"#define RANS_CENTER_B {cen_b}")
    def arr(name, vals):
        lines = [f"static const uint16_t {name}[{len(vals)}] = {{"]
        for i in range(0, len(vals), 12):
            lines.append("    " + ", ".join(map(str, vals[i:i+12])) + ",")
        lines.append("};")
        return lines
    out += arr("rans_freq_z1", f1)
    out += arr("rans_freq_zs", fs)
    out += arr("rans_freq_zb", fb)
    out.append("#endif")
    out.append("")
    def H(sigma, L, nb, base):
        ps = bucket_probs(nb, base, L, sigma)
        return sum(float(p) * (L - math.log2(float(p) + 1e-300)) for p in ps if p > 0)
    n, k, l = m["n"], m["k"], m["l"]; d = k + l - 1
    h1, hs = H(sy, L1, nz1, B0), H(s_rest, LS, nzs, cen_s)
    ct = level // 4
    est = (n * h1 + d * n * hs) / 8 + ct + 4
    print(f"mode {level}: d={d} s_rest={s_rest:.0f} | bits z1={h1:.2f} zrest={hs:.2f} | "
          f"ctilde={ct} -> est sig ~{est:.0f} B")

out = ["/* Auto-generated by gen/rans_tables.py - rANS split-coding freq tables (cut-F) */",
       "#ifndef RHYME_RANS_FREQS_H", "#define RHYME_RANS_FREQS_H",
       "#include <stdint.h>", ""]
for lv in (128, 256, 384, 512):
    emit_mode(lv, out)
out.append("#endif")
open("include/rans_freqs.h", "w").write("\n".join(out))
print("written include/rans_freqs.h")
