#!/usr/bin/env python3
# Generate incomplete-NTT constants for Rhyme (degree-4 base blocks).
# n=512: x^512+1 = prod_{128 blocks} (x^4 - r), r ranges over elements with r^128 = -1 (order 256)
# n=1024: 256 blocks, r^256 = -1 (order 512)
# Layered CT butterflies, 7 (resp. 8) layers, zetas in bit-reversed order, Montgomery R=2^32.

import sys

def is_prime(p):
    if p < 2: return False
    for i in range(2, int(p**0.5)+1):
        if p % i == 0: return False
    return True

def primitive_root(q):
    fact = []
    t = q-1
    d = 2
    while d*d <= t:
        if t % d == 0:
            fact.append(d)
            while t % d == 0: t //= d
        d += 1
    if t > 1: fact.append(t)
    for g in range(2, q):
        if all(pow(g, (q-1)//f, q) != 1 for f in fact):
            return g
    raise RuntimeError

def brv(i, bits):
    return int(format(i, f'0{bits}b')[::-1], 2)

def gen(q, n):
    assert is_prime(q)
    blocks = n // 4              # number of degree-4 blocks
    order = 2 * blocks           # required root order (= 2n/4)
    assert (q - 1) % order == 0, (q, order)
    g = primitive_root(q)
    psi = pow(g, (q-1)//order, q)        # primitive `order`-th root of unity
    assert pow(psi, order, q) == 1 and pow(psi, order//2, q) == q-1
    layers = blocks.bit_length() - 1 + 1  # log2(blocks) ... layers = log2(n) - 2
    L = (n.bit_length()-1) - 2            # log2(n)-2 layers
    # zetas[k] for k=1..blocks-1 in standard CT order: zeta = psi^{brv(k, L)} ... like Kyber:
    # tree: at node k (1-indexed heap), twiddle = psi^{brv(k,L)}; leaves give r values psi^{2*brv+1}? 
    # We'll emulate Kyber exactly: zetas[k] = psi^{brv(k, L)} where psi is primitive 2*blocks-th root.
    zetas = [pow(psi, brv(k, L), q) for k in range(blocks)]
    # base block constants: after full CT recursion, block i (0..blocks-1) is mod (x^4 - zeta_base[i])
    # zeta_base[i] = psi^{2*brv(i,L)+1}? verify empirically below instead of trusting formula.
    return g, psi, L, zetas

def ntt_ref(a, q, n, zetas):
    # Kyber-style iterative CT, stopping at length-4 blocks
    a = list(a)
    k = 1
    length = n // 2
    while length >= 4:
        start = 0
        while start < n:
            zeta = zetas[k]; k += 1
            for j in range(start, start+length):
                t = (zeta * a[j+length]) % q
                a[j+length] = (a[j] - t) % q
                a[j] = (a[j] + t) % q
            start += 2*length
        length //= 2
    return a

def intt_ref(ah, q, n, zetas):
    a = list(ah)
    k = n//4 - 1
    length = 4
    while length <= n//2:
        start = 0
        while start < n:
            zeta = zetas[k]; k -= 1
            for j in range(start, start+length):
                t = a[j]
                a[j] = (t + a[j+length]) % q
                a[j+length] = (zeta * (a[j+length] - t)) % q
            start += 2*length
        length *= 2
    ninv = pow(n//4, q-2, q)   # we divided n/4 times by 2 => multiply by (blocks)^-1? verify below
    return a

def poly_mul_school(a, b, q, n):
    c = [0]*(2*n)
    for i, ai in enumerate(a):
        if ai == 0: continue
        for j, bj in enumerate(b):
            c[i+j] = (c[i+j] + ai*bj) % q
    for i in range(n, 2*n):
        c[i-n] = (c[i-n] - c[i]) % q
    return c[:n]

def base_roots(q, n, zetas):
    # After NTT, block i (i=0..blocks-1) holds a(x) mod (x^4 - r_i).
    # Determine r_i empirically: NTT of x^4 should give r_i in coeff slot 0 of each block.
    blocks = n//4
    a = [0]*n; a[4] = 1
    ah = ntt_ref(a, q, n, zetas)
    r = [ah[4*i] for i in range(blocks)]
    # sanity: block must be exactly (r_i, 0,0,0)
    for i in range(blocks):
        assert ah[4*i+1] == 0 and ah[4*i+2] == 0 and ah[4*i+3] == 0
        assert pow(r[i], blocks, q) == q-1, "r^blocks must be -1"
    return r

def basemul(a, b, r, q):
    c = [0]*4
    for i in range(4):
        for j in range(4):
            k = i+j
            if k < 4: c[k] = (c[k] + a[i]*b[j]) % q
            else: c[k-4] = (c[k-4] + a[i]*b[j]*r) % q
    return c

def full_check(q, n):
    import random
    random.seed(1)
    g, psi, L, zetas = gen(q, n)
    blocks = n//4
    r = base_roots(q, n, zetas)
    a = [random.randrange(q) for _ in range(n)]
    b = [random.randrange(q) for _ in range(n)]
    ah, bh = ntt_ref(a, q, n, zetas), ntt_ref(b, q, n, zetas)
    ch = []
    for i in range(blocks):
        ch += basemul(ah[4*i:4*i+4], bh[4*i:4*i+4], r[i], q)
    # inverse: invert the CT loop properly (GS butterflies with inverse zetas)
    c = gs_intt(ch, q, n, zetas)
    ref = poly_mul_school(a, b, q, n)
    assert c == ref, "NTT-based mul mismatch"
    return g, psi, zetas, r

def gs_intt(ah, q, n, zetas):
    a = list(ah)
    blocks = n//4
    # rebuild forward (length,start)->k mapping
    kmap = {}
    k = 1
    length = n//2
    while length >= 4:
        start = 0
        while start < n:
            kmap[(length, start)] = k; k += 1
            start += 2*length
        length //= 2
    length = 4
    while length <= n//2:
        start = 0
        while start < n:
            zinv = pow(zetas[kmap[(length, start)]], q-2, q)
            for j in range(start, start+length):
                t = a[j]
                a[j] = (t + a[j+length]) % q
                a[j+length] = ((t - a[j+length]) * zinv) % q
            start += 2*length
        length *= 2
    f = pow(blocks, q-2, q)
    return [(x*f) % q for x in a]

if __name__ == "__main__":
    # Build verified Montgomery-domain tables for the cut-F parameter set and
    # emit src/ntt_tables.c.  Forward zetas / izetas / broots scale as plain*R;
    # the final inverse factor intt_f scales as (blocks^-1)*R^2 (matches ntt.c).
    R = 1 << 32
    def brv2(i, bits): return int(format(i, f"0{bits}b")[::-1], 2)
    def prim_root(qq):
        fact = []; t = qq - 1; dd = 2
        while dd * dd <= t:
            if t % dd == 0:
                fact.append(dd)
                while t % dd == 0: t //= dd
            dd += 1
        if t > 1: fact.append(t)
        for gg in range(2, qq):
            if all(pow(gg, (qq - 1) // f, qq) != 1 for f in fact): return gg
        raise RuntimeError
    def center(x, qq):
        x %= qq
        if x > qq // 2: x -= qq
        return x
    def build(qq, nn):
        blocks = nn // 4; order = 2 * blocks
        g = prim_root(qq); psi = pow(g, (qq - 1) // order, qq)
        Lb = blocks.bit_length() - 1
        zfull = [pow(psi, brv2(k, Lb), qq) for k in range(blocks)]
        def fwd(a, z):
            a = [x % qq for x in a]; k = 1; length = nn // 2
            while length >= 4:
                start = 0
                while start < nn:
                    zeta = z[k]; k += 1; j = start
                    while j < start + length:
                        t = (zeta * a[j + length]) % qq
                        a[j + length] = (a[j] - t) % qq; a[j] = (a[j] + t) % qq; j += 1
                    start = j + length
                length >>= 1
            return a
        xa = [0] * nn; xa[4] = 1
        r = [fwd(xa, zfull)[4 * i] for i in range(blocks)]
        kmap = {}; k = 1; length = nn // 2
        while length >= 4:
            start = 0
            while start < nn:
                kmap[(length, start)] = k; k += 1; start += 2 * length
            length //= 2
        iz = [0] * blocks; length = 4
        while length <= nn // 2:
            kb = nn // (2 * length); start = 0
            while start < nn:
                iz[kb + start // (2 * length)] = pow(zfull[kmap[(length, start)]], qq - 2, qq)
                start += 2 * length
            length *= 2
        intf = pow(blocks, qq - 2, qq)
        return zfull, iz, r, intf, blocks
    out = ["/* Auto-generated (gen/ntt_tables.py). Incomplete NTT, degree-4 base blocks. */",
           "#include <stdint.h>", ""]
    for qq, nn in [(3329, 256), (9473, 512), (11777, 512), (18433, 1024)]:
        z, iz, r, intf, blk = build(qq, nn)
        tag = f"{qq}_{nn}"
        def arr(name, vals):
            s = [f"const int32_t rhyme_{name}_{tag}[{len(vals)}] = {{"]
            for i in range(0, len(vals), 8):
                s.append("    " + ", ".join(str(v) for v in vals[i:i + 8]) + ",")
            s.append("};")
            return s
        zm = [center((x * R) % qq, qq) for x in z]
        izm = [center((x * R) % qq, qq) for x in iz]
        rm = [center((x * R) % qq, qq) for x in r]
        intfm = center((intf * R * R) % qq, qq)
        out.append(f"/* q={qq} n={nn}: {blk} blocks (x^4 - r_i) */")
        out += arr("zetas", zm); out += arr("izetas", izm); out += arr("broots", rm)
        out.append(f"const int32_t rhyme_intt_f_{tag} = {intfm};")
        out.append("")
    open("src/ntt_tables.c", "w").write("\n".join(out))
    print("wrote src/ntt_tables.c")
