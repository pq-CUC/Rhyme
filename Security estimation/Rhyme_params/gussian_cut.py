#!/usr/bin/env python3
"""
Rhyme Parameter Derivation (F-Cut Version) — Derives all parameters from core inputs.
L2 bounds are compared using three methods: Standard Chi-Square / Generalized Chi-Square / Monte Carlo.
(Updated: Decoupled the strict k+l binding, introduced D_MATRIX to prevent algebraic distinguisher attacks.)
"""

import math
import numpy as np
from scipy.stats import chi2
from scipy.signal import fftconvolve
from scipy.special import erfcinv

# ==================== Core Parameters (Sole Inputs) ====================
TARGET_Z1_SUCCESS_RATE = 0.5
TARGET_L2_SUCCESS_RATE = 0.7

K = 2
L = 3
N_DIM = 512
ACC_RATE_BOTTOM = 0.99
NTT_DEGREE = 4

SECURITY_LEVEL_CHALLENGE = 256

LAMBDA_SECURITY = 256
QS_SIGNATURES = 2 ** 64
EPSILON = 1.0 / math.sqrt(QS_SIGNATURES * LAMBDA_SECURITY)

ETA_F = 4

# ----- [NEW] Introduce D_MATRIX to satisfy the security constraint D > 4k + 1 -----
if K == 2 and L == 3:
    D_MATRIX = 10
elif K == 3 and L == 4:
    D_MATRIX = 14
else:
    D_MATRIX = 4 * K + 2  # Fallback logic
# -----------------------------------------------------------

# Internal Settings
GRAM_TRIALS = 30
MC_SAMPLES = 20000
RNG_SEED = 12342
# ============================================================

SIGMA_F_SQ = ETA_F / 2.0


# ---------- Basic Utilities ----------
def is_prime(num):
    """Check if a number is prime."""
    if num <= 1: return False
    if num <= 3: return True
    if num % 2 == 0 or num % 3 == 0: return False
    i = 5
    while i * i <= num:
        if num % i == 0 or num % (i + 2) == 0: return False
        i += 6
    return True


def find_q_incomplete_ntt(linf, n, deg, count=20):
    """Find incomplete NTT-friendly primes."""
    order = (2 * n) // deg
    c = max(1, math.ceil((2 * linf) / order))
    qs = []
    while len(qs) < count:
        q = order * c + 1
        if q > 2 * linf and is_prime(q):
            qs.append(q)
        c += 1
    return qs  # Returns a list; the first element is the minimum satisfying value


def calculate_min_tau(n, bits):
    """Calculate minimum tau for required security bits."""
    for tau in range(1, n + 1):
        ent = math.log2(math.comb(n, tau)) + tau
        if ent >= bits: return tau, ent
    raise ValueError("n is too small")


# ---------- True Acceptance Rate for z1 Rejection Sampling (Theorem 3.1) ----------
def calculate_z1_true_rate(b0, sigma_y, sigma_base, n):
    """Calculate the true rejection sampling rate for z1."""
    sigma_tar = 2.0 * sigma_base
    R = math.ceil(5.0 * sigma_tar)
    B0 = int(round(b0))
    Vs, ps = {}, {}

    for c in (0, 1):
        V = [v for v in range(-R, R + 1) if (v - c) % 2 == 0]
        w = [math.exp(-(v * v) / (2.0 * sigma_tar * sigma_tar)) for v in V]
        S = sum(w)
        Vs[c], ps[c] = V, [x / S for x in w]

    M = 0.0
    for c in (0, 1):
        for y in range(-(B0 + R), B0 + R + 1):
            ry = math.exp(-(y * y) / (2.0 * sigma_y * sigma_y))
            s = sum(ps[c][k] * math.exp(-((y + Vs[c][k]) ** 2) / (2.0 * sigma_y * sigma_y))
                    for k in range(len(Vs[c]))) / ry
            M = max(M, s)

    num = sum(math.exp(-(z * z) / (2.0 * sigma_y * sigma_y)) for z in range(-B0, B0 + 1))
    den = sum(math.exp(-(y * y) / (2.0 * sigma_y * sigma_y)) for y in range(-(B0 + R), B0 + R + 1))
    return ((num / den) / M) ** n


def find_inflection_point_fast(sigma_y, sigma_base, n):
    """Quickly find the inflection point for optimal parameters."""
    best_b0, best_rate = None, -1.0
    for b0 in range(int(sigma_y * 2.0), int(sigma_y * 5.0)):
        rate = calculate_z1_true_rate(b0, sigma_y, sigma_base, n)
        if rate > best_rate:
            best_rate, best_b0 = rate, b0
    return best_b0, best_rate


def search_optimal_parameters_z1(sigma_base, n):
    """Search for the optimal z1 parameters based on target success rate."""
    for sigma_y in range(50, 2000, 2):
        b0, z1_rate = find_inflection_point_fast(sigma_y, sigma_base, n)
        if z1_rate >= TARGET_Z1_SUCCESS_RATE:
            return float(sigma_y), float(b0), z1_rate
    raise ValueError("Failed to find parameters satisfying z1!")


# ---------- Ring / CBD / G0 Spectra (F-Cut Version) ----------
def cbd(eta, size, rng):
    """Centered Binomial Distribution sampling."""
    a = rng.integers(0, 2, (size, eta)).sum(1)
    b = rng.integers(0, 2, (size, eta)).sum(1)
    return (a - b).astype(float)


def negacyclic(c):
    """Generate a negacyclic matrix."""
    n = len(c)
    M = np.empty((n, n))
    for k in range(n):
        col = np.empty(n)
        col[k:] = c[:n - k]
        col[:k] = -c[n - k:]
        M[:, k] = col
    return M


# [MODIFIED] Added D parameter. Replaced d+1 with D.
def build_B_short(d, D, n, eta, rng):
    """
    Construct short matrix B:
    First d rows × D columns, independent CBD_eta; the long row F is cut.
    """
    return np.block([[negacyclic(cbd(eta, n, rng)) for _ in range(D)] for _ in range(d)])


# [MODIFIED] Passed through D parameter.
def sample_G0_spectra(d, D, n, eta, trials, rng):
    """Sample the spectra of G0 over multiple trials."""
    spectra = []
    for _ in range(trials):
        Bs = build_B_short(d, D, n, eta, rng)
        ev = np.clip(np.linalg.eigvalsh(Bs @ Bs.T), 0.0, None)
        spectra.append(ev)
    return spectra


# ---------- L2 Method 1: Standard Chi-Square ----------
# [MODIFIED] Added D parameter, replacing (d+1) in degrees of freedom and expected variance.
def l2_standard_chi2(d, D, n, sigma_y, sigma_base, eta, target):
    """Calculate L2 bound using the standard Chi-Square distribution."""
    v_marg = (sigma_base ** 2) * (1.0 + 4.0 * D * n * (eta / 2.0))
    scale_1, df_1 = sigma_y ** 2, n
    scale_2, df_2 = v_marg, d * n

    mu = df_1 * scale_1 + df_2 * scale_2
    sd = math.sqrt(2 * df_1 * scale_1 ** 2 + 2 * df_2 * scale_2 ** 2)
    grid = np.linspace(0, mu + 8 * sd, 400000)

    p1 = np.nan_to_num(chi2.pdf(grid / scale_1, df=df_1) / scale_1)
    p1 /= p1.sum()
    p2 = np.nan_to_num(chi2.pdf(grid / scale_2, df=df_2) / scale_2)
    p2 /= p2.sum()

    tot = fftconvolve(p1, p2)[:len(grid)]
    tot = np.clip(tot, 0, None)
    tot /= tot.sum()
    cdf = np.cumsum(tot)

    ti = int(np.argmax(cdf >= target))
    return math.sqrt(grid[ti]), float(cdf[ti]), v_marg


# ---------- L2 Method 2: Generalized Chi-Square (True Spectra) ----------
def weighted_chi2_pdf(weights, grid):
    """Helper for Generalized Chi-Square: Compute weighted Chi-Square PDF."""
    nb = min(80, max(4, int(np.sqrt(len(weights)))))
    edges = np.quantile(weights, np.linspace(0, 1, nb + 1))
    edges[0] -= 1e-9
    edges[-1] += 1e-9

    total = None
    for b in range(nb):
        m = (weights > edges[b]) & (weights <= edges[b + 1])
        cnt = int(m.sum())
        if cnt == 0: continue

        wbar = float(weights[m].mean())
        pb = np.nan_to_num(chi2.pdf(grid / wbar, df=cnt) / wbar)
        s = pb.sum()
        if s <= 0: continue
        pb /= s

        total = pb if total is None else fftconvolve(total, pb)[:len(grid)]
        total = np.clip(total, 0, None)
        ss = total.sum()
        if ss > 0: total /= ss
    return total


def l2_generalized_chi2(d, n, sigma_y, sigma_base, spectrum, target):
    """Calculate L2 bound using the Generalized Chi-Square distribution."""
    w = (sigma_base ** 2) * (1.0 + 4.0 * spectrum)
    mu_b = w.sum()
    var_b = (2 * w ** 2).sum()

    scale_1, df_1 = sigma_y ** 2, n
    mu_1, var_1 = df_1 * scale_1, 2 * df_1 * scale_1 ** 2
    grid = np.linspace(0, (mu_1 + mu_b) + 8 * math.sqrt(var_1 + var_b), 400000)

    pb = weighted_chi2_pdf(w, grid)
    p1 = np.nan_to_num(chi2.pdf(grid / scale_1, df=df_1) / scale_1)
    p1 /= p1.sum()

    tot = fftconvolve(p1, pb)[:len(grid)]
    tot = np.clip(tot, 0, None)
    tot /= tot.sum()
    cdf = np.cumsum(tot)

    ti = int(np.argmax(cdf >= target))
    return math.sqrt(grid[ti])


# ---------- L2 Method 3: Monte Carlo ----------
def l2_monte_carlo(d, n, sigma_y, sigma_base, spectrum, target, n_samples, rng):
    """Calculate L2 bound using Monte Carlo simulations."""
    w = (sigma_base ** 2) * (1.0 + 4.0 * spectrum)
    dim_b = len(w)
    norms = np.empty(n_samples)
    for i in range(n_samples):
        z1_sq = (sigma_y ** 2) * (rng.standard_normal(n) ** 2).sum()
        zb_sq = (w * rng.standard_normal(dim_b) ** 2).sum()
        norms[i] = z1_sq + zb_sq
    return math.sqrt(np.quantile(norms, target))


# ---------- L_inf (Union Bound) ----------
def linf_union(d, n, v_marg, acc):
    """Calculate L_infinity limit using the union bound."""
    per = (1.0 - acc) / (d * n)
    return math.ceil(math.sqrt(2 * v_marg) * erfcinv(per)) + 1


# ==================== Main Workflow ====================
def derive():
    """Main parameter derivation routine."""
    rng = np.random.default_rng(RNG_SEED)
    m = K + L
    d = m - 1
    d_lattice = d * N_DIM
    base_eta = math.sqrt(math.log(2 * d_lattice * (1 + 1 / EPSILON)) / math.pi)
    sigma_base = base_eta / math.sqrt(2 * math.pi)

    print("=" * 72)
    print(" Rhyme Parameter Derivation (F-Cut Version) - Security Columns D Fix")
    print("=" * 72)

    # [0] Challenge Space
    tau, ent = calculate_min_tau(N_DIM, SECURITY_LEVEL_CHALLENGE)
    print(f" [0] Challenge Space: tau = {tau}, entropy = {ent:.2f} bits")
    # [MODIFIED] Added D_MATRIX to the printed info
    print(
        f"     d = K+L-1 = {d}, D (Matrix Columns) = {D_MATRIX}, N = {N_DIM}, ETA = {ETA_F}, sigma_base = {sigma_base:.6f}")

    # [1] z1: Search for sigma_y and B0 to reach target acceptance rate
    print("-" * 72)
    sigma_y, b0, z1_rate = search_optimal_parameters_z1(sigma_base, N_DIM)
    print(f" [1] z1 Rejection Sampling (Derived via search)")
    print(f"     sigma_y = {sigma_y}, B0 = {b0}, True Acceptance Rate = {z1_rate * 100:.2f}%")

    # [2] G0 Spectra (F-Cut)
    print("-" * 72)
    # [MODIFIED] Updated print info with D_MATRIX
    print(f" [2] Covariance Derivation for Residual z_bottom (F-Cut: {d} rows × {D_MATRIX} columns, F has been cut)")
    sigma_e_sq = sigma_base ** 2  # e_bottom variance (I term)
    sigma_sign_sq = sigma_base ** 2  # sigma for 4G term
    print(f"     Sampled {GRAM_TRIALS} B_short matrices to estimate G0 = B_short·B_short* spectra")

    # [MODIFIED] Passed D_MATRIX
    spectra = sample_G0_spectra(d, D_MATRIX, N_DIM, ETA_F, GRAM_TRIALS, rng)
    eigs_all = np.concatenate(spectra)
    cond = eigs_all.max() / max(eigs_all.min(), 1e-12)

    # Variance for each eigen-direction (for Generalized Chi-Sq/Monte Carlo) vs Marginal variance (diagonal approx)
    w_dir = (sigma_base ** 2) * (1.0 + 4.0 * spectra[0])

    # [MODIFIED] Replaced d+1 with D_MATRIX when calculating marginal variance
    v_marg = (sigma_base ** 2) * (1.0 + 4.0 * D_MATRIX * N_DIM * (ETA_F / 2.0))

    print(
        f"     Base smooth Gaussian std dev sigma_base = {sigma_base:.4f}  <-- Aligns with HAWK theoretical lower bound")
    print(f"     e_bottom variance (I term)       = {sigma_e_sq:.4f}")
    print(f"     4G term sigma^2                  = {sigma_sign_sq:.4f}")
    print(
        f"     G0 eigenvalues min/mean/max      = {eigs_all.min():.3f} / {eigs_all.mean():.3f} / {eigs_all.max():.3f}")
    print(f"     Condition Number (Anisotropy)    = {cond:.1f}")
    print(f"     >> Marginal Variance (Isotropic Approx) = {v_marg:.4f} (Equivalent std = {math.sqrt(v_marg):.4f}) <<")
    print(f"     >> Eigen-direction Variance (True Anisotropy) min/max = {w_dir.min():.2f} / {w_dir.max():.2f} <<")

    # [3] Three Methods for L2
    print("-" * 72)
    print(f" [3] Anti-forgery Bound and L2 Weight Bridge (Target Success Rate {TARGET_L2_SUCCESS_RATE * 100:.0f}%)")

    # [MODIFIED] Passed D_MATRIX
    L2_std, ach_std, v_marg = l2_standard_chi2(d, D_MATRIX, N_DIM, sigma_y, sigma_base, ETA_F, TARGET_L2_SUCCESS_RATE)

    gen_list = [l2_generalized_chi2(d, N_DIM, sigma_y, sigma_base, sp, TARGET_L2_SUCCESS_RATE) for sp in spectra]
    L2_gen = float(np.mean(gen_list))

    mc_list = [l2_monte_carlo(d, N_DIM, sigma_y, sigma_base, sp, TARGET_L2_SUCCESS_RATE, MC_SAMPLES, rng) for sp in
               spectra]
    L2_mc = float(np.mean(mc_list))

    print(
        f"     Method 1 Standard Chi-Square   : L2 = {L2_std:.2f}  (L2SQ={int(L2_std ** 2)}, Boundary Success Rate {ach_std * 100:.4f}%)")
    print(f"     Method 2 Generalized Chi-Square: L2 = {L2_gen:.2f}  (L2SQ={int(L2_gen ** 2)})")
    print(f"     Method 3 Monte Carlo           : L2 = {L2_mc:.2f}  (L2SQ={int(L2_mc ** 2)})")
    print(
        f"     Standard vs MC Diff = {abs(L2_std - L2_mc) / L2_mc * 100:.2f}%   Generalized vs MC Diff = {abs(L2_gen - L2_mc) / L2_mc * 100:.3f}%")
    print(f"     >> Recommended L2 Weight Bridge (Adopting Generalized Chi-Square/MC) = {L2_mc:.2f} <<")

    # [4] L_inf + Q
    print("-" * 72)
    print(f" [4] L_inf Truncation Bound and Modulus Q")
    p_single = ACC_RATE_BOTTOM ** (1 / (d * N_DIM))  # Original script style (for reference)
    B1 = linf_union(d, N_DIM, v_marg, ACC_RATE_BOTTOM)
    max_linf = max(int(round(b0)), B1)
    q_list = find_q_incomplete_ntt(max_linf, N_DIM, NTT_DEGREE, count=20)

    print(f"     Bottom Height Limit B1 (Union Bound) = {B1}")
    print(f"     Global Safe Maximum Bound L_inf = max(B0, B1) = {max_linf}")
    print(f"     >> Top 10 Q values satisfying conditions (Ascending): <<")
    for i, q in enumerate(q_list, 1):
        if i > 10: break
        print(f"        {i:2d}. Q = {q}")
    print(f"     Minimum usable Q = {q_list[0]}; Larger values can be chosen to satisfy SIS security margin")
    print("=" * 72)


if __name__ == "__main__":
    derive()