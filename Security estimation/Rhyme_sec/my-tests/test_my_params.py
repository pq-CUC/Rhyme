import sys
import math
from sage.all import oo, log, RR

# Replace with the path to your local lattice-estimator environment
sys.path.append('/home/sage/project/estimator') 
from estimator import *

# Explicitly import DiscreteGaussian and CenteredBinomial
from estimator.nd import DiscreteGaussian, CenteredBinomial 

print("=" * 70)
print(" FINAL: Security Assessment for V3 Scheme (Rhyme / Const 3 Architecture)")
print(" (Features: Exact l-1 s-Split | mod q MSIS | Dynamic L_2 bound via Chi-square)")
print("=" * 70)

# ==============================================================================
# V3 Parameter Configuration (Rhyme) —— [Perfectly Decoupled Interface Version]
# ==============================================================================
v3_scheme_params = {
    "Rhyme-V3-128": {
        "n_ring": 1024, 
        "k": 2, 
        "l": 3, 
        "q": 18433, 
        
        # [Core Interface 1]: Equivalent standard deviation of a single coefficient for the large polynomial F
        "sigma_F": 34.5640,  
        
        # [Core Interface 2]: Global secure anti-forgery bound
        "L2_bound":  30935.16
    }
}

for name, params in v3_scheme_params.items():
    print(f"\n--- Evaluating parameter set: {name} ---\n")

    # ==============================================================================
    # A. Security evaluation for key recovery attack (MLWE -> LWE)
    # ==============================================================================
    print(">>> (A) Key Recovery Attack (MLWE -> LWE)")

    lwe_n = params["n_ring"] * (params["l"] - 1)
    lwe_m = params["n_ring"] * params["k"]

    # 1. Secret s (corresponding to small polynomial f) uses a Centered Binomial Distribution.
    # (Note: eta=6 is used here. If your scheme uses a different eta, modify this directly.)
    dist_binary_s = CenteredBinomial(eta=6, n=lwe_n)
    
    # 2. Error e (!!! Modification: Uniformly using a Centered Binomial Distribution to 
    # calculate the security lower bound for a purely small distribution !!!)
    dist_binary_e = CenteredBinomial(eta=6, n=lwe_m)

    mlwe_params = LWE.Parameters(
        n=lwe_n,
        q=params["q"],
        Xs=dist_binary_s,    
        Xe=dist_binary_e,    # <--- Modification: Use purely small distribution
        m=lwe_m,
        tag=f"{name}-KeyRecovery"
    )

    print("LWE equivalent parameters:", mlwe_params)

    # --- Classical security evaluation ---
    print("\n--- [Classical] Key Recovery Evaluation Results ---")
    lwe_results_classical = LWE.estimate(mlwe_params, quiet=True, red_cost_model=RC.MATZOV, red_shape_model="cn11")

    min_lwe_classical = oo
    for attack, cost in lwe_results_classical.items():
        print(f"Attack '{attack}': {cost!r}")
        if cost is not None and cost.get("rop") is not None:
            current_rop = cost.get("rop", oo)
            if current_rop != oo and current_rop < min_lwe_classical:
                min_lwe_classical = current_rop

    if min_lwe_classical != oo:
        print(f"==> [Classical] Key Recovery Security Level: {float(log(min_lwe_classical, 2)):.2f} bits")
    else:
        print("==> [Classical] Key Recovery Security Level: Not determined")

    # --- Quantum security evaluation ---
    print("\n--- [Quantum] Key Recovery Evaluation Results ---")
    lwe_results_quantum = LWE.estimate(mlwe_params, quiet=True, red_cost_model=RC.ChaLoy21, red_shape_model="cn11")

    min_lwe_quantum = oo
    for attack, cost in lwe_results_quantum.items():
        print(f"Attack '{attack}': {cost!r}")
        if cost is not None and cost.get("rop") is not None:
            current_rop = cost.get("rop", oo)
            if current_rop != oo and current_rop < min_lwe_quantum:
                min_lwe_quantum = current_rop

    if min_lwe_quantum != oo:
        print(f"==> [Quantum] Key Recovery Security Level: {float(log(min_lwe_quantum, 2)):.2f} bits\n")
    else:
        print("==> [Quantum] Key Recovery Security Level: Not determined\n")


    # ==============================================================================
    # B. Security evaluation for forgery attack (MSIS -> SIS)
    # ==============================================================================
    print(">>> (B) Forgery Attack (MSIS -> SIS)")

    sis_n = params["n_ring"] * params["k"]
    sis_m = params["n_ring"] * (params["k"] + params["l"])
    
    sis_bound = params["L2_bound"]

    msis_params = SIS.Parameters(
        n=sis_n,
        q=params["q"],
        length_bound=sis_bound,
        m=sis_m,
        norm=2, 
        tag=f"{name}-Forgery"
    )

    print("SIS equivalent parameters:", msis_params)

    # --- Classical security evaluation ---
    print("\n--- [Classical] Forgery Resistance Evaluation Results ---")
    sis_results_classical = SIS.estimate(msis_params, quiet=True, red_cost_model=RC.MATZOV)

    min_sis_classical = oo
    for attack, cost in sis_results_classical.items():
        print(f"Attack '{attack}': {cost!r}")
        if cost is not None and cost.get("rop") is not None:
            current_rop = cost.get("rop", oo)
            if current_rop != oo and current_rop < min_sis_classical:
                min_sis_classical = cost.get("rop", oo)

    if min_sis_classical != oo:
        print(f"==> [Classical] Forgery Resistance Security Level: {float(log(min_sis_classical, 2)):.2f} bits")
    else:
        print("==> [Classical] Forgery Resistance Security Level: Not determined")

    # --- Quantum security evaluation ---
    print("\n--- [Quantum] Forgery Resistance Evaluation Results ---")
    sis_results_quantum = SIS.estimate(msis_params, quiet=True, red_cost_model=RC.ChaLoy21)

    min_sis_quantum = oo
    for attack, cost in sis_results_quantum.items():
        print(f"Attack '{attack}': {cost!r}")
        if cost is not None and cost.get("rop") is not None:
            current_rop = cost.get("rop", oo)
            if current_rop != oo and current_rop < min_sis_quantum:
                min_sis_quantum = cost.get("rop", oo)

    if min_sis_quantum != oo:
        print(f"==> [Quantum] Forgery Resistance Security Level: {float(log(min_sis_quantum, 2)):.2f} bits\n")
    else:
        print("==> [Quantum] Forgery Resistance Security Level: Not determined\n")


print("=" * 70)
print(" All parameter sets evaluated")
print("=" * 70)

# ==============================================================================
# C. Security evaluation for Gram Matrix Factorization (Module-LIP -> uSVP)
# ==============================================================================
print(">>> (C) Gram Factorization Attack (Module-LIP -> Unusual-SVP)")

d_matrix = v3_scheme_params["Rhyme-V3-128"]["k"] + v3_scheme_params["Rhyme-V3-128"]["l"] - 1
lip_d = d_matrix * v3_scheme_params["Rhyme-V3-128"]["n_ring"] 

target_norm = 1.0

print(f"LIP Equivalent Dimension (d*n) = {lip_d}")
print(f"Target Unusual Vector Length: {target_norm}")

def estimate_lip_beta(d, target_len):
    """Estimate the required BKZ blocksize (beta) to find a vector of target_len in dimension d."""
    for beta in range(50, d):
        delta_beta = (beta / (2 * math.pi * math.e)) ** (1.0 / (2 * beta - 1))
        lhs = target_len * math.sqrt(beta / d)
        rhs = delta_beta ** (2 * beta - d + 1)
        if lhs <= rhs:
            return beta
    return oo

lip_beta = estimate_lip_beta(lip_d, target_norm)

if lip_beta != oo:
    print(f"==> Required BKZ Blocksize (beta): {lip_beta}")
    lip_classical_bits = 0.292 * lip_beta
    lip_quantum_bits = 0.265 * lip_beta
    
    print(f"==> [Classical] Gram Factorization (LIP) Security Level: {lip_classical_bits:.2f} bits")
    print(f"==> [Quantum] Gram Factorization (LIP) Security Level: {lip_quantum_bits:.2f} bits\n")
else:
    print("==> Required BKZ Blocksize (beta): Exceeds dimension (Secure)")
    print("==> Gram Factorization (LIP) Security Level: Not determined (Extremely High)\n")