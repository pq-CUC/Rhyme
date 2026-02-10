# --- Code required for the solution ---
import sys
import math
# Add the directory containing the estimator library to Python's search path
sys.path.append('/lattice-estimator')
# --- Code required for the solution ---
from sage.all import oo, log
# 1. Import all functionalities from the estimator library
from estimator import *
from estimator.reduction import ADPS16

print("="*60)
print(" Security Assessment for Rhyme Scheme (Binary {0, 1} distribution)")
print("="*60)

# 2. Define the three parameter sets for Rhyme
rhyme_params = {
    "Rhyme-128": {"n_ring": 256, "k": 2, "l": 3, "B": 30, "q": 257},
    "Rhyme-192": {"n_ring": 256, "k": 3, "l": 4, "B": 42, "q": 257},
    "Rhyme-256": {"n_ring": 256, "k": 4, "l": 5, "B": 56, "q": 257}
}

# 3. Loop through each parameter set and evaluate it
for name, params in rhyme_params.items():
    print(f"\n--- Evaluating parameter set: {name} ---\n")

    # A. Security evaluation for key recovery attack (MLWE -> LWE)
    print(">>> (A) Key Recovery Attack (MLWE -> LWE)")

    lwe_n = params["n_ring"] * params["k"]
    lwe_m = params["n_ring"] * params["l"]

    # --- Code modification start ---
    # Change to Uniform distribution over {0, 1}
    # This corresponds to Binary-LWE where s, e <- {0, 1}
    # Mean = 0.5, Stddev = 0.5
    dist = ND.Uniform(0, 1)
    # --- Code modification end ---

    mlwe_params = LWE.Parameters(
        n=lwe_n,
        q=params["q"],
        Xs=dist,
        Xe=dist,
        m=lwe_n,
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


    # B. Security evaluation for forgery attack (MSIS -> SIS)
    # (Note: MSIS part is generally independent of the secret distribution s, depending mostly on dimensions and bound B)
    print(">>> (B) Forgery Attack (MSIS -> SIS)")

    sis_n = params["n_ring"] * params["k"]
    sis_m = params["n_ring"] * (params["k"] + params["l"])
    sis_bound = params["B"]

    msis_params = SIS.Parameters(
        n=sis_n,
        q=params["q"],
        length_bound=sis_bound,
        m=sis_m,
        norm=oo,
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


print("="*60)
print(" All parameter sets evaluated")
print("="*60)