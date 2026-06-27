# Rhyme: A Fiat-Shamir Lattice-based Signature

---

### Overview

The core innovation of Rhyme is its "cut-F" compressed-unimodular construction combined with a widened matrix dimension and HAWK-style Lattice Isomorphism Problem (LIP) hiding. In lattice-based Fiat-Shamir signatures, preventing the signature from leaking the secret key is a central challenge.

Rhyme constructs the signature vector using a hybrid approach: it employs a deterministic, lightweight rejection sampler (Algorithm 5) solely for the first component to decouple it from the shared mask, while the remaining components rely on HAWK-style LIP hiding. The hiding distribution depends on the secret only through the Gram matrix of the short rows ($G_0 = B_{short} \cdot B_{short}^*$). 

Furthermore, widening the unimodular column dimension to $D$ and discarding the long NTRU-solved rows completely defeats GA-MLWE distinguishers by spreading each Gram entry over a D-term superposition. This design ensures that the final signature distribution is inherently independent of the secret key without relying on the expensive, full-dimensional rejection sampling used in schemes like Dilithium.

### Directory Structure

This repository is organized as follows:

-   `Reference implementation/ref/`: The core source code for the Rhyme signature scheme.
-   `Security estimation/`: A copy of the [Lattice Estimator](https://github.com/malb/lattice-estimator) tool, used for security analysis of the Rhyme parameter sets. Please see the `README.md` inside this directory for original author and license information.
-   `my-tests/`: Contains the specific scripts used to evaluate Rhyme's security with the Lattice Estimator.

### Build and Run

The project uses CMake for building.

#### Prerequisites

-   CMake (version 3.18 or later)
-   A C compiler (e.g., GCC, Clang) with AES-NI / SSSE3 support (`-maes -mssse3`)

#### Building the Project

1.  **Navigate to the implementation directory:**
    ```bash
    cd "Reference implementation/ref"
    ```

2.  **Create a build directory:**
    *(Append `-DCMAKE_C_FLAGS="-DUSE_DP1_CHECK"` to enable the fast DP1 keygen optimization)*
    ```bash
    cmake -S . -B build
    ```

3.  **Compile the code:**
    ```bash
    cmake --build build --clean-first
    ```

    The compiled binaries will be located in the `Reference implementation/ref/build/bin/` directory.

### Security Analysis

The classical and quantum hardness of each parameter set (Forgery, Key-Recovery, and LIP tasks) has been analyzed using the Lattice Estimator. Specifically, classical and quantum security have been evaluated using the **MATZOV (2022)** and **ChaLoy21** models, respectively, as the current standard benchmarks. The analysis script is `my-tests/test_my_params.py`.

### License

The code in the `Reference implementation/ref` directory is licensed under the MIT License. See the `LICENSE` file for details. The `Security estimation` directory contains third-party code under its own license.