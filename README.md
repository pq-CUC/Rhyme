# Rhyme: A Fiat-Shamir Lattice-based Signature with 3C Sampling

This repository contains the reference implementation of the **Rhyme** digital signature scheme, as presented in the paper "Rhyme: A Fiat-Shamir Lattice-based Signature with 3C Sampling".

The core innovation of Rhyme is the **Convolution of Constant Counts Sampling (3C Sampling)** method. This technique allows for the generation of the masking vector `y` without using rejection sampling, thereby mitigating secret key leakage through side-channel attacks.

## Directory Structure

This repository is organized as follows:

- `reference implemention/ref/`: The core source code for the Rhyme signature scheme.
- `security estimation/`: A copy of the [Lattice Estimator](https://github.com/malb/lattice-estimator) tool, used for security analysis of the Rhyme parameter sets. Please see the `README.md` inside this directory for original author and license information.
- `my-tests/`: Contains the specific scripts used to evaluate Rhyme's security with the Lattice Estimator.

## Build and Run

The project uses CMake for building.

### Prerequisites

- CMake (version 3.18 or later)
- A C compiler (e.g., GCC, Clang)
- OpenSSL (for random number generation)

### Building the Project

1. **Navigate to the implementation directory:**
    
    ```bash
    cd "reference implemention/ref"
    ```
2. **Create a build directory:**
    
    ```bash
    cmake -S . -B build
    ```
3. **Compile the code:**
    
    ```bash
    cmake --build build --clean-first
    ```
    
    The compiled binaries will be located in the `reference implemention/ref/build/bin/` directory.

## Security Analysis

The security of the Rhyme parameter sets against known lattice attacks has been analyzed using the Lattice Estimator. The analysis script is `my-tests/test_my_params.py`.

## License

The code in the `reference implemention/ref` directory is licensed under the MIT License. See the [LICENSE](LICENSE) file for details. The `security estimation` directory contains third-party code under its own license.

