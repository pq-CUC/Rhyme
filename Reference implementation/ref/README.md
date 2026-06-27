# Rhyme — Reference Implementation

Reference implementation of **Rhyme**, a Fiat–Shamir lattice-based signature scheme. This implementation features the "cut-F" compressed-unimodular design combined with HAWK-style LIP hiding.

The signature algorithm relies solely on integer arithmetic over bounded domains. Exact integer products are computed via a dedicated Montgomery 31-bit NTT, ensuring constant-time operations for core polynomial multiplications. The implementation utilizes a deterministic rejection sampler with fixed-point arithmetic for the first coordinate, eliminating the need for floating-point sampling, while using an OpenSSL-free AES256-CTR-DRBG for symmetric expansions.

### Layout

-   `include/`: Header files, auto-generated parameter tables, and rANS frequencies.
-   `src/`: Core implementations including polynomial arithmetic (`poly.c`), NTT (`ntt.c`), symmetric primitives (`symmetric-shake.c`, `aes256ctr.c`), the sampler (`sampler.c`), and entropy coding (`encoding.c`).
-   `src/keygen/`: Pure-C unimodular-basis solver for key generation. Includes the fast DP1 "truncated tail columns" coprimality check.
-   `test/`: Functional tests and golden NTT vectors.
-   `kat/`: NIST PQCgenKAT harness and deterministic KAT vectors (OpenSSL-free).
-   `benchmark/`: Performance benchmarking tools (`speed.c`, `cpucycles`).

### Build Instructions

You can build the implementation using either CMake or plain GNU Make.

#### Prerequisites

-   CMake (3.18 or later) or GNU make
-   A C compiler (GCC or Clang) with AES-NI / SSSE3 support (`-maes -mssse3`)

#### Option 1: Using Make (Recommended for testing)

```bash
make            # Builds tests, KAT generators, and benchmarks for all modes
make check      # Runs NTT golden tests and functional tests
```

**Fast Keygen Optimization:**
You can enable the DP1 "truncated tail columns" check for key generation, which speeds up the unimodular extension check by ~4-5x:
```bash
make DP1=1
```

#### Option 2: Using CMake

To compile with CMake and properly enable the DP1 fast keygen check, you must pass the compiler definition directly during configuration:

```bash
cmake -S . -B build -DCMAKE_C_FLAGS="-DUSE_DP1_CHECK"
cmake --build build
ctest --test-dir build
```

Compiled binaries (tests, benchmarks, and KAT generators) will be placed in the `build/bin/` directory.

### Running Tests & Benchmarks

Once built, you can run the executables directly from the `build/bin/` directory:

-   **Functional Test:** `./build/bin/test_rhyme256`
-   **Performance Benchmark:** `./build/bin/speed_rhyme256`
-   **KAT Generation:** `./build/bin/kat_rhyme256`

*(Replace `256` with `128`, `384`, or `512` to test different security levels).*

### License

This reference implementation is licensed under the MIT License. See the `LICENSE` file in the root directory for details.