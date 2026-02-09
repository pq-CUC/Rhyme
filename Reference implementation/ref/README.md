# Rhyme

A lattice-based digital signature scheme from the paper "Rhyme: A Fiat-Shamir Lattice-based Signature with 3C Sampling".

This scheme is based on the Fiat-Shamir transform and introduces a new sampling method called "Convolution of Constant Counts Sampling (3C Sampling)" to prevent secret key leakage, avoiding expensive full-dimensional rejection sampling.

This repository contains the reference implementation for the Rhyme signature scheme.

---

## Build

build with cmake
```
cmake -S ./ -B build
cmake --build build --clean-first
```


