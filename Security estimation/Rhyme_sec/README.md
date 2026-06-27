# Lattice Estimator (Third-Party Tool)

**This directory contains a copy of the Lattice Estimator, an open-source tool for estimating the security of lattice-based cryptographic schemes.**

This tool is used by the parent project (Rhyme Signature Scheme) for security analysis.

- **Original Repository:** [https://github.com/malb/lattice-estimator](https://github.com/malb/lattice-estimator)
- **Original Authors & Contributors:** Please see the `README.rst` file in this directory.
- **License:** The Lattice Estimator is licensed under the LGPLv3+ license.

All credit for the code within this directory goes to its original authors. It is included here for convenience to reproduce the security analysis results for the Rhyme signature scheme.

## Contributions to this Directory

The test script located at `my-tests/test_my_params.py` in the parent directory was written to specifically evaluate the security parameters of the Rhyme scheme using this estimator.

