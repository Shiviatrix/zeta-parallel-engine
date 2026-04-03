# Riemann Zeta Parallel Engine (MPFR + OpenMP)

This repository contains the high-performance, massively parallel engine used for the high-altitude evaluation of the Riemann zeta function $\zeta(s)$, as described in the manuscript: 

**"A Massively Parallel Architecture for Arbitrary-Precision Evaluation of the Riemann Zeta Function: Performance and Heuristic Limits"** (SIURO Submission).

## Overview

The engine utilizes the **Riemann-Siegel formula** for pointwise evaluation at extreme altitudes ($T > 10^{30}$). Key features include:
- **Thread-Local Private Workspace (TLPW)**: Lock-free parallel summation using OpenMP.
- **Arbitrary Precision**: Configurable precision (default 500 bits) using the GNU MPFR library.
- **Memory Efficient**: $O(\log T)$ space complexity, bypassing the hardware limitations of FFT-based pre-computation methods.

## Prerequisites

- **GCC 12.1.0+** (with OpenMP 5.0 support)
- **GNU MPFR 4.1.0+**
- **GNU GMP 6.2.1+**

## Installation (Magus/HPC Cluster)

To compile the engine with optimization flags for many-core architectures:

```bash
g++-12 -O3 -fopenmp fast_10_37_mpfr.cpp -o zeta_engine -lmpfr -lgmp
```

## Usage

### Simple Evaluation
To evaluate a single point or a small interval:
```bash
./zeta_engine <T_START> <T_END> <STEP_SIZE>
```

### Distributed Probing
Use the provided `zeta_worker.py` to manage larger batches across multiple nodes:
```bash
python3 zeta_worker.py --start 1e37 --nodes 32
```

## Academic Context

This code was developed as part of an undergraduate research study at the **Shiv Nadar Institution of Eminence (SNIoE)**. It is intended for investigating the statistical distribution of zeros at extreme heights where deterministic leading-term heuristics begin to break down.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
