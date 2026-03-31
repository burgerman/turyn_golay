# Turyn-Golay Sequence Search (Parallel Solution: MPI Optimized)

This project implements a high-performance parallel search for **Turyn-Golay sequences** of length $4n+2$. It uses the Turyn construction to find four sequences (A, B, C, D) that satisfy the Golay property: the sum of their Non-periodic Autocorrelation Functions (NPAF) is zero for all non-zero shifts.

## Background & Applications

Golay sequences (and their Turyn-constructed variants) are sets of sequences whose combined autocorrelation side lobes cancel out perfectly. This "perfect autocorrelation" property makes them invaluable in several high-tech fields:

- **Radar & Sonar:** Used for pulse compression to achieve high-resolution range detection. The zero-sidelobe property ensures that small targets are not masked by the "ringing" or side lobes of larger nearby targets.
- **Telecommunications (OFDM):** In systems like Wi-Fi and 5G, Golay sequences are used to reduce the **Peak-to-Average Power Ratio (PAPR)**, allowing for more efficient power amplifier usage and longer battery life in mobile devices.
- **Optical Coherence Tomography (OCT):** Improves the signal-to-noise ratio in medical imaging, allowing for clearer scans of biological tissues.
- **Cryptography & Coding Theory:** Used in the construction of error-correcting codes and robust synchronization headers in data packets.

Finding these sequences is a computationally intensive combinatorial problem. As the length $n$ increases, the search space grows exponentially, requiring advanced algorithmic pruning and parallel computing techniques.

## Algorithm Features

### 1. Turyn Construction
The algorithm searches for four sequences:
- **A, B**: Length $n+1$
- **C, D**: Length $n$

The total length of the resulting Golay sequence is $(n+1) + (n+1) + n + n = 4n+2$.

### 2. NSOKS (Sum of Squares) with Memoization
Before searching for sequences, the algorithm finds all quadruples of integers $(k_{11}, r_{11}, p_{11}, q_{11})$ such that:
$$k_{11}^2 + r_{11}^2 + p_{11}^2 + q_{11}^2 = 4n + 2$$
This is implemented using a recursive `nsoks` (n-squares-of-k-sums) algorithm with **memoization** to avoid redundant partitioning calculations.

### 3. Multi-Stage Filtering
To handle the massive search space, the implementation applies multiple theoretical constraints:
- **Step 2 Condition 3:** Prunes candidate sequences based on parity and specific element-wise sum constraints derived from Turyn's theorem.
- **Step 2 Condition 5:** Ensures the combined Periodic Autocorrelation Function (PAF) contribution of the four sequences is zero.

## Key Optimizations

### 🚀 Meet-in-the-Middle (Hash Map Lookup)
The most significant optimization is the transformation of the search from a nested loop ($O(N^2)$) into a near-linear lookup ($O(N)$).
- We precompute all valid pairs of (K, R) sequences and store their combined PAF contribution in an `std::unordered_map`.
- For every valid (P, Q) pair, we simply look up the negative of its PAF contribution in the map.
- This effectively "meets in the middle," drastically reducing the number of required comparisons.

### 📊 PAF Precomputation
Instead of calculating autocorrelation values repeatedly in the inner loops, the algorithm precomputes the **PAF contribution vector** for every candidate sequence immediately after generation. This turns an $O(m^2)$ operation into a simple vector addition during the search phase.

### 🌐 MPI Parallelization
The search space is distributed across multiple processes using **MPI (Message Passing Interface)**.
- Each rank processes a subset of the candidate `PQ` pairs.
- **Early Exit:** Ranks periodically synchronize using `MPI_Allreduce`. If any rank finds a valid solution, all ranks terminate immediately to save computational resources.

## Requirements

- **C++17** compatible compiler.
- **MPI** library (e.g., OpenMPI or MPICH).

## Build and Run

### Compilation
Using a standard C++ compiler with MPI support:
```bash
mpic++ -O3 main.cpp utilities.cpp -o turyn_golay -I.
```

### Running
To search for sequences of length $4n+2$ (e.g., $n=7$ results in length 30):
```bash
mpirun -n <number_of_processes> ./turyn_golay <n>
```
Example for $n=15$:
```bash
mpirun -n 8 ./turyn_golay 15
```

## Implementation Details

- **`main.cpp`**: Entry point, handles MPI initialization and parameter parsing.
- **`utilities.h`**: Interface for sequence generation, verification, and optimized data structures (`PAFVector`, `SequenceCandidate`).
- **`utilities.cpp`**: Core implementation of the Turyn construction, `nsoks` partitioning, and the optimized Meet-in-the-Middle search logic.
