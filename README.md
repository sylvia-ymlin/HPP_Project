# High‑Performance N‑Body Simulation (HPP_Project)

Gravitational N‑body simulation suffers from O(N²) complexity in naïve pairwise force computation. 
This project implements the **Barnes–Hut treecode** to reduce that cost, combined with parallelization and memory optimizations.

## Key Aspects
- Implementation in C with **OpenMP** for parallel force computation  
- Two integrators: **symplectic Euler** and **Velocity Verlet**  
- Use of **K‑means clustering** to improve spatial locality and load balancing  
- Optimizations like precomputing reciprocals, minimizing branching  
- Performance and accuracy benchmarks via accompanying python scripts

**Code** | **Report**  
[Code](https://github.com/sylvia-ymlin/HPP_Project/tree/main/project) | [Project Report (PDF)](https://github.com/sylvia-ymlin/HPP_Project/blob/main/project/report.pdf)

---

## Build & Run

```bash
cd project
make
```

### Run simulation

```bash
./simulate -n 1000 -t 500 -dt 0.001 -theta 0.5 -method verlet
```

**Parameters:**  
- `-n` — number of particles  
- `-t` — number of time steps  
- `-dt` — time increment  
- `-theta` — Barnes–Hut opening angle  
- `-method` — `euler` or `verlet`  

### Performance analysis & plotting

```bash
cd python
python code_performance.py
python code_complexity.py
```

### Validation tests

```bash
cd test
./run_tests.sh
```

This compares outputs against reference data in `ref_output_data`.

---
