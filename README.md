# Imaginary Time QMC for H2SQ Model

This repository implements an **Imaginary Time Quantum Monte Carlo (QMC)** simulation for the **H2SQ model**, which involves five distinct interaction types. The code employs a hybrid approach combining the **Swendsen-Wang (SW) cluster algorithm** along interlayer and Trotter layers (\( J_3 \), \( K \)), and the **heatbath algorithm** for other interactions (\( J_1 \), \( J_2 \)). A novel loop model was developed, enabling efficient and exact computations at low temperatures and external fields.

---

## Features

- **Hybrid QMC Framework:**
  - SW algorithm for interlayer and Trotter layer interactions (\( J_3 \), \( K \)).
  - Heatbath algorithm for in-plane interactions (\( J_1 \), \( J_2 \)).
- **Loop Model for Low Temperatures:** A specialized loop algorithm for small temperature and field regimes, ensuring exact results.
- **Multi-Interaction Support:** Handles complex interactions characteristic of the H2SQ model.
- **Scalable:** Adaptable to various lattice sizes and Trotter discretizations.
- **Observables Calculation:** Computes thermodynamic quantities and field-dependent observables.

---

## Getting Started

### **Prerequisites**

Ensure you have the following:

- **C++ Compiler** (C++17 or later, e.g., GCC, Clang)
- **CMake** (optional for build automation)
- **Python 3.7+** (for data analysis and plotting):
  - **NumPy**
  - **Matplotlib**

Install Python dependencies with:

```bash
pip install numpy matplotlib
