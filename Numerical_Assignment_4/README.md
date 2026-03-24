# README — Numerical Assignment 2

## Nearly Free Electron Model (Intro to CMP, 2026)

---

## 1. Overview

This project implements a **numerical solution of the Nearly Free Electron (NFE) model** in one dimension using a **plane-wave (Fourier) basis**. The goal is to compute and analyze:

* Electronic band structure
* Density of states (DOS)
* Convergence with respect to basis size
* Fermi energy and its convergence
* Effect of modifying the periodic potential
* Approximate total energy per unit cell

The implementation follows the standard **Bloch theorem + Fourier expansion approach** used in solid state physics.

---

## 2. Physical Model

We consider a 1D periodic system with lattice constant ( a ), and a weak periodic potential:

[
V(x) = 2A \cos(G_0 x), \quad G_0 = \frac{2\pi}{a}
]

In Fourier space:

* ( V_{\pm G_0} = A )
* All other Fourier components are zero

---

### Hamiltonian in plane-wave basis

[
H_{G,G'}(k) = \frac{\hbar^2}{2m}(k + G)^2 \delta_{G,G'} + V_{G-G'}
]

This is diagonalized for each ( k )-point.

---

## 3. Numerical Method

### Steps:

1. Choose a finite set of reciprocal lattice vectors:
   [
   G = n G_0, \quad n = -N, ..., N
   ]

2. For each ( k ):

   * Construct Hamiltonian matrix ( H(k) )
   * Diagonalize to obtain eigenvalues ( E_n(k) )

3. Repeat over a grid of ( k )-points

---

## 4. Units and Assumptions

* Atomic units used:
  [
  \hbar = 1, \quad m = 1, \quad a = 1
  ]

* Spin degeneracy included:

  * Each band holds **2 electrons per k**

* System is:

  * Non-interacting
  * No electron-electron interactions
  * No self-consistency (no Poisson solver)

---

## 5. Description of Outputs

---

### **Q1_band_structure.png**

* Band structure plotted across:

  * First Brillouin Zone: ( [-\pi/a, \pi/a] )
  * Second Brillouin Zone: ( [-2\pi/a, 2\pi/a] )

* Features:

  * Multiple bands (at least 5)
  * Brillouin zone boundaries marked

**Physics:**

* Band gaps open at zone boundaries due to Bragg reflection

---

### **Q2_DOS.png**

* Density of states computed using histogram method

[
DOS(E) = \frac{\text{counts}}{\Delta E \cdot N_k}
]

**Physics:**

* Peaks correspond to **van Hove singularities**
* Flat bands → high DOS

---

### **Q3_convergence.png**

* Band structure plotted for increasing number of plane waves

**Physics:**

* Increasing basis size improves accuracy
* Convergence indicates sufficient Fourier truncation

---

### **Q4_fermi.png**

* Band structure with Fermi energy marked

**Fermi Energy Determination:**

* 1 electron per unit cell
* With spin degeneracy:

  * First band is half-filled

[
E_F = \max_k E_1(k)
]

---

### **Q4_fermi_convergence.png**

* Fermi energy vs k-mesh density

**Physics:**

* Converges as k-sampling improves
* Demonstrates numerical stability

---

### **Q5_modified_band.png**

* Band structure after adding extra Fourier component:

  * ( V_{\pm 2G_0} = A/2 )

**Physics:**

* Additional harmonics modify band gaps
* More complex dispersion

---

### **Q5_modified_DOS.png**

* DOS for modified potential

**Physics:**

* Redistribution of states due to modified periodicity

---

### **Q6_total_energy.png**

* Total energy per unit cell vs potential strength ( A )

[
E_{\text{total}} = \frac{1}{N_k} \sum_k E_1(k)
]

(only occupied states included)

**Physics:**

* Shows how periodic potential affects ground-state energy
* This is a **non-interacting approximation**

---

## 6. Limitations

This model is simplified:

* No electron-electron interaction
* No exchange-correlation effects
* No self-consistent potential

A more accurate method would:

* Solve Poisson equation for charge density
* Iterate to self-consistency

---

## 7. How to Run

1. Install dependencies:

   ```
   pip install numpy matplotlib
   ```

2. Run the script:

   ```
   python filename.py
   ```

3. Output plots will be saved automatically in the working directory.

---

## 8. Key Learning Outcomes

This assignment demonstrates:

* How periodic potentials create band structures
* Origin of band gaps in solids
* Numerical diagonalization of Hamiltonians
* Role of basis truncation
* Relationship between dispersion and DOS

---

## 9. Suggested Extensions

* Implement Gaussian broadening for DOS
* Extend to 2D lattice
* Add self-consistent field (Poisson solver)
* Compare with tight-binding model

---

## 10. Author’s Note

This implementation prioritizes **physical correctness and numerical clarity** over optimization. The results are sufficient for conceptual understanding and coursework submission.

---
