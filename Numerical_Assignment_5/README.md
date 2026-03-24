---
## Tight Binding Model Assignment

This project numerically implements and analyzes tight-binding models corresponding to the three problems in the assignment. The code is written in Python using `NumPy` for numerical computation and `Matplotlib` for visualization.

---

# 1. Structure of the Code

The script is organized into functions corresponding to each question:

* `ladder_bands()` → Question 1
* `ladder_diagonal_case_a()` → Question 2(a)
* `ladder_diagonal_case_b()` → Question 2(b)
* `square_lattice_band()` → Question 3(a) band contour
* `fermi_surface()` → Fermi surface extraction
* `square_lattice_with_diagonal()` → Question 3(b)

All plots are saved in the `plots/` directory.

---

# 2. Question 1: Two-Leg Ladder Model

## Physical Model

The system consists of a two-leg ladder with:

* nearest-neighbor hopping along each leg: ( t )
* inter-leg (rung) hopping: ( t' )

Using a two-site basis per unit cell, the Hamiltonian in momentum space is:

[
H(k) =
\begin{pmatrix}
-2t\cos k & -t' \
-t' & -2t\cos k
\end{pmatrix}
]

Diagonalization yields two bands:

[
E_\pm(k) = -2t\cos k \pm t'
]

## Numerical Implementation

The code:

* discretizes ( k \in [-\pi, \pi] )
* computes both bands for different values of ( t'/t )
* plots them on the same figure

## Physical Interpretation

* The splitting between bands is controlled by ( t' )
* A gap opens when the bands no longer overlap
* This identifies the metal–insulator transition

---

# 3. Question 2: Effect of Diagonal Hopping

## (a) Same-direction diagonal hopping

The off-diagonal term becomes momentum-dependent:

[
t' \rightarrow t' + 2t''\cos k
]

This modifies the band structure asymmetrically.

### Code Behavior

* Computes:
  [
  E_\pm(k) = -2t\cos k \pm (t' + 2t''\cos k)
  ]
* Plots for several values of ( t'' )

### Physical Effect

* Band splitting becomes ( k )-dependent
* Transition point shifts

---

## (b) Alternating diagonal hopping

The off-diagonal term becomes complex:

[
\sim t' + 2i t'' \sin k
]

Leading to:

[
E_\pm(k) = -2t\cos k \pm \sqrt{t'^2 + (2t''\sin k)^2}
]

### Code Behavior

* Computes square-root expression for dispersion
* Plots modified bands

### Physical Effect

* Gap becomes momentum-dependent
* More robust gap opening compared to case (a)

---

# 4. Question 3: Square Lattice

## (a) Nearest-neighbor hopping

The dispersion is:

[
E(k_x, k_y) = -2t(\cos k_x + \cos k_y)
]

### Numerical Method

1. Construct a 2D grid in the first Brillouin zone:
   [
   k_x, k_y \in [-\pi, \pi]
   ]

2. Compute energy at each point

3. Flatten and sort energies to determine Fermi energy:

   * For a given filling ( n ), the number of occupied states is proportional to ( n )
   * Spin degeneracy is included via factor of 2

4. Extract Fermi surface using contour:
   [
   E(k_x, k_y) = E_F
   ]

### Output

* Band contour plot
* Fermi surfaces for selected fillings

---

## (b) Inclusion of diagonal hopping

The dispersion becomes:

[
E(k_x,k_y) = -2t(\cos k_x + \cos k_y) - 4t_d \cos k_x \cos k_y
]

### Numerical Implementation

* Same grid as part (a)
* Modify energy using additional term
* Repeat Fermi surface extraction

### Physical Interpretation

* Breaks particle-hole symmetry
* Alters topology of Fermi surface
* Produces anisotropic deformation depending on sign of ( t_d )

---

# 5. Fermi Surface Calculation

The function `fermi_surface()` performs:

1. Flattening of energy grid
2. Sorting energies
3. Determination of Fermi energy index:
   [
   \text{index} = \left\lfloor \frac{n}{2} \cdot N \right\rfloor
   ]
4. Extraction of contour at that energy

A boundary check ensures the index does not exceed array size.

---

# 6. Output Files

All plots are saved in the `plots/` directory with names indicating:

* question number
* subpart (a/b)
* parameter values
* filling (for Fermi surfaces)

Example:

```
Q1_band_structure.png
Q2a_diagonal.png
Q3a_fermi_n_1.00.png
Q3b_td_-1.0_fermi_n_1.50.png
```

---

# 7. How to Run

Install required packages:

```
pip install numpy matplotlib
```

Run the script:

```
python your_script_name.py
```

All results will be generated automatically.

---
