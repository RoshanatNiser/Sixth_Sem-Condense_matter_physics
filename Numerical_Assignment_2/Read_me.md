# Sommerfeld Model: Bloch Hamiltonian in the Free-Electron Limit

## Assignment Description

This repository contains a numerical study of the Sommerfeld (free-electron) model
derived from Bloch’s theorem in one dimension.

Starting from the Bloch Hamiltonian

H_k u_k = E_k u_k,

the amplitude of the periodic potential is taken to zero while the lattice periodicity
is retained. The goal is to visualize the resulting energy bands and understand their
dependence on lattice periodicity.

---

## Theory Background

For a one-dimensional lattice with periodicity `a`, Bloch’s theorem gives

ψ_k(x) = exp(i k x) u_k(x),

where u_k(x) is periodic with the lattice.

In the free-electron (Sommerfeld) limit, the Bloch Hamiltonian reduces to

H_k = (1/2m) ( -iħ d/dx + ħk )².

Choosing a plane-wave basis consistent with lattice periodicity,

u_m(x) = exp(i 2πm x / a),

the energy eigenvalues are

E_k(m) = (ħ² / 2m) ( k + 2πm/a )².

Each integer `m` corresponds to a distinct energy band.

---

## Numerical Tasks

### Q1: Periodicity = a

Energy bands are computed and plotted for

k ∈ [ -6π/a , +6π/a ]

using

E_k(m) = (ħ² / 2m) ( k + 2πm/a )².

The integer `m` labels the band index.

### Q2: Periodicity = 2a

The lattice periodicity is doubled to `2a`, leading to

E_k(m) = (ħ² / 2m) ( k + πm/a )².

The band structure is plotted over the same k-range and compared with Q1.

---

## Code Structure

- `sommerfeld_band_structure.py`
  - Defines physical parameters
  - Computes energy eigenvalues for different band indices
  - Generates and saves plots automatically

- Output plots:
  - `Ek_periodicity_a.png`
  - `Ek_periodicity_2a.png`

---

## Normalized Momentum and Brillouin Zone Boundaries

The crystal momentum axis is normalized as k/(π/a), rendering it dimensionless.
With this choice, Brillouin zone boundaries occur at integer values of the
normalized momentum.

Vertical dashed lines are drawn in the plots at these integer positions to explicitly
indicate zone boundaries. This representation makes zone folding and the
effect of changing lattice periodicity immediately visible.

---
## How to Use

1. Clone or download this repository.
2. Ensure Python 3 is installed with the following packages:
   - numpy
   - matplotlib
3. Run the script:
   ```bash
   python sommerfeld_band_structure.py

4. The generated plots will be saved in the same directory.

