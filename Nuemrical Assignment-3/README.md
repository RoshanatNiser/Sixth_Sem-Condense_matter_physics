# Nearly Free Electron Model – Numerical Assignment

## Objective

Using the argument of plane-wave mixing under a weak periodic potential, we numerically demonstrate:

1. Energy bands have zero slope at:
   - Γ point (k = 0)
   - Brillouin zone boundaries (k = nG/2)

2. Band gaps depend on the strength of Fourier harmonics of the periodic potential:
   
   Δₙ ≈ 2 |VₙG|

---

# Theoretical Background

## 1. Nearly Free Electron Model

Unperturbed dispersion:

E₀(k) = ħ²k² / 2m

The periodic potential is written as:

V(x) = Σ VₙG e^{i nG x}

This couples plane waves:

|k⟩ ↔ |k − nG⟩

---

## 2. Zero Slope at BZ Center and Boundaries

Degeneracy condition:

E₀(k) = E₀(k − nG)

This gives:

k = nG/2

At these points, the Hamiltonian reduces to a 2×2 form:

| E₀  VₙG |
| VₙG E₀  |

Eigenvalues:

E± = E₀ ± |VₙG|

Near degeneracy:

E±(q) = E₀ ± √[(αq)² + |VₙG|²]

Derivative at q = 0:

dE/dk = 0

Numerically, slopes at k = 0, G/2, G, 3G/2 are ≈ 0.

---

## 3. Band Gap Formation

Energy splitting:

Δₙ = 2 |VₙG|

Thus each Brillouin zone boundary corresponds to a specific harmonic:

- k = G/2 → gap depends on V_G
- k = G → gap depends on V_2G
- etc.

---

# Numerical Results

## 1. Band Structure

File: figures/band_structure.png

- Shows multiple Brillouin zones.
- Flat extrema at Γ and zone boundaries.

## 2. Gap Visualization

File: figures/band_gap.png

- Clear avoided crossing at k = G/2.
- Horizontal lines mark gap.
- Numerical Δ ≈ 2V₁.

## 3. Gap Scaling

File: figures/gap_scaling.png

- Δ plotted vs V₁.
- Linear fit slope ≈ 2.
- Confirms perturbation theory.

---

# Conclusion

The numerical simulation confirms:

1. Bands have zero slope at Γ and Bragg planes.
2. Band gaps scale linearly with Fourier harmonics.
3. Δₙ ≈ 2|VₙG| holds in weak potential regime.

This matches the perturbative plane-wave mixing argument of the Nearly Free Electron model.
