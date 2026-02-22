"""
===========================================================
Nearly Free Electron Model - Numerical Demonstration
===========================================================

This script numerically demonstrates:

1) show that the bands will have zero slope at the middle and edges of the BZs. Recall in class how we defined 1st, 2nd, 3rd ... Brillouin zones.

Solution: Energy bands have zero slope at:
            - Γ point (k = 0)
            - Brillouin zone boundaries (k = nG/2)


2) Show that the band-gaps between the bands primarily depend on the strength of the harmonics in the decomposition of the periodic potential.:

Solution: Band gaps satisfy:
            Δ_n ≈ 2 |V_{nG}|

All plots are saved as image files.
===========================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eigvalsh
import os

# =========================================================
# SECTION 1: PHYSICAL CONSTANTS (Dimensionless Units)
# =========================================================

hbar = 1.0
m = 1.0
a = 1.0
G = 2*np.pi/a

# Weak periodic potential harmonics
V1 = 0.3     # First harmonic V_G
V2 = 0.15    # Second harmonic V_2G


# =========================================================
# SECTION 2: FOURIER DECOMPOSED PERIODIC POTENTIAL
# =========================================================

def V_fourier(n, V1, V2):
    if abs(n) == 1:
        return V1
    elif abs(n) == 2:
        return V2
    else:
        return 0.0


# =========================================================
# SECTION 3: HAMILTONIAN (Plane Wave Mixing)
# =========================================================

def Hamiltonian(k, V1, V2, N=7):

    n_vals = np.arange(-N, N+1)
    size = len(n_vals)
    H = np.zeros((size, size))

    for i, n in enumerate(n_vals):
        k_n = k + n*G
        H[i,i] = (hbar**2 * k_n**2)/(2*m)

        for j, m_ in enumerate(n_vals):
            if i != j:
                H[i,j] = V_fourier(n-m_, V1, V2)

    return H


# =========================================================
# SECTION 4: COMPUTE BAND STRUCTURE
# =========================================================

def compute_bands(V1, V2, N=7, Nk=2000):

    k_vals = np.linspace(-2*G, 2*G, Nk)
    bands = []

    for k in k_vals:
        H = Hamiltonian(k, V1, V2, N)
        eigvals = eigvalsh(H)
        bands.append(eigvals)

    return k_vals, np.array(bands)


# Create output folder
os.makedirs("figures", exist_ok=True)

k_vals, bands = compute_bands(V1, V2)


# =========================================================
# QUESTION 1: ZERO SLOPE VERIFICATION
# =========================================================

def slope_at(k_target, band):

    idx = np.argmin(np.abs(k_vals - k_target))
    window = 6

    k_local = k_vals[idx-window:idx+window]
    E_local = band[idx-window:idx+window]

    coeff = np.polyfit(k_local, E_local, 2)
    return 2*coeff[0]*k_target + coeff[1]


print("\n========== QUESTION 1 ==========")
for point in [0, G/2, G, 3*G/2]:
    slope = slope_at(point, bands[:,0])
    print(f"Slope at k = {point:.3f} : {slope:.6e}")


# Save full band structure
plt.figure(figsize=(10,7))
for i in range(5):
    plt.plot(k_vals, bands[:,i], linewidth=2)

for n in range(-4,5):
    plt.axvline(n*G/2, linestyle='--', color='black', alpha=0.3)

plt.xlabel("Crystal Momentum k (1/a)")
plt.ylabel("Energy (ħ²/2m)")
plt.title("Band Structure Across Brillouin Zones")
plt.grid(True)
plt.tight_layout()
plt.savefig("figures/1_band_structure.png", dpi=300)
plt.close()


# =========================================================
# QUESTION 2: GAP MEASUREMENT
# =========================================================

def measure_gap(k_target, lower_band, upper_band):
    idx = np.argmin(np.abs(k_vals - k_target))
    gap = bands[idx, upper_band] - bands[idx, lower_band]
    return gap, bands[idx, lower_band], bands[idx, upper_band]


gap1, E_low, E_high = measure_gap(G/2, 0, 1)

print("\n========== QUESTION 2 ==========")
print(f"Gap at k = G/2 : {gap1:.6f}")
print(f"Theoretical 2V1: {2*V1:.6f}")


# Save gap visualization
boundary = G/2
idx = np.argmin(np.abs(k_vals - boundary))
window = 120

plt.figure(figsize=(8,6))
k_zoom = k_vals[idx-window:idx+window]
band1 = bands[idx-window:idx+window, 0]
band2 = bands[idx-window:idx+window, 1]

plt.plot(k_zoom, band1, label="Band 1")
plt.plot(k_zoom, band2, label="Band 2")

plt.axvline(boundary, linestyle='--')
plt.hlines(E_low, boundary-0.4, boundary+0.4, linestyles='--')
plt.hlines(E_high, boundary-0.4, boundary+0.4, linestyles='--')

plt.text(boundary-1.2, E_high,
         f"Δ ≈ {gap1:.4f}\n2V1 = {2*V1:.4f}",
         bbox=dict(facecolor='white', alpha=0.8))

plt.xlabel("Crystal Momentum k (1/a)")
plt.ylabel("Energy (ħ²/2m)")
plt.title("Band Gap at First BZ Boundary")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("figures/2_band_gap.png", dpi=300)
plt.close()


# =========================================================
# GAP SCALING TEST
# =========================================================

V_values = np.linspace(0.05, 0.5, 15)
gaps = []

for V in V_values:
    k_temp, bands_temp = compute_bands(V, 0.0)
    idx = np.argmin(np.abs(k_temp - G/2))
    gaps.append(bands_temp[idx,1] - bands_temp[idx,0])

gaps = np.array(gaps)
coeff = np.polyfit(V_values, gaps, 1)
slope = coeff[0]

print("\n========== GAP SCALING ==========")
print(f"Fitted slope ≈ {slope:.6f}")
print("Theoretical slope = 2")


plt.figure(figsize=(8,6))
plt.plot(V_values, gaps, 'o', label="Numerical Δ")
plt.plot(V_values, 2*V_values, '--', label="Theory Δ = 2V")
plt.plot(V_values, slope*V_values,
         label=f"Fit slope ≈ {slope:.3f}")

plt.xlabel("Fourier Harmonic V₁ (energy)")
plt.ylabel("Band Gap (energy)")
plt.title("Gap Scaling with Potential Harmonic")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("figures/3_gap_scaling.png", dpi=300)
plt.close()


"""
Results:

========== QUESTION 1 ==========
Slope at k = 0.000 : 3.146025e-09
Slope at k = 3.142 : -1.007320e-02
Slope at k = 6.283 : 2.842286e-09
Slope at k = 9.425 : -8.281337e-03

========== QUESTION 2 ==========
Gap at k = G/2  : 0.596139
Theoretical 2V1 : 0.600000

========== GAP SCALING TEST ==========
Fitted slope from numerical data ≈ 1.993686
Theoretical slope = 2.000000


"""
