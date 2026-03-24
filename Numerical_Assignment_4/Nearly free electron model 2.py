import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eigh

# =========================
# Constants (atomic units)
# =========================
hbar = 1.0
m = 1.0
a = 1.0
G0 = 2*np.pi/a

# =========================
# Parameters
# =========================
A = 3.0
k_points = 400

# =========================
# Generate G vectors
# =========================
def generate_G(N):
    n = np.arange(-N, N+1)
    return n * G0

# =========================
# Fourier potential
# =========================
def V_fourier(Gdiff, A):
    if np.isclose(abs(Gdiff), G0):
        return A
    return 0.0

# =========================
# Hamiltonian
# =========================
def build_H(k, G_list, A):
    N = len(G_list)
    H = np.zeros((N, N))

    for i, G in enumerate(G_list):
        for j, Gp in enumerate(G_list):
            if i == j:
                H[i, j] = (hbar**2/(2*m))*(k + G)**2
            H[i, j] += V_fourier(G - Gp, A)
    return H

# =========================
# Band calculation
# =========================
def compute_bands(A, N_G, k_vals):
    G_list = generate_G(N_G)
    bands = []

    for k in k_vals:
        H = build_H(k, G_list, A)
        eigvals, _ = eigh(H)
        bands.append(eigvals)

    return np.array(bands)

# =========================
# k-space grids
# =========================
k_1BZ = np.linspace(-np.pi/a, np.pi/a, k_points)
k_2BZ = np.linspace(-2*np.pi/a, 2*np.pi/a, k_points)

# =========================
# Q1: Band structure
# =========================
N_G = 5
bands_2BZ = compute_bands(A, N_G, k_2BZ)

plt.figure(figsize=(8,6))
for i in range(5):
    plt.plot(k_2BZ, bands_2BZ[:, i], label=f'Band {i+1}')

# Mark BZ boundaries
plt.axvline(np.pi/a, linestyle='--', color='k')
plt.axvline(-np.pi/a, linestyle='--', color='k')
plt.axvline(2*np.pi/a, linestyle='--', color='k')
plt.axvline(-2*np.pi/a, linestyle='--', color='k')

plt.xlabel("k")
plt.ylabel("Energy")
plt.title("Band Structure (1st & 2nd BZ)")
plt.legend()
plt.grid()
plt.savefig("Q1_band_structure.png")
plt.close()

# =========================
# Q2: Density of States
# =========================
bands_1BZ = compute_bands(A, N_G, k_1BZ)
energies = bands_1BZ.flatten()

bins = 300
counts, edges = np.histogram(energies, bins=bins)
dE = edges[1] - edges[0]

DOS = counts / (dE * len(k_1BZ))   # proper normalization
E_centers = 0.5*(edges[:-1] + edges[1:])

plt.figure()
plt.plot(E_centers, DOS)
plt.xlabel("Energy")
plt.ylabel("DOS")
plt.title("Density of States")
plt.grid()
plt.savefig("Q2_DOS.png")
plt.close()

# =========================
# Q3: Convergence vs G
# =========================
plt.figure()
for NG in [2, 3, 4, 5]:
    bands_temp = compute_bands(A, NG, k_1BZ)
    for i in range(2):
        plt.plot(k_1BZ, bands_temp[:, i])

plt.title("Convergence vs Plane Waves")
plt.xlabel("k")
plt.ylabel("Energy")
plt.grid()
plt.savefig("Q3_convergence.png")
plt.close()

# =========================
# Q4: Fermi Energy (1 electron/cell)
# =========================
# Each band holds 2 electrons per k (spin)
# So 1 electron → half filling of first band

Ef = np.max(bands_1BZ[:, 0])

plt.figure()
for i in range(5):
    plt.plot(k_1BZ, bands_1BZ[:, i])

plt.axhline(Ef, linestyle='--', label="Fermi Energy")
plt.legend()
plt.title("Band Structure with Fermi Level")
plt.savefig("Q4_fermi.png")
plt.close()

print("Fermi Energy =", Ef)

# =========================
# Fermi convergence vs k mesh
# =========================
k_mesh_list = [50, 100, 200, 400, 800]
Ef_list = []

for kN in k_mesh_list:
    k_temp = np.linspace(-np.pi/a, np.pi/a, kN)
    bands_temp = compute_bands(A, N_G, k_temp)
    Ef_temp = np.max(bands_temp[:, 0])
    Ef_list.append(Ef_temp)

plt.figure()
plt.plot(k_mesh_list, Ef_list, marker='o')
plt.xlabel("k-mesh")
plt.ylabel("Fermi Energy")
plt.title("Fermi Energy Convergence")
plt.grid()
plt.savefig("Q4_fermi_convergence.png")
plt.close()

# =========================
# Q5: Modified potential
# =========================
def V_modified(Gdiff, A):
    if np.isclose(abs(Gdiff), G0):
        return A
    if np.isclose(abs(Gdiff), 2*G0):
        return A/2
    return 0

def build_H_mod(k, G_list, A):
    N = len(G_list)
    H = np.zeros((N, N))

    for i, G in enumerate(G_list):
        for j, Gp in enumerate(G_list):
            if i == j:
                H[i, j] = (k + G)**2 / 2
            H[i, j] += V_modified(G - Gp, A)
    return H

def compute_bands_mod(A, N_G, k_vals):
    G_list = generate_G(N_G)
    bands = []
    for k in k_vals:
        H = build_H_mod(k, G_list, A)
        eigvals, _ = eigh(H)
        bands.append(eigvals)
    return np.array(bands)

bands_mod = compute_bands_mod(A, N_G, k_1BZ)

plt.figure()
for i in range(5):
    plt.plot(k_1BZ, bands_mod[:, i])

plt.title("Band Structure (Modified Potential)")
plt.savefig("Q5_modified_band.png")
plt.close()

# DOS for modified case
energies_mod = bands_mod.flatten()
counts, edges = np.histogram(energies_mod, bins=300)
DOS_mod = counts / ((edges[1]-edges[0]) * len(k_1BZ))
E_centers = 0.5*(edges[:-1] + edges[1:])

plt.figure()
plt.plot(E_centers, DOS_mod)
plt.title("DOS (Modified Potential)")
plt.savefig("Q5_modified_DOS.png")
plt.close()

# =========================
# Q6: Total energy
# =========================
A_vals = np.linspace(0.5, 5, 12)
E_total = []

for A_val in A_vals:
    bands_temp = compute_bands(A_val, N_G, k_1BZ)
    # half-filled first band
    energy = np.sum(bands_temp[:, 0]) / len(k_1BZ)
    E_total.append(energy)

plt.figure()
plt.plot(A_vals, E_total, marker='o')
plt.xlabel("A")
plt.ylabel("Energy per unit cell")
plt.title("Total Energy vs A")
plt.grid()
plt.savefig("Q6_total_energy.png")
plt.close()
