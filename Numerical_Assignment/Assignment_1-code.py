import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------
# 1. Inputs
# ----------------------------------
print("\n--- Inputs ---")
n = int(input("Enter number of atomic sites n: "))
E0 = float(input("Enter onsite energy E0: "))
beta = float(input("Enter hopping parameter beta: "))
bin_width = float(input("Enter energy bin width: "))

# ----------------------------------
# 2. Hamiltonian 
# ----------------------------------
H = np.zeros((n, n))

for i in range(n):
    H[i, i] = E0
    if i < n - 1:
        H[i, i + 1] = beta
        H[i + 1, i] = beta

# ----------------------------------
# 3. Eigenvalues and eigenvectors
# ----------------------------------
eigenvalues, eigenvectors = np.linalg.eigh(H)

# ----------------------------------
# 4. Analysing Eigenvalues
# ----------------------------------
num_electrons = n                  # half filling
num_occupied_states = num_electrons // 2

E_min = eigenvalues[0]
E_max = eigenvalues[-1]
E_HOMO = eigenvalues[num_occupied_states - 1]
E_LUMO = eigenvalues[num_occupied_states]
E_F = 0.5 * (E_HOMO + E_LUMO) # Fermi line.

print("\n--- Energy Levels ---")
print(f"Minimum energy        : {E_min:.5f}")
print(f"Maximum energy        : {E_max:.5f}")
print(f"HOMO energy           : {E_HOMO:.5f}")
print(f"LUMO energy           : {E_LUMO:.5f}")
print(f"Fermi energy (Ef)     : {E_F:.5f}")

# ----------------------------------
# 5. Plotting Eigenvalues vs Site 
# ----------------------------------
sites = np.arange(1, n + 1)


plt.figure(figsize=(10, 6))
plt.plot(sites,eigenvalues, label="Eigenvalues vs sites", linewidth=2, color='pink')            # plotting Eigenvalue vs Sites
plt.axhline(E_min, label=f"Minimum energy ({E_min:.5f}eV)", linewidth=2,color='blue')           # Marking Minimum Energy
plt.axhline(E_max, label=f"Maximum energy ({E_max:.5f}eV)", linewidth=2, color='red')           # Marking Maximum Energy 
plt.axhline(E_HOMO, label=f"HOMO ({E_HOMO:.5f}eV)", linewidth=2, color='green')                 # Marking HUMO
plt.axhline(E_LUMO, label=f"LUMO ({E_LUMO:.5f}eV)",linewidth=2, color='orange')                 # Marking LUMO  
plt.axhline(E_F, linestyle='--', linewidth=2, label=f"Fermi level ({E_F:.5f}eVs)", color='red') # Marking the Fermi line


plt.xlabel("Sites")
plt.ylabel("Eigenvalues in eV")
plt.title(f"Eigenstates vs Atomic Site for 1D chain (n={n}, E0={E0}eV, beta={beta}eV, ΔE={bin_width}eV)")
plt.legend()
plt.grid(True)
plt.show()

# ----------------------------------
# 6. Density of States
# ----------------------------------
bins = np.arange(E_min, E_max + bin_width, bin_width)
counts, bin_edges = np.histogram(eigenvalues, bins=bins)

DOS = counts / bin_width
energy_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# ----------------------------------
# 7. Plotting DOS vs Energy
# ----------------------------------
plt.figure(figsize=(8, 6))
plt.bar(
    energy_centers,
    DOS,
    width=bin_width,
    align='center',
    edgecolor='black'
)

plt.axvline(E_F, linestyle='--', linewidth=2, label=f"Fermi level ({E_F:.5f}eV)", color='red')  # Marking the Fermi line

plt.xlabel("Energy in eV")
plt.ylabel("Density of States")
plt.title(f"DOS vs Energy for 1D chain (n={n}, E0={E0}eV, beta={beta}eV, ΔE={bin_width}eV)")
plt.legend()
plt.grid(True)
plt.show()
