import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Physical parameters (Set to 1 for simplicity)
# -----------------------------
hbar = 1.0 
m_e = 1.0
a = 1.0

prefactor = hbar**2 / (2*m_e)

# k-range
k = np.linspace(-6*np.pi/a, 6*np.pi/a, 3000)

# Normalized crystal momentum
k_norm = k / (np.pi/a)

# -----------------------------
# Q1: Periodicity = a
# -----------------------------
plt.figure(figsize=(8, 6))

for m in range(-4, 5):
    E_k = prefactor * (k + 2*np.pi*m/a)**2
    plt.plot(k_norm, E_k, label=f"m = {m}")

# Brillouin zone boundaries
for n in range(-6, 7):
    plt.axvline(x=n, color='k', linestyle='--', linewidth=0.6)

plt.xlabel(r"Normalized crystal momentum $k/(\pi/a)$")
plt.ylabel(r"Energy $E_k$")
plt.title("Sommerfeld Model: Band Structure (Periodicity = a)")
plt.legend()
plt.grid(True)

plt.savefig("Ek_periodicity_a_normalized.png", dpi=300)
plt.close()

# -----------------------------
# Q2: Periodicity = 2a
# -----------------------------
plt.figure(figsize=(8, 6))

for m in range(-8, 9):
    E_k = prefactor * (k + np.pi*m/a)**2
    plt.plot(k_norm, E_k)

# Brillouin zone boundaries
for n in range(-6, 7):
    plt.axvline(x=n, color='k', linestyle='--', linewidth=0.6)

plt.xlabel(r"Normalized crystal momentum $k/(\pi/a)$")
plt.ylabel(r"Energy $E_k$")
plt.title("Sommerfeld Model: Band Structure (Periodicity = 2a)")
plt.grid(True)

plt.savefig("Ek_periodicity_2a_normalized.png", dpi=300)
plt.close()



