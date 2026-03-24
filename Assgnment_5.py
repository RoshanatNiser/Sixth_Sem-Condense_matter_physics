import numpy as np
import matplotlib.pyplot as plt
import os

# ==============================
# SETUP
# ==============================
plt.rcParams["figure.figsize"] = (8,6)
plt.rcParams["font.size"] = 12

os.makedirs("plots", exist_ok=True)

# ==============================
# Q1: TWO-LEG LADDER
# ==============================
def ladder_bands(t=1.0, t_prime_list=[0.5, 1.5, 2.5]):
    k = np.linspace(-np.pi, np.pi, 500)

    plt.figure()
    for t_prime in t_prime_list:
        E_plus = -2*t*np.cos(k) + t_prime
        E_minus = -2*t*np.cos(k) - t_prime

        plt.plot(k, E_plus, label=f"E+ (t'/t={t_prime/t:.2f})")
        plt.plot(k, E_minus, linestyle='--',
                 label=f"E- (t'/t={t_prime/t:.2f})")

    plt.title("Q1: Two-leg ladder band structure")
    plt.xlabel("k")
    plt.ylabel("Energy")
    plt.legend()
    plt.grid()

    plt.savefig("plots/Q1_band_structure.png", dpi=300, bbox_inches='tight')
    plt.close()


# ==============================
# Q2(a): DIAGONAL HOPPING
# ==============================
def ladder_diagonal_case_a(t=1.0, t_prime=1.0, tpp_list=[-0.5, 0.0, 0.5]):
    k = np.linspace(-np.pi, np.pi, 500)

    plt.figure()
    for tpp in tpp_list:
        E_plus = -2*t*np.cos(k) + (t_prime + 2*tpp*np.cos(k))
        E_minus = -2*t*np.cos(k) - (t_prime + 2*tpp*np.cos(k))

        plt.plot(k, E_plus, label=f"E+ (t''={tpp})")
        plt.plot(k, E_minus, linestyle='--',
                 label=f"E- (t''={tpp})")

    plt.title("Q2(a): Diagonal hopping effect")
    plt.xlabel("k")
    plt.ylabel("Energy")
    plt.legend()
    plt.grid()

    plt.savefig("plots/Q2a_diagonal.png", dpi=300, bbox_inches='tight')
    plt.close()


# ==============================
# Q2(b): DIAGONAL HOPPING
# ==============================
def ladder_diagonal_case_b(t=1.0, t_prime=1.0, tpp_list=[-0.5, 0.0, 0.5]):
    k = np.linspace(-np.pi, np.pi, 500)

    plt.figure()
    for tpp in tpp_list:
        E_plus = -2*t*np.cos(k) + np.sqrt(t_prime**2 + (2*tpp*np.sin(k))**2)
        E_minus = -2*t*np.cos(k) - np.sqrt(t_prime**2 + (2*tpp*np.sin(k))**2)

        plt.plot(k, E_plus, label=f"E+ (t''={tpp})")
        plt.plot(k, E_minus, linestyle='--',
                 label=f"E- (t''={tpp})")

    plt.title("Q2(b): Diagonal hopping effect")
    plt.xlabel("k")
    plt.ylabel("Energy")
    plt.legend()
    plt.grid()

    plt.savefig("plots/Q2b_diagonal.png", dpi=300, bbox_inches='tight')
    plt.close()


# ==============================
# Q3(a): SQUARE LATTICE
# ==============================
def square_lattice_band(t=-2.0, N=200):
    kx = np.linspace(-np.pi, np.pi, N)
    ky = np.linspace(-np.pi, np.pi, N)
    KX, KY = np.meshgrid(kx, ky)

    E = -2*t*(np.cos(KX) + np.cos(KY))

    plt.figure()
    plt.contourf(KX, KY, E, levels=50)
    plt.colorbar(label="Energy (eV)")
    plt.title("Q3(a): Band contour")
    plt.xlabel("kx")
    plt.ylabel("ky")

    plt.savefig("plots/Q3a_band_contour.png", dpi=300, bbox_inches='tight')
    plt.close()

    return E, KX, KY


# ==============================
# FERMI SURFACE (FIXED)
# ==============================
def fermi_surface(E, KX, KY, fillings, prefix):
    flat_E = E.flatten()
    sorted_E = np.sort(flat_E)
    N = len(sorted_E)

    for n in fillings:
        idx = int(n * N / 2)
        if idx >= N:
            idx = N - 1

        Ef = sorted_E[idx]

        plt.figure()
        plt.contour(KX, KY, E, levels=[Ef])
        plt.title(f"{prefix}: Fermi surface (n={n:.2f})")
        plt.xlabel("kx")
        plt.ylabel("ky")
        plt.grid()

        filename = f"plots/{prefix}_fermi_n_{n:.2f}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"Saved: {filename}")


# ==============================
# Q3(b): WITH DIAGONAL HOPPING
# ==============================
def square_lattice_with_diagonal(t=-2.0, td=1.0, N=200):
    kx = np.linspace(-np.pi, np.pi, N)
    ky = np.linspace(-np.pi, np.pi, N)
    KX, KY = np.meshgrid(kx, ky)

    E = -2*t*(np.cos(KX) + np.cos(KY)) - 4*td*np.cos(KX)*np.cos(KY)

    plt.figure()
    plt.contourf(KX, KY, E, levels=50)
    plt.colorbar(label="Energy (eV)")
    plt.title(f"Q3(b): Band (td={td})")
    plt.xlabel("kx")
    plt.ylabel("ky")

    plt.savefig(f"plots/Q3b_td_{td}_band.png", dpi=300, bbox_inches='tight')
    plt.close()

    return E, KX, KY


# ==============================
# MAIN
# ==============================
if __name__ == "__main__":

    # Q1
    ladder_bands()

    # Q2
    ladder_diagonal_case_a()
    ladder_diagonal_case_b()

    # Q3(a)
    E, KX, KY = square_lattice_band()
    fillings = [0.5, 1.0, 1.5, 2.0]
    fermi_surface(E, KX, KY, fillings, prefix="Q3a")

    # Q3(b)
    for td in [-1.0, 1.0]:
        E, KX, KY = square_lattice_with_diagonal(td=td)
        fermi_surface(E, KX, KY, fillings, prefix=f"Q3b_td_{td}")