from AuxEnergies import extract_energies, normalize_energies, chemPot_correction, extract_energies_chain
from AuxEnergies import plot_spectrum, plot_TOS, plot_spectrum_chain, plot_multiple_spectra, plot_CTES, analyze_Nmax
import numpy as np

lattice =  "kagome" # triangular  # kagome
Nmax = 2   # 7
t = -1
Ns = 12

df = extract_energies(Nmax, t, lattice)
# print("U: ", np.sort(df['U'].unique()))
# print("Nb: ", np.sort(df['Nb'].unique()))
# print("Ns: ", np.sort(df['Ns'].unique()))
# print("Rep: ", np.sort(df['Rep'].unique()))
df = chemPot_correction(df)
df = normalize_energies(df)

if lattice == "kagome":
    if t == 1:   Umin, Umax = (10,20)  # (0,20)
    if t == -1:  Umin, Umax = (0,3)
elif lattice == "triangular":
    if t == 1:   Umin, Umax = (0,3)  #  (25.5,27.5)   (0,30)
    if t == -1:  Umin, Umax = (0,5)   #  (4,8)   (0,10)

# # Plot the spectrum 
# plot_spectrum(df, Nmax, Ns, t, lattice, Umin, Umax, energy_cutoff=50, Rep_cutoff=False, save_plot=True, bareGaps=False)
# plot_spectrum(df, Nmax, Ns, t, lattice, Umin, Umax, energy_cutoff=5, Rep_cutoff=False, save_plot=False, bareGaps=True)

# Ns_list = [9, 12, 21] 
# for Ns in Ns_list:
#     plot_spectrum(df, Nmax, Ns, t, lattice, Umin, Umax, energy_cutoff=50, Rep_cutoff=False, save_plot=True, bareGaps=False)
#     plot_spectrum(df, Nmax, Ns, t, lattice, Umin, Umax, energy_cutoff=50, Rep_cutoff=False, save_plot=True, bareGaps=True)

# # Plot multiple spectra
# plot_multiple_spectra(df, Nmax, t, lattice, Umin, Umax, energy_cutoff=10, Rep_cutoff=True, save_plot=False)

# Plot the TOS
U = 0.1
plot_TOS(df, U, Nmax, Ns, t, lattice, energy_cutoff=5 , save_plot=False)

# # Plot CTES
# Uc = 1.3
# plot_CTES(df, Uc, Nmax, t, lattice, energy_cutoff=7, save_plot=False)



# # Plot spectrum vs Nmax
# U_vals = np.arange(0.0, 1.5, 0.5)
# U_vals = [6.5]
# Nmax_vals = np.arange(2, 10)
# Ns=12
# for U in U_vals:
#     analyze_Nmax(Nmax_vals, t, lattice, U, Ns, energy_cutoff=20)
