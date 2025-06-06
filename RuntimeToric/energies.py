from AuxEnergies import extract_energies, normalize_energies, extract_energies_chain
from AuxEnergies import plot_spectrum, plot_TOS, plot_spectrum_chain, plot_multiple_spectra, plot_CTES, analyze_Nmax
import numpy as np

lattice =  "square" 
Ns = 18
scan_type = "dual"    # "fix_hz"   "dual"
h_min = 0.00
h_max = 0.60
h_fix = 0.10
scan_interval = (h_min, h_max) if scan_type == "dual" else (h_min, h_max, h_fix)

df = extract_energies(lattice)
# print(df['hx'].unique())
# print(df['hz'].unique())
df = normalize_energies(df)

# Plot the spectrum 
plot_spectrum(df, Ns, scan_type, scan_interval, lattice, energy_cutoff=20, bareGaps=False)
# plot_spectrum(df, Ns, scan_type, scan_interval, lattice, energy_cutoff=20, bareGaps=True)

# Ns_list = [10, 16, 18, 20] 
# for Ns in Ns_list:
#     plot_spectrum(df, Ns, scan_type, scan_interval, lattice, energy_cutoff=30, bareGaps=False)
#     plot_spectrum(df, Ns, scan_type, scan_interval, lattice, energy_cutoff=30, bareGaps=True)

# Plot multiple spectra
# plot_multiple_spectra(df, scan_type, scan_interval, lattice, energy_cutoff=30, Rep_cutoff=False)

# # Plot the TOS
# U = 0.1
# plot_TOS(df, U, Nmax, Ns, t, lattice, energy_cutoff=5 , save_plot=False)

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
