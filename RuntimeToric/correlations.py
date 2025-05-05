import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os
import re
import glob
from AuxCorrelations import analyze_structure_factors, analyze_corr_length, plot_corr_length
from AuxCorrelations import analyze_structure_factors_U, analyze_structure_factors_Nmax

lattice = "triangular"
operator= "bdagb"

# Plot correlation lenght estimate as a function of U
Ns_vals = [12]
Nmax = 2
t = 1
if lattice == "kagome":
    if t == 1:
        reference_momentum = "Gamma"
        U_vals = np.arange(0.0, 20.5, 0.5)
        U_vals_add = np.concatenate((np.arange(21,25,1), np.arange(25, 85, 5)))
        U_vals = np.concatenate((U_vals,U_vals_add))
    if t == -1:
        reference_momentum = "Gamma*"  # K* Gamma* C1.4*
        U_vals = np.arange(0.00, 3.10, 0.10)
        U_vals_add = np.concatenate((np.arange(3.25, 10.25, 0.25), np.arange(11.00, 15.00, 1), np.arange(15, 42.5, 2.5)))   # , np.arange(15, 42.5, 2.5)
        U_vals = np.concatenate((U_vals,U_vals_add))
        # print(U_vals)
elif lattice == "triangular":
    if t == 1: 
        reference_momentum = "Gamma"
        U_vals = np.arange(0.0, 50.1, 0.1)  #  40.5
        # U_vals_add = np.arange(45, 105, 5)   #  (45, 105, 5)
        # U_vals = np.concatenate((U_vals,U_vals_add))
    if t == -1:
        reference_momentum = "K"
        U_vals = np.arange(0.0, 21, 1)
        U_vals_add = np.arange(22.5, 42.5, 2.5)
        U_vals = np.concatenate((np.concatenate((U_vals,U_vals_add)), np.arange(45.00, 55.00, 5)))

df_k = analyze_corr_length(lattice, Nmax, Ns_vals, t, U_vals, operator, reference_momentum)
plot_corr_length(df_k, lattice, Nmax, Ns_vals, t, reference_momentum)


# # Plot structure factors peaks as a function of U
# Ns_vals = [9,12]
# Nmax = 2
# t = -1
# analyze_structure_factors_U(lattice, Nmax, Ns_vals, t, U_vals, operator)




# # Plot structure factors on the BZ
# Ns = 12
# t = 1
# U_vals = [0.0, 2.0, 10.0, 100.0]
# for U in [100.00]:   
#     analyze_structure_factors(lattice, Nmax, Ns, t, U, operator)




# # Plot structure factors peaks as a function of Nmax
# t = 1
# U = 0.5
# Ns_vals = [12]
# analyze_structure_factors_Nmax(lattice, Ns_vals, t, U, operator)