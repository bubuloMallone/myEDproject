import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os
import re
import glob
from AuxTimeEvolution import extract_correlations, analyze_structure_factors, analyze_corr_length, plot_corr_length
from AuxTimeEvolution import analyze_structure_factors_U, analyze_structure_factors_Nmax, analyze_Dt_convergence, plot_Dt_convergence

lattice = "triangular"
operator= "bdagb"
Nmax = 2
t = -1
Ns = 12
Umax=100.0

# tot_time = 5.0
# Dt = 0.5

# df = extract_correlations(Nmax, Ns, t, Umax, lattice, operator, tot_time, Dt)
# df.to_csv('filename.csv', index=False)

if lattice == "kagome":
    if t == 1:   
        reference_momentum = "Gamma"
    if t == -1:  
        reference_momentum = "Gamma*"
elif lattice == "triangular":
    if t == 1:   
        reference_momentum = "Gamma"
    if t == -1:  
        reference_momentum = "K"
# reference_momentum = "C1.4*"

# Plot correlation lenght estimate as a function of time
Ns_vals = [9, 12]
tot_times= [0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0]     # 0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0
Dts = [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10.0]    # 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10.0
df = analyze_corr_length(lattice, Nmax, Ns_vals, t, Umax, tot_times, Dts, operator, reference_momentum)
plot_corr_length(df, lattice, Nmax, t, Umax, reference_momentum, logscale=True)

# Ns_vals = [12]
# tot_times= [1000.0]     # 0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0
# Dts = [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10.0]    #   0.001, 0.01, 0.1
# df = analyze_Dt_convergence(lattice, Nmax, Ns_vals, t, Umax, tot_times, Dts, operator, reference_momentum)
# plot_Dt_convergence(df, lattice, Nmax, t, Umax, reference_momentum, logscale=False)



# # Plot structure factors peaks as a function of U
# Ns_vals = [9,12]
# Nmax = 2
# t = -1
# analyze_structure_factors_U(lattice, Nmax, Ns_vals, t, U_vals, operator)



# # Plot structure factors on the BZ
# Ns = 12
# t = -1
# U_vals = [0.0, 2.0, 10.0, 100.0]
# for U in [0.0]:   # 2.0, 10.0, 100.0
#     analyze_structure_factors(lattice, Nmax, Ns, t, U, operator)




# # Plot structure factors peaks as a function of Nmax
# t = -1
# U = 0.5
# Ns_vals = [12]
# analyze_structure_factors_Nmax(lattice, Ns_vals, t, U, operator)