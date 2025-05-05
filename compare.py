from EDroutine.AuxEnergies import extract_energies as extract_energies_1
from EDroutine.AuxEnergies import normalize_energies as normalize_energies_1
from EDroutine.AuxEnergies import duality_color

from RuntimeToric.AuxEnergies import extract_energies as extract_energies_2
from RuntimeToric.AuxEnergies import normalize_energies as normalize_energies_2
from RuntimeToric.AuxEnergies import Rep_sort_key

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os

def compare_spectrum(bareGaps=False):
    """
    Compare the spectrum of two different implementations of the same function.
    """ 

    lattice =  "square" 
    Ns = 18
    type = "dual"    # "fix_hz"   "dual"
    h_min = 0.00
    h_max = 0.60
    h_fix = 0.01
    scan_interval = (h_min, h_max) if type == "dual" else (h_min, h_max, h_fix)
    energy_cutoff = 20


    df1 = extract_energies_1(lattice)
    df1 = normalize_energies_1(df1)

    df2 = extract_energies_2(lattice)
    df2 = normalize_energies_2(df2)

    plt.figure(figsize=(18, 10))

    # ---- Plot routine for df1 ----
    if type in ['fix_hx', 'fix_hz']:
        if not isinstance(scan_interval, tuple) or len(scan_interval) != 3:
            raise ValueError("For 'fix_hx' and 'fix_hz', the 'interval' parameter must be a tuple with exactly three elements.")
    elif type == 'dual':
        if not isinstance(scan_interval, tuple) or len(scan_interval) != 2:
            raise ValueError("For 'dual', the 'interval' parameter must be a tuple with exactly two elements.")

    if type == 'fix_hx':
        df1_filtered = df1[(df1['Ns'] == Ns) & (df1['hx'] == scan_interval[2]) & (df1['hz'] >= scan_interval[0]) & (df1['hz'] <= scan_interval[1])].copy()
    elif type == 'fix_hz':
        df1_filtered = df1[(df1['Ns'] == Ns) & (df1['hz'] == scan_interval[2]) & (df1['hx'] >= scan_interval[0]) & (df1['hx'] <= scan_interval[1])].copy()
    elif type == 'dual':
        df1_filtered = df1[(df1['Ns'] == Ns) & (df1['hx'] == df1['hz']) & (df1['hx'] >= scan_interval[0]) & (df1['hx'] <= scan_interval[1])].copy()


    df1_filtered['scaled_energy'] = df1_filtered['Energy'] * np.sqrt(df1_filtered['Ns'])
    df1_filtered = df1_filtered[(df1_filtered['scaled_energy'] <= energy_cutoff)]
    size = 50
    linewidth = 1

    df1_filtered['color'] = df1_filtered['Duality'].apply(duality_color)

    for (h), group in df1_filtered.groupby(['hz'] if type == 'fix_hx' else ['hx']):
        if bareGaps:
            plt.scatter(group['hz'] if type == 'fix_hx' else group['hx'], group['Energy'], marker='o', edgecolors=group['color'], facecolors='none', linewidths=linewidth, s=size, label=None)
        else:
            plt.scatter(group['hz'] if type == 'fix_hx' else group['hx'], group['scaled_energy'], marker='o', edgecolors=group['color'], facecolors='none', linewidths=linewidth, s=size, label=None)
    

    # ---- Plot routine for df2 ----
    markers = ['o', 's', '^', 'v', '<', '>', 'p', '*', 'd', 'D', 'P', 'X', 'h', 'H', '8']
    legend_handles_reps = {}  # To store legend handles for Rep
    epsilon = 1e-10  # A small number close to machine precision


    if type == 'fix_hx':
        df2_filtered = df2[(df2['Ns'] == Ns) & (df2['hx'] == scan_interval[2]) & (df2['hz'] >= scan_interval[0]) & (df2['hz'] <= scan_interval[1])].copy()
    elif type == 'fix_hz':
        df2_filtered = df2[(df2['Ns'] == Ns) & (df2['hz'] == scan_interval[2]) & (df2['hx'] >= scan_interval[0]) & (df2['hx'] <= scan_interval[1])].copy()
    elif type == 'dual':
        df2_filtered = df2[(df2['Ns'] == Ns) & (df2['hx'] == df2['hz']) & (df2['hx'] >= scan_interval[0]) & (df2['hx'] <= scan_interval[1])].copy()

    unique_reps = np.array(sorted(df2_filtered['Rep'].unique(), key=Rep_sort_key))
    marker_dict = {rep: markers[i % len(markers)] for i, rep in enumerate(unique_reps)}

    df2_filtered['scaled_energy'] = df2_filtered['Energy'] * np.sqrt(df2_filtered['Ns'])
    df2_filtered = df2_filtered[(df2_filtered['scaled_energy'] <= energy_cutoff)]
    size = 150

    for (h, rep), group in df2_filtered.groupby(['hz', 'Rep'] if type == 'fix_hx' else ['hx', 'Rep']):
        marker = marker_dict[rep]
        if bareGaps:
            plt.scatter(group['hz'] if type == 'fix_hx' else group['hx'], group['Energy'], marker=marker, edgecolors='lime', facecolors='none', linewidths=1, s=size, label=None)
        else:
            plt.scatter(group['hz'] if type == 'fix_hx' else group['hx'], group['scaled_energy'], marker=marker, edgecolors='lime', facecolors='none', linewidths=1, s=size, label=None)
        # Create legend handles for reps if not already created
        if rep not in legend_handles_reps:
            handle_rep = mlines.Line2D([], [], marker=marker, linestyle='None', markeredgecolor='lime',
                                    markerfacecolor='none', markersize=10, label=rep)
            legend_handles_reps[rep] = handle_rep

    # ---- Position the legends on the plot ----
    legend_reps = plt.legend(handles=list(legend_handles_reps.values()), title='Space group irreps', loc='upper left', bbox_to_anchor=(1.005, 1))
    plt.gca().add_artist(legend_reps)  # This keeps the first legend visible when adding the second
    # Create handles for the Duality legend
    handle_plus1 = mlines.Line2D([], [], marker='s', linestyle='None', color='blue',
                                markersize=10, label='-1')
    handle_minus1 = mlines.Line2D([], [], marker='s', linestyle='None', color='red',
                                markersize=10, label='+1')
    handle_other = mlines.Line2D([], [], marker='s', linestyle='None', color='gray',
                                markersize=10, label='Other')
    # Add the legend to the plot
    plt.legend(handles=[handle_plus1, handle_minus1, handle_other], title='Duality', loc='lower left', bbox_to_anchor=(1.005, 0.1))
    plt.grid(True)
    plt.subplots_adjust(left=0.05, right=0.85, bottom=0.05, top=0.95)

    if type == 'fix_hx':
        plt.xlabel(r'$h_z$')
        plt.title(f'Energy Spectrum for Ns={Ns}, hx={scan_interval[2]}')
    elif type == 'fix_hz':
        plt.xlabel(r'$h_x$')
        plt.title(f'Energy Spectrum for Ns={Ns}, hz={scan_interval[2]}')
    elif type == 'dual':
        plt.xlabel(r'$h_x(h_z)$')
        plt.title(f'Energy Spectrum for Ns={Ns}, hx=hz')

    if bareGaps: plt.ylabel('Energy')
    else: plt.ylabel(r'$E\sqrt{N_s}$')
    
    # Decide to save the plot as a file or display it
    directory = 'plots'   # Directory path
    resolution = 400
    if not os.path.exists(directory): os.makedirs(directory)   # Create the directory if it doesn't exist
    
    file_suffix = f"_Ns{Ns}_hx{scan_interval[2]}_fix" if type == 'fix_hx' else \
                  f"_Ns{Ns}_hz{scan_interval[2]}_fix" if type == 'fix_hz' else \
                  f"_Ns{Ns}_hx(hz)_dual"
    file_prefix = f"{lattice}_spectrumCompare"
    plt.savefig(f'{directory}/{file_prefix}{file_suffix}.png', dpi=resolution)
    plt.show()



compare_spectrum()