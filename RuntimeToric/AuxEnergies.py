import os
import re
import glob
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

def extract_energies(lattice):
    directory = f"/Users/pietro/myEDproject/RuntimeToric/OutToric"
    path_pattern = f"Energies.Toric_*.{lattice}.hx.*.hz.*.rep.*.dat"
    files = glob.glob(os.path.join(directory, path_pattern))
    files = sorted(files)
    all_data = []
    
    for file_path in files:
        # print(file_path)
        Ns_match = re.search(fr'Toric_(\d+)', file_path)  # Matches Ns after the specified lattice type
        hx_match = re.search(r'hx\.(\d+\.\d+)', file_path) 
        hz_match = re.search(r'hz\.(\d+\.\d+)', file_path)
        Rep_match = re.search(r'rep\.([^.]+.*?)(?=\.dat)', file_path)  # Matches string Rep
        
        if hx_match and Rep_match and hz_match and Ns_match:
            hx = float(hx_match.group(1))
            hz = float(hz_match.group(1))
            Rep = Rep_match.group(1)
            Ns = int(Ns_match.group(1))
            
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith('Eigenvalue'):
                        energy = float(line.split('=')[1])
                        all_data.append({'Ns': Ns, 'hx': hx, 'hz': hz, 'Rep': Rep, 'Energy': energy})
    
    df = pd.DataFrame(all_data)
    return df


def normalize_energies(df):
    """
    Normalizes the energy values in the DataFrame for each 'Ns', by subtracting the minimum energy
    for each group of "U", where minimum energies are calculated from a subset filtered by 'Nb'
    equal to 'Ns'.

    Parameters:
    df (pd.DataFrame): DataFrame containing 'Ns', 'U', 'Nb', 'Rep', and 'Energy' columns.

    Returns:
    pd.DataFrame: DataFrame with normalized energy values.
    """
    # Iterate over each unique 'Ns' in the DataFrame
    for Ns in df['Ns'].unique():
        for hz in df['hz'].unique():
            for hx in df['hx'].unique():
                # Filter the DataFrame to include only the rows with the specified 'hz' and 'hx'
                df_filtered = df[(df['Ns'] == Ns) & (df['hz'] == hz) & (df['hx'] == hx)].copy()
                # Calculate the minimum energy for each 'hz' and 'hx' in the filtered DataFrame
                min_energy = df_filtered['Energy'].min()
                # min_rep = df_filtered[df_filtered['Energy'] == min_energy]['Rep'].values[0]
                # print(f"Ns={Ns}, hz={hz}, hx={hx}, min_energy={min_energy}, rep={min_rep}")

                # Adjust all rows in the original DataFrame for this 'Ns', 'hz', and 'hx'
                df.loc[(df['Ns'] == Ns) & (df['hz'] == hz) & (df['hx'] == hx), 'Energy'] -= min_energy

    return df


def find_first_excited(group):
    # Filter out the ground state energy
    min_energy = group['Energy'].min()
    excited_energies = group[group['Energy'] > min_energy]['Energy']
    if not excited_energies.empty:
        return excited_energies.min()  # Return the minimum of the excited energies
    return pd.NA  # Return NA if no excited state is found


def Rep_sort_key(s):
    # Check if the first character is a letter and return a tuple that puts letters first
    return (s[0].isdigit(), s)


def plot_spectrum(df, Ns, type, interval, lattice, energy_cutoff=20, bareGaps=False):
    # Define marker styles for different 'Rep' values and colors for 'Nb'
    markers = ['o', 's', '^', 'v', '<', '>', 'p', '*', 'd', 'D', 'P', 'X', 'h', 'H', '8']
    # colors = ['r', 'g', 'b', 'c', 'magenta', 'gold', 'blueviolet', 'orange', 'saddlebrown']
    # colors = ['black', 'firebrick', 'g', 'goldenrod', 'b', 'gold', 'lime', 'r', 'dimgray']
    plt.figure(figsize=(18, 10))

    legend_handles_reps = {}  # To store legend handles for Rep
    # legend_handles_dNb = {}    # To store legend handles for Nb
    epsilon = 1e-10  # A small number close to machine precision

    if type in ['fix_hx', 'fix_hz']:
        if not isinstance(interval, tuple) or len(interval) != 3:
            raise ValueError("For 'fix_hx' and 'fix_hz', the 'interval' parameter must be a tuple with exactly three elements.")
    elif type == 'dual':
        if not isinstance(interval, tuple) or len(interval) != 2:
            raise ValueError("For 'dual', the 'interval' parameter must be a tuple with exactly two elements.")

    if type == 'fix_hx':
        df_filtered = df[(df['Ns'] == Ns) & (df['hx'] == interval[2]) & (df['hz'] >= interval[0]) & (df['hz'] <= interval[1])].copy()
        # df_filtered = df_filtered[((10 * df_filtered['hz']) % 5).abs() < epsilon]
    elif type == 'fix_hz':
        df_filtered = df[(df['Ns'] == Ns) & (df['hz'] == interval[2]) & (df['hx'] >= interval[0]) & (df['hx'] <= interval[1])].copy()
        # df_filtered = df_filtered[((10 * df_filtered['hx']) % 5).abs() < epsilon]
    elif type == 'dual':
        df_filtered = df[(df['Ns'] == Ns) & (df['hx'] == df['hz']) & (df['hx'] >= interval[0]) & (df['hx'] <= interval[1])].copy()
        # df_filtered = df_filtered[((10 * df_filtered['hx']) % 5).abs() < epsilon]

    unique_reps = np.array(sorted(df_filtered['Rep'].unique(), key=Rep_sort_key))
    print(f"Unique representations: {unique_reps}")
    marker_dict = {rep: markers[i % len(markers)] for i, rep in enumerate(unique_reps)}
    # color_dict = {rep: colors[i % len(colors)] for i, rep in enumerate(unique_reps)}

    df_filtered['scaled_energy'] = df_filtered['Energy'] * np.sqrt(df_filtered['Ns'])
    df_filtered = df_filtered[(df_filtered['scaled_energy'] <= energy_cutoff)]
    size = 150

    for (h, rep), group in df_filtered.groupby(['hz', 'Rep'] if type == 'fix_hx' else ['hx', 'Rep']):
        # rep = rep[0] if isinstance(rep, tuple) else rep  # Ensure rep is a string
        # print(h,rep)
        marker = marker_dict[rep]
        # color = color_dict[rep]
        if bareGaps:
            plt.scatter(group['hz'] if type == 'fix_hx' else group['hx'], group['Energy'], marker=marker, edgecolors='g', facecolors='none', linewidths=1.5, s=size, label=None)
        else:
            plt.scatter(group['hz'] if type == 'fix_hx' else group['hx'], group['scaled_energy'], marker=marker, edgecolors='b', facecolors='none', linewidths=1.5, s=size, label=None)
        # Create legend handles for reps if not already created
        if rep not in legend_handles_reps:
            handle_rep = mlines.Line2D([], [], marker=marker, linestyle='None', markeredgecolor='black',
                                    markerfacecolor='none', markersize=10, label=rep)
            legend_handles_reps[rep] = handle_rep

    if type == 'fix_hx':
        plt.xlabel(r'$h_z$')
        plt.title(f'Energy Spectrum for Ns={Ns}, hx={interval[2]}')
    elif type == 'fix_hz':
        plt.xlabel(r'$h_x$')
        plt.title(f'Energy Spectrum for Ns={Ns}, hz={interval[2]}')
    elif type == 'dual':
        plt.xlabel(r'$h_x(h_z)$')
        plt.title(f'Energy Spectrum for Ns={Ns}, hx=hz')

    if bareGaps: plt.ylabel('Energy')
    else: plt.ylabel(r'$E\sqrt{N_s}$')
    
    # Position the legends on the plot
    legend_reps = plt.legend(handles=list(legend_handles_reps.values()), title='Space group irreps', loc='upper left', bbox_to_anchor=(1.005, 1))
    plt.gca().add_artist(legend_reps)  # This keeps the first legend visible when adding the second
    # plt.legend(handles=list(legend_handles_dNb.values()), title='Particle Number', loc='upper left', bbox_to_anchor=(1.005, 0.4))
    plt.grid(True)
    plt.subplots_adjust(left=0.05, right=0.85, bottom=0.05, top=0.95)
    # Decide to save the plot as a file or display it
    directory = 'plots'   # Directory path
    resolution = 400
    if not os.path.exists(directory): os.makedirs(directory)   # Create the directory if it doesn't exist
    
    file_suffix = f"_Ns{Ns}_hx{interval[2]}_fix" if type == 'fix_hx' else \
                  f"_Ns{Ns}_hz{interval[2]}_fix" if type == 'fix_hz' else \
                  f"_Ns{Ns}_hx(hz)_dual"
    file_prefix = f"{lattice}_spectrumBareGaps" if bareGaps else f"{lattice}_spectrumNs"
    plt.savefig(f'{directory}/{file_prefix}{file_suffix}.png', dpi=resolution)
    plt.show()


def plot_TOS(df, selected_U, Nmax, Ns, t, lattice, energy_cutoff=20, save_plot=False):
    # Ensure selected_U is the correct type, e.g., if 'U' is stored as strings:
    selected_U = float(selected_U)  # Use float(selected_U) if 'U' is stored as floating points in the DataFrame

    # Define marker styles for different 'Rep' values and colors for 'Nb'
    # markers = ['o', 's', '^', 'v', '<', '>', '8', 'p', '*', 'd', 'D', 'H', 'h', 'P', 'X']
    # colors = ['red', 'blue', 'green', 'purple', 'darkorange', 'darkcyan', 'magenta', 'lime', 'lightpink', 'saddlebrown', 'gray', 'olive', 'deepskyblue', 'gold', 'black']
    markers = ['o', 's', '^', 'v', '<', '>', 'p', '*', 'd', 'D', 'P', 'X', 'h', 'H', '8']
    colors = ['black', 'firebrick', 'g', 'gold', 'b', 'yellow', 'lime', 'r', 'dimgray']
    legend_handles_dNb = {}  # To store legend handles for dNb
    legend_handles_reps = {}  # To store legend handles for Rep

    df_filtered = df[(df['Ns'] == Ns) & (df['U'] == selected_U)].copy()
    # Pre-calculate deltaNb and scaled energy for all filtered entries
    df_filtered['deltaNb'] = df_filtered['Nb'] - df_filtered['Ns']
    df_filtered['dNb2'] = df_filtered['deltaNb']**2
    df_filtered['scaled_energy'] = df_filtered['Energy'] * np.sqrt(df_filtered['Ns'])

    # Create dictionaries for markers and colors paired with each Rep
    unique_reps = np.array(sorted(df_filtered['Rep'].unique(), key=Rep_sort_key))
    unique_dNb = np.sort(df_filtered['deltaNb'].unique())
    marker_dict = {rep: markers[i % len(markers)] for i, rep in enumerate(unique_reps)}
    color_dict = {dNb: colors[i % len(colors)] for i, dNb in enumerate(unique_dNb)}
    # Apply additional filtering on scaled_energy if needed
    df_filtered = df_filtered[df_filtered['scaled_energy'] <= energy_cutoff]

    # Debugging output
    print(f"Filtered DataFrame for U={selected_U} has {len(df_filtered)} rows.")
    if df_filtered.empty:
        print("No data to plot.")
        return
    
    plt.figure(figsize=(10, 8))

    # Group by 'Nb' and 'Rep' and plot for the selected 'U'
    for (dNb, rep), group in df_filtered.groupby(['deltaNb', 'Rep']):
        marker = marker_dict[rep]
        color = color_dict[dNb]
        if dNb <= 0 : size=100 
        else: size=200
        plt.scatter(group['dNb2'], group['Energy'], marker=marker, edgecolors=color, facecolors='none', linewidths=1.5, s=size, label=None)

        if rep not in legend_handles_reps:
                handle_reps = mlines.Line2D([], [], marker=marker, linestyle='None', markeredgecolor='black',
                                        markerfacecolor='none', markersize=10, label=rep)
                legend_handles_reps[rep] = handle_reps
        if dNb not in legend_handles_dNb:
            handle_dnb = mlines.Line2D([], [], color=color, marker='s', linestyle='None', markeredgecolor=color, 
                                    markerfacecolor='none', markersize=10, markeredgewidth=2, label=f'ΔNb = {dNb}')
            legend_handles_dNb[dNb] = handle_dnb

    plt.xlabel(r'$(\Delta N_b)^2$')
    plt.ylabel('Energy')
    plt.title(f'TOS Analysis for Nmax={Nmax}, Ns={Ns}, t={t}, U={selected_U}')

    # # Adjust the x-axis tick labels to show Nb instead of Nb^2
    # squared_values = df_filtered['dNb2'].unique()
    # original_values = df_filtered['deltaNb'].unique()
    # plt.xticks(squared_values, original_values)  # Set x-ticks to original Nb values

    # Position the legends on the plot
    legend_reps = plt.legend(handles=list(legend_handles_reps.values()), title='Representations', loc='upper left', bbox_to_anchor=(1.005, 1))
    plt.gca().add_artist(legend_reps)
    plt.legend(handles=list(legend_handles_dNb.values()), title='Particle Number', loc='upper left', bbox_to_anchor=(1.005, 0.4))
    plt.grid(True)
    plt.subplots_adjust(left=0.05, right=0.8, bottom=0.1, top=0.95)
    
    # Decide to save the plot as a file or display it
    directory = 'plots'   # Directory path
    resolution = 400

    if not os.path.exists(directory): os.makedirs(directory)   # Create the directory if it doesn't exist

    plt.savefig(f'{directory}/{lattice}_TOS_Nmax{Nmax}_Ns{Ns}_t{t}_U{selected_U}_Rep_All.png', dpi = resolution)
    plt.show()


def plot_CTES(df, selected_U, Nmax, t, lattice, energy_cutoff=20, save_plot=False):
    # Ensure selected_U is the correct type, e.g., if 'U' is stored as strings:
    selected_U = float(selected_U)  # Use float(selected_U) if 'U' is stored as floating points in the DataFrame
    legend_handles_dNb = {}  # To store legend handles for Rep
    legend_handles_reps = {}    # To store legend handles for Nb

    if t == -1 : 
        # df_filtered = df[(df['U'] == selected_U) & (df['Rep'].str.contains(r'^(?:Gamma|K)\.'))].copy()   
        df_filtered = df[(df['U'] == selected_U) & ((df['Rep'] == 'Gamma.D6.A1') | (df['Rep'] == 'Gamma.C6.A') | (df['Rep'] == 'Gamma.D6.A2') |
                                                    (df['Rep'] == 'Gamma.D6.B1') | (df['Rep'] == 'Gamma.C6.B') | # (df['Rep'] == 'Gamma.D6.B2') |   # triang
                                                    # (df['Rep'] == 'K.D3.A1') | (df['Rep'] == 'K.C3.A') | (df['Rep'] == 'K.D3.A2') # triang
                                                    (df['Rep'] == 'Gamma.C6.E2a') | (df['Rep'] == 'Gamma.C6.E2b') | (df['Rep'] == 'Gamma.D6.E2')   # kagome
                                                    )].copy()
        #   (df['Rep'] == 'Gamma.C6.E1a') | (df['Rep'] == 'Gamma.C6.E1b') 
        #   (df['Rep'] == 'Gamma.C6.E2a') | (df['Rep'] == 'Gamma.C6.E2b')
    if t == 1 : 
        df_filtered = df[(df['U'] == selected_U) & (df['Rep'].str.contains(r'^(?:Gamma)\.'))].copy() 
        df_filtered = df[(df['U'] == selected_U) & ((df['Rep'] == 'Gamma.D6.A1') | (df['Rep'] == 'Gamma.D6.B1') |
                                                    (df['Rep'] == 'Gamma.D6.A2') | (df['Rep'] == 'Gamma.D6.B2') |
                                                    (df['Rep'] == 'Gamma.C6.A') | (df['Rep'] == 'Gamma.C6.B')
                                                    )].copy()  

    # Define marker styles for different 'Rep' values and colors for 'Nb'
    markers = ['o', 's', '^', 'v', '<', '>', 'p', '*', 'd', 'D', 'P', 'X', 'h', 'H', '8']
    # colors = ['r', 'g', 'b', 'c', 'magenta', 'gold', 'blueviolet', 'orange', 'saddlebrown']
    colors = ['black', 'firebrick', 'g', 'goldenrod', 'b', 'gold', 'lime', 'r', 'dimgray']
    # Debugging output
    print(f"Filtered DataFrame for U={selected_U} has {len(df_filtered)} rows.")
    if df_filtered.empty:
        print("No data to plot.")
        return
    
    # Pre-calculate deltaNb and scaled energy for all filtered entries
    df_filtered['deltaNb'] = df_filtered['Nb'] - df_filtered['Ns']
    df_filtered['scaled_energy'] = df_filtered['Energy'] * np.sqrt(df_filtered['Ns'])
    df_filtered['1/Ns'] = 1.0 / df_filtered['Ns']   

    # Create dictionaries for markers and colors paired with each Rep
    unique_reps = np.array(sorted(df_filtered['Rep'].unique(), key=Rep_sort_key))
    unique_dNb = np.sort(df_filtered['deltaNb'].unique())
    marker_dict = {rep: markers[i % len(markers)] for i, rep in enumerate(unique_reps)}
    color_dict = {dNb: colors[i % len(colors)] for i, dNb in enumerate(unique_dNb)}
    # Apply additional filtering on scaled_energy if needed
    df_filtered = df_filtered[df_filtered['scaled_energy'] <= energy_cutoff]

    plt.figure(figsize=(10, 8))
    # Group by 'Nb' and 'Rep' and plot for the selected 'U'
    for (Ns, rep, dNb), group in df_filtered.groupby(['Ns', 'Rep', 'deltaNb']):
        marker = marker_dict[rep]
        color = color_dict[dNb]
        if dNb <= 0 : size=100 
        else: size=200
        if t==-1:
            if (rep == 'Gamma.D6.B1'): group['1/Ns'] = group['1/Ns'] - 0.002
            elif (rep == 'Gamma.C6.B'): group['1/Ns'] = group['1/Ns'] - 0.002
            if (rep == 'K.D3.A1'): group['1/Ns'] = group['1/Ns'] + 0.002
            elif (rep == 'K.C3.A'): group['1/Ns'] = group['1/Ns'] +0.002
            if (rep == 'Gamma.D6.A2'): group['1/Ns'] = group['1/Ns'] - 0.002
            if (rep == 'Gamma.D6.E2'): group['1/Ns'] = group['1/Ns'] + 0.002
            elif (rep == 'Gamma.C6.E2a'): group['1/Ns'] = group['1/Ns'] + 0.002
            elif (rep == 'Gamma.C6.E2b'): group['1/Ns'] = group['1/Ns'] - 0.002
        if t==1:
            if (rep == 'Gamma.D6.E2'): group['1/Ns'] = group['1/Ns'] + 0.002
            elif (rep == 'Gamma.C6.E2a'): group['1/Ns'] = group['1/Ns'] + 0.002
            elif (rep == 'Gamma.C6.E2b'): group['1/Ns'] = group['1/Ns'] + 0.002

        plt.scatter(group['1/Ns'], group['scaled_energy'], marker=marker, edgecolors=color, facecolors='none', linewidths=2, s=size, label=None)
        # Create legend handles for reps if not already created
        if rep not in legend_handles_reps:
            handle_reps = mlines.Line2D([], [], marker=marker, linestyle='None', markeredgecolor='black',
                                    markerfacecolor='none', markersize=10, label=rep)
            legend_handles_reps[rep] = handle_reps
        # Create legend handles for Nb if not already created
        if dNb not in legend_handles_dNb:
            handle_dnb = mlines.Line2D([], [], color=color, marker='s', linestyle='None', markeredgecolor=color, 
                                    markerfacecolor='none', markersize=10, markeredgewidth=2, label=f'ΔNb = {dNb}')
            legend_handles_dNb[dNb] = handle_dnb

    plt.xlabel('1/Ns')
    plt.ylabel(r'$E\sqrt{N_s}$')
    plt.title(f'CTES for Nmax={Nmax}, t={t}, U={selected_U}')

    # Position the legends on the plot
    legend_reps = plt.legend(handles=list(legend_handles_reps.values()), title='Representations', loc='upper left', bbox_to_anchor=(1.005, 1))
    plt.gca().add_artist(legend_reps)
    plt.legend(handles=list(legend_handles_dNb.values()), title='Particle Number', loc='upper left', bbox_to_anchor=(1.005, 0.4))
    plt.grid(True)
    plt.subplots_adjust(left=0.1, right=0.8, bottom=0.1, top=0.95)
    
    # Decide to save the plot as a file or display it
    directory = 'plots'   # Directory path
    resolution = 400
    if not os.path.exists(directory): os.makedirs(directory)   # Create the directory if it doesn't exist
    plt.savefig(f'{directory}/{lattice}_CTES_Nmax{Nmax}_Ns{Ns}_t{t}_U{selected_U}.png', dpi = resolution)
    plt.show()


def plot_multiple_spectra(df, type, interval, lattice, energy_cutoff=40, Rep_cutoff=False):
    # Define marker styles for different 'Rep' values and colors for 'Ns'
    markers = ['o', 's', '^', 'v', '<', '>', 'p', '*', 'd', 'D', 'P', 'X', 'h', 'H', '8']
    colors = ['red', 'green', 'blue', 'cyan', 'magenta', 'gold', 'blueviolet', 'orange', 'saddlebrown']

    if type in ['fix_hx', 'fix_hz']:
        if not isinstance(interval, tuple) or len(interval) != 3:
            raise ValueError("For 'fix_hx' and 'fix_hz', the 'interval' parameter must be a tuple with exactly three elements.")
    elif type == 'dual':
        if not isinstance(interval, tuple) or len(interval) != 2:
            raise ValueError("For 'dual', the 'interval' parameter must be a tuple with exactly two elements.")

    if type == 'fix_hx':
        df_filtered = df[(df['hx'] == interval[2]) & (df['hz'] >= interval[0]) & (df['hz'] <= interval[1])].copy()
        # df_filtered = df_filtered[((10 * df_filtered['hz']) % 5).abs() < epsilon]
    elif type == 'fix_hz':
        df_filtered = df[(df['hz'] == interval[2]) & (df['hx'] >= interval[0]) & (df['hx'] <= interval[1])].copy()
    elif type == 'dual':
        df_filtered = df[(df['hx'] == df['hz']) & (df['hx'] >= interval[0]) & (df['hx'] <= interval[1])].copy()

    # Apply Rep cutoff filter if True
    if Rep_cutoff: df_filtered = df_filtered[df_filtered['Rep'].str.contains(r'^(?:Gamma)\.')]
    # Pre-calculate deltaNb and scaled energy for all filtered entries
    df_filtered['scaled_energy'] = df_filtered['Energy'] * np.sqrt(df_filtered['Ns'])
    unique_Ns = np.sort(df_filtered['Ns'].unique())
    unique_reps = np.array(sorted(df_filtered['Rep'].unique(), key=Rep_sort_key))
    # Create dictionaries for markers and colors
    marker_dict = {rep: markers[i % len(markers)] for i, rep in enumerate(unique_reps)}
    color_dict = {Ns: colors[i % len(colors)] for i, Ns in enumerate(unique_Ns)}
    # Apply additional filtering on scaled_energy if needed
    df_filtered = df_filtered[df_filtered['scaled_energy'] <= energy_cutoff]
    legend_handles_reps = {}  # To store legend handles for Rep
    legend_handles_Ns = {}

    plt.figure(figsize=(16, 8))
    for (h, rep, Ns), group in df_filtered.groupby(['hz', 'Rep', 'Ns'] if type == 'fix_hx' else ['hx', 'Rep', 'Ns']):
        marker = marker_dict[rep]
        color = color_dict[Ns]
        scaled_energy = group['Energy'] * np.sqrt(Ns)
        plt.scatter(group['hz'] if type == 'fix_hx' else group['hx'], scaled_energy, marker=marker, edgecolors=color, facecolors='none', linewidths=2, s=100, label=None)
        # Create legend handles for reps if not already created
        if rep not in legend_handles_reps:
            handle_rep = mlines.Line2D([], [], marker=marker, linestyle='None', markeredgecolor='black',
                                        markerfacecolor='none', markersize=10, label=rep)
            legend_handles_reps[rep] = handle_rep
        if Ns not in legend_handles_Ns:
            handle_Ns = mlines.Line2D([], [], marker='s', linestyle='None', color = color,
                                        markersize=10, label=f'Ns = {Ns}')
            legend_handles_Ns[Ns] = handle_Ns

    plt.grid(True)

    # Adjust layout
    plt.subplots_adjust(left=0.05, right=0.85, bottom=0.05, top=0.95)
    legend_reps = plt.legend(handles=list(legend_handles_reps.values()), title='Space group irreps', loc='upper left', bbox_to_anchor=(1.005, 1))
    plt.gca().add_artist(legend_reps)  # This keeps the first legend visible when adding the second
    plt.legend(handles=list(legend_handles_Ns.values()), title='Particle Number', loc='upper left', bbox_to_anchor=(1.005, 0.4))

    if type == 'fix_hx':
        plt.xlabel(r'$h_z$')
        plt.title(f'Gap scaling with Ns, hx={interval[2]}')
    elif type == 'fix_hz':
        plt.xlabel(r'$h_x$')
        plt.title(f'Gap scaling with Ns, hz={interval[2]}')
    elif type == 'dual':
        plt.xlabel(r'$h_x(h_z)$')
        plt.title(f'Gap scaling with system size, hx=hz')
    plt.ylabel(r'$E\sqrt{N_s}$')

    # Decide to save the plot as a file or display it
    directory = 'plots'   # Directory path
    resolution = 400
    if not os.path.exists(directory): os.makedirs(directory)   # Create the directory if it doesn't exist
    
    file_suffix = f"_hx{interval[2]}_fix" if type == 'fix_hx' else \
                  f"_hz{interval[2]}_fix" if type == 'fix_hz' else \
                  f"_hx(hz)_dual"
    file_prefix = f"{lattice}_spectrumMultiScale"
    plt.savefig(f'{directory}/{file_prefix}{file_suffix}.png', dpi=resolution)
    plt.show()


def analyze_Nmax(Nmax_vals, t, lattice, U, Ns, energy_cutoff=20):

    results_df = pd.DataFrame()
    for Nmax in Nmax_vals:
        df = extract_energies(Nmax, t, lattice)
        # df = chemPot_correction(df)
        df = normalize_energies(df)
        df['Nmax'] = Nmax
        results_df = pd.concat([results_df, df], ignore_index=True)

    print("Results DataFrame constructed for Nmax = ", Nmax_vals, "\n")

    # filter DataFrame for U values 
    # results_df = results_df[(results_df['U'] <= Umin) & (results_df['U'] <= Umax)]
    # print("Result DataFrame filtered for [Umin, Umax] = ", Umin, ", ", Umax, "\n")
    # Ensure selected_U is the correct type, e.g., if 'U' is stored as strings:
    U = float(U)  # Use float(selected_U) if 'U' is stored as floating points in the DataFrame
    results_df = results_df[(results_df['U'] == U)]
    print("Result DataFrame filtered for U = ", U, "\n")
    results_df = results_df[(results_df['Ns'] == Ns)]

    # Pre-calculate deltaNb and scaled energy for all filtered entries
    results_df['deltaNb'] = results_df['Nb'] - results_df['Ns']
    # Define marker styles for different 'Rep' values and colors for 'Nb'
    markers = ['o', 's', '^', 'v', '<', '>', 'p', '*', 'd', 'D', 'P', 'X', 'h', 'H', '8']
    colors = ['black', 'firebrick', 'g', 'goldenrod', 'b', 'gold', 'lime', 'r', 'dimgray']
    legend_handles_dNb = {}  # To store legend handles for dNb
    legend_handles_reps = {}  # To store legend handles for Rep

    # Create dictionaries for markers and colors paired with each Rep
    unique_reps = np.array(sorted(results_df['Rep'].unique(), key=Rep_sort_key))
    print(unique_reps)
    unique_dNb = np.sort(results_df['deltaNb'].unique())
    marker_dict = {rep: markers[i % len(markers)] for i, rep in enumerate(unique_reps)}
    color_dict = {dNb: colors[i % len(colors)] for i, dNb in enumerate(unique_dNb)}

    results_df['scaled_energy'] = results_df['Energy'] * np.sqrt(results_df['Ns'])
    # Apply additional filtering on scaled_energy if needed
    results_df = results_df[results_df['scaled_energy'] <= energy_cutoff]

    # Debugging output
    print(f"Filtered DataFrame for U={U} has {len(results_df)} rows.")
    if results_df.empty:
        print("No data to plot.")
        return
    
    plt.figure(figsize=(10, 8))
    # Group by 'Nb' and 'Rep' and plot for the selected 'U'
    for (dNb, rep), group in results_df.groupby(['deltaNb', 'Rep']):
        marker = marker_dict[rep]
        color = color_dict[dNb]
        if dNb <= 0 : size=100 
        else: size=200
        plt.scatter(group['Nmax'], group['Energy'], marker=marker, edgecolors=color, facecolors='none', linewidths=2, s=size, label=None)

        # Create legend handles for reps if not already created
        if rep not in legend_handles_reps:
            handle_rep = mlines.Line2D([], [], marker=marker, linestyle='None', markeredgecolor='black',
                                       markerfacecolor='none', markersize=10, label=rep)
            legend_handles_reps[rep] = handle_rep

        # Create legend handles for Nb if not already created
        if dNb not in legend_handles_dNb:
            handle_nb = mlines.Line2D([], [], color=color, marker='s', linestyle='None',
                                      markersize=10, markeredgewidth=2, label=f'ΔNb = {dNb}')
            legend_handles_dNb[dNb] = handle_nb

    plt.xlabel(r'$N_{max}$')
    plt.ylabel('Energy')
    plt.title(f'Spectrum vs Nmax, Ns={Ns}, t={t}, U={U}')
    
    # Position the legends on the plot
    legend_reps = plt.legend(handles=list(legend_handles_reps.values()), title='Representation', loc='upper left', bbox_to_anchor=(1.005, 1))
    plt.gca().add_artist(legend_reps)  # This keeps the first legend visible when adding the second
    plt.legend(handles=list(legend_handles_dNb.values()), title='Particle Number', loc='upper left', bbox_to_anchor=(1.005, 0.4))

    plt.grid(True)
    plt.subplots_adjust(left=0.07, right=0.8, bottom=0.07, top=0.95)
    
    # Decide to save the plot as a file or display it
    directory = 'plots'   # Directory path
    resolution = 400

    if not os.path.exists(directory): os.makedirs(directory)   # Create the directory if it doesn't exist
    plt.savefig(f'{directory}/{lattice}_SpectrumVsNmax_U{U}_Ns{Ns}_t{t}.png', dpi = resolution)
    plt.show()

    return









def plot_spectrum_chain(df, Nmax, Nn, t):
    
    # Define marker styles for different 'Rep' values to distinguish them in the plot
    markers = ['o', 's', '^', 'v', '<', '>', 'p', '*', '+', 'x', 'D']
    unique_reps = df['Rep'].unique()
    marker_dict = {rep: markers[i % len(markers)] for i, rep in enumerate(unique_reps)}

    plt.figure(figsize=(10, 6))
    
    # Track which representations have been plotted for the legend
    plotted_reps = set()

    # Group by 'U' and 'Rep' and plot for the selected 'Nb'
    for (U, rep), group in df.groupby(['U', 'Rep']):
        if rep not in plotted_reps:
            plt.scatter([U] * len(group), group['Energy'], label=rep, marker=marker_dict[rep], facecolors='none', edgecolors='black', s=100)
            plotted_reps.add(rep)
        else:
            plt.scatter([U] * len(group), group['Energy'], marker=marker_dict[rep], facecolors='none', edgecolors='black', s=100)
    
    plt.xlabel('U/t')
    plt.ylabel('Energy')
    title = f'Energy Spectrum for Nmax={Nmax}, Nn={Nn}, t={t}'
    plt.title(title)
    plt.legend(title='Representation', loc='upper left')  # bbox_to_anchor=(1.01, 1)
    plt.grid(True)
    plt.show()

def extract_energies_chain(Nmax, Nn, t):
    directory = "/Users/pietro/Runtime/Out"
    path_pattern = f"Energies.Bosons_{Nmax}.chain.{Nn}.t.{t}.U.*.N.*.k.*.dat"
    files = glob.glob(os.path.join(directory, path_pattern))
    # print(files)
    
    # Prepare a list to hold all data rows
    all_data = []
    
    for file_path in files:
        U_match = re.search(r'U\.(\d+\.\d+|\d+)', file_path)
        Rep_match = re.search(r'k\.([^.]+.*?)(?=\.dat)', file_path)
        if Rep_match:
            U = U_match.group(1)
            Rep = Rep_match.group(1)
            
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith('Eigenvalue'):
                        energy = float(line.split('=')[1])
                        # Append each found energy along with its Nb and r values to the list
                        all_data.append({'U': U, 'Rep': Rep, 'Energy': energy})
    
    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(all_data)
    
    return df