import os
import re
import glob
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

def extract_energies(Nmax, t, lattice):
    directory = f"/Users/pietro/Runtime/Out{lattice}"
    path_pattern = f"Energies.Bosons_{Nmax}.{lattice}.*.t.{t}.U.*.Nb.*.rep.*.dat"
    files = glob.glob(os.path.join(directory, path_pattern))
    
    all_data = []
    
    for file_path in files:
        Ns_match = re.search(fr'{lattice}\.(\d+)', file_path)  # Matches Ns after the specified lattice type
        U_match = re.search(r'U\.(\d+\.\d+)', file_path)  # Matches floating-point U
        Nb_match = re.search(r'Nb\.(\d+)', file_path)     # Matches integer Nb
        Rep_match = re.search(r'rep\.([^.]+.*?)(?=\.dat)', file_path)  # Matches string Rep
        
        if Nb_match and Rep_match and U_match and Ns_match:
            U = float(U_match.group(1))
            Nb = int(Nb_match.group(1))
            Rep = Rep_match.group(1)
            Ns = int(Ns_match.group(1))
            
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith('Eigenvalue'):
                        energy = float(line.split('=')[1])
                        all_data.append({'Ns': Ns, 'U': U, 'Nb': Nb, 'Rep': Rep, 'Energy': energy})
    
    df = pd.DataFrame(all_data)
    return df


def chemPot_correction(df):
    """
    Adjusts the energy values in the DataFrame by subtracting a calculated chemical potential for each Ns,
    directly modifying the original DataFrame.
    
    Parameters:
    df (pd.DataFrame): DataFrame containing 'Ns', 'U', 'Nb', 'Rep', and 'Energy' columns.
    """
    for Ns in df['Ns'].unique():
        # Select the subset for the specific Ns
        df_subset = df[df['Ns'] == Ns]

        # Compute minimum energies for each Nb and U combination within the subset
        E0 = df_subset.groupby(['Nb', 'U'])['Energy'].min().unstack()

        # Calculate chemical potentials mu1 and mu2 for each U
        mu1 = E0.loc[Ns-1] - E0.loc[Ns] if (Ns-1) in E0.index else pd.Series(0, index=E0.columns)
        mu2 = E0.loc[Ns] - E0.loc[Ns+1] if (Ns+1) in E0.index else pd.Series(0, index=E0.columns)
        # print(mu1)
        # print(mu2)
        # Calculate the mean chemical potential for each U
        mu = 0.5 * (mu1 + mu2)
        # print(mu)
        # Calculate adjusted energies and update directly in df
        adjust_energy = df_subset.apply(lambda row: row['Energy'] + mu.get(row['U'], 0) * row['Nb'], axis=1)
        df.loc[adjust_energy.index, 'Energy'] = adjust_energy

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
        # Filter the DataFrame to include only the rows with the specified 'Nb' equal to 'Ns'
        df_filtered = df[(df['Ns'] == Ns) & (df['Nb'] == Ns)]
        
        # Calculate the minimum energy for each 'U' in the filtered DataFrame
        min_energies = df_filtered.groupby('U')['Energy'].min()
        
        # Map these minimum energies back onto the original DataFrame based on 'U' to normalize all energies within this 'Ns' subset
        # Adjust all rows in the original DataFrame for this 'Ns', not just those with 'Nb' = Ns
        for U, min_energy in min_energies.items():
            df.loc[(df['Ns'] == Ns) & (df['U'] == U), 'Energy'] -= min_energy

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


def plot_spectrum(df, Nmax, Ns, t, lattice, Umin, Umax, energy_cutoff=20, Rep_cutoff=False, save_plot=False, bareGaps=False):
    # Define marker styles for different 'Rep' values and colors for 'Nb'
    markers = ['o', 's', '^', 'v', '<', '>', 'p', '*', 'd', 'D', 'P', 'X', 'h', 'H', '8']
    # colors = ['r', 'g', 'b', 'c', 'magenta', 'gold', 'blueviolet', 'orange', 'saddlebrown']
    colors = ['black', 'firebrick', 'g', 'goldenrod', 'b', 'gold', 'lime', 'r', 'dimgray']

    plt.figure(figsize=(18, 10))

    legend_handles_reps = {}  # To store legend handles for Rep
    legend_handles_dNb = {}    # To store legend handles for Nb

    # Filter the DataFrame to include only energy values less than or equal to the cutoff
    df_filtered = df[(df['Ns'] == Ns) & (df['U'] >= Umin) & (df['U'] <= Umax)].copy()   #  & (df['U'] % 0.1 == 0)
    epsilon = 1e-10  # A small number close to machine precision
    df_filtered = df_filtered[((10 * df_filtered['U']) % 5).abs() < epsilon]

    df_filtered['deltaNb'] = df_filtered['Nb'] - df_filtered['Ns']
    unique_reps = np.array(sorted(df_filtered['Rep'].unique(), key=Rep_sort_key))
    # print(unique_reps)
    unique_dNb = np.sort(df_filtered['deltaNb'].unique())
    # Create dictionaries for markers and colors
    marker_dict = {rep: markers[i % len(markers)] for i, rep in enumerate(unique_reps)}
    color_dict = {dNb: colors[i % len(colors)] for i, dNb in enumerate(unique_dNb)}

    # Apply Rep cutoff filter if True
    if Rep_cutoff:
        # df_filtered = df_filtered[(df_filtered['Rep'].str.contains(r'^(?:Gamma|K)\.'))] 
        df_filtered = df_filtered[((df['Rep'] == 'Gamma.D6.A1') | (df['Rep'] == 'Gamma.D6.B1') | (df['Rep'] == 'K.D3.A1'))]
    
    df_filtered['scaled_energy'] = df_filtered['Energy'] * np.sqrt(df_filtered['Ns'])
    df_filtered = df_filtered[(df_filtered['scaled_energy'] <= energy_cutoff)]

    # Plot each combination of 'U', 'Rep', and 'Nb' for the filtered DataFrame
    for (U, rep, dNb), group in df_filtered.groupby(['U', 'Rep', 'deltaNb']):
        marker = marker_dict[rep]
        color = color_dict[dNb]
        if dNb <= 0 : size=100 
        else: size=200
        # Plot data
        if bareGaps:
            plt.scatter(group['U'], group['Energy'], marker=marker, edgecolors=color, facecolors='none', linewidths=2, s=size, label=None)
        else:
            plt.scatter(group['U'], group['scaled_energy'], marker=marker, edgecolors=color, facecolors='none', linewidths=2, s=size, label=None)

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

    plt.xlabel('U/t')
    if bareGaps: plt.ylabel('Energy')
    else: plt.ylabel(r'$E\sqrt{N_s}$')
    plt.title(f'Energy Spectrum for Nmax={Nmax}, Ns={Ns}, t={t}')

    # Position the legends on the plot
    legend_reps = plt.legend(handles=list(legend_handles_reps.values()), title='Representation', loc='upper left', bbox_to_anchor=(1.005, 1))
    plt.gca().add_artist(legend_reps)  # This keeps the first legend visible when adding the second
    plt.legend(handles=list(legend_handles_dNb.values()), title='Particle Number', loc='upper left', bbox_to_anchor=(1.005, 0.4))
    
    plt.grid(True)
    plt.subplots_adjust(left=0.05, right=0.85, bottom=0.05, top=0.95)
    
    # Decide to save the plot as a file or display it
    directory = 'plots'   # Directory path
    resolution = 400

    if not os.path.exists(directory): os.makedirs(directory)   # Create the directory if it doesn't exist

    if save_plot:
        if bareGaps:
            plt.savefig(f'{directory}/{lattice}_spectrumBareGaps_Nmax{Nmax}_Ns{Ns}_t{t}_{"_Rep_GK" if Rep_cutoff else "_Rep_All"}.png', dpi = resolution)
        else: 
            plt.savefig(f'{directory}/{lattice}_spectrumNs_Nmax{Nmax}_Ns{Ns}_t{t}_{"_Rep_GK" if Rep_cutoff else "_Rep_All"}.png', dpi = resolution)
    else:
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


def plot_multiple_spectra(df, Nmax, t, lattice, Umin, Umax, energy_cutoff=40, Rep_cutoff=False, save_plot=False):
    # Define marker styles for different 'Rep' values and colors for 'Ns'
    markers = ['o', 's', '^', 'v', '<', '>', 'p', '*', 'd', 'D', 'P', 'X', 'h', 'H', '8']
    colors = ['red', 'green', 'blue', 'cyan', 'magenta', 'gold', 'blueviolet', 'orange', 'saddlebrown']

    # Filter the DataFrame to include only energy values less than or equal to the cutoff
    df_filtered = df[(df['U'] >= Umin) & (df['U'] <= Umax)].copy()  # (df['U'] % 0.5 == 0) &
    # epsilon = 1e-10  # A small number close to machine precision
    # df_filtered = df[(df['U'] % 0.5).abs() < epsilon]

    # Apply Rep cutoff filter if True
    if Rep_cutoff: df_filtered = df_filtered[(df['Rep'] == 'Gamma.D6.A1') | (df['Rep'] == 'Gamma.C6.A') | (df['Rep'] == 'Gamma.D6.A2') | 
                                             (df['Rep'] == 'Gamma.D6.B1') | (df['Rep'] == 'Gamma.C6.B') | # (df['Rep'] == 'Gamma.D6.B2') |   # triang 
                                            # (df['Rep'] == 'K.D3.A1') | (df['Rep'] == 'K.C3.A') # | (df['Rep'] == 'K.D3.A2')  # triang
                                             (df['Rep'] == 'Gamma.C6.E2a') | (df['Rep'] == 'Gamma.C6.E2b') | (df['Rep'] == 'Gamma.D6.E2')   # kagome
                                                    ]   # df_filtered['Rep'].str.contains(r'^(?:Gamma|K)\.')

    # Pre-calculate deltaNb and scaled energy for all filtered entries
    df_filtered['deltaNb'] = df_filtered['Nb'] - df_filtered['Ns']
    df_filtered['scaled_energy'] = df_filtered['Energy'] * np.sqrt(df_filtered['Ns'])

    unique_Ns = np.sort(df_filtered['Ns'].unique())
    unique_reps = np.sort(df_filtered['Rep'].unique())
    # Create dictionaries for markers and colors
    marker_dict = {rep: markers[i % len(markers)] for i, rep in enumerate(unique_reps)}
    color_dict = {Ns: colors[i % len(colors)] for i, Ns in enumerate(unique_Ns)}

    # Apply additional filtering on scaled_energy if needed
    df_filtered = df_filtered[df_filtered['scaled_energy'] <= energy_cutoff]
    unique_Nb = np.sort(df_filtered['deltaNb'].unique())

    legend_handles_reps = {}  # To store legend handles for Rep
    legend_handles_Ns = {}

    # Setup the plot layout with two columns per row
    cols = 3
    rows = (len(unique_Nb) + 1) // cols # Calculate number of rows needed
    fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(20, 6 * rows), sharex = False, sharey = False)
    axes = axes.flatten()  # Flatten the axes array for easier indexing
    if len(unique_Nb) == 1: axes = [axes]  # Make it iterable if only one subplot

    for index, dNb in enumerate(unique_Nb):
        ax = axes[index]
        df_nb = df_filtered[df_filtered['deltaNb'] == dNb]
        for (U, rep, Ns), group in df_nb.groupby(['U', 'Rep', 'Ns']):
            marker = marker_dict[rep]
            color = color_dict[Ns]
            scaled_energy = group['Energy'] * np.sqrt(Ns)
            ax.scatter(group['U'], scaled_energy, marker=marker, edgecolors=color, facecolors='none', s=100, label=f'{rep} (Ns={Ns})')
            # Create legend handles for reps if not already created
            if rep not in legend_handles_reps:
                handle_rep = mlines.Line2D([], [], marker=marker, linestyle='None', markeredgecolor='black',
                                            markerfacecolor='none', markersize=10, label=rep)
                legend_handles_reps[rep] = handle_rep
            if Ns not in legend_handles_Ns:
                handle_Ns = mlines.Line2D([], [], marker='s', linestyle='None', color = color,
                                            markersize=10, label=f'Ns = {Ns}')
                legend_handles_Ns[Ns] = handle_Ns

        ax.text(0.80, 0.10, f'ΔNb = {dNb}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        ax.grid(True)
    
    # Adjust unused subplots if any
    for j in range(index + 1, len(axes)):
        axes[j].axis('off')

    # Adjust layout
    fig.subplots_adjust(left=0.05, right=0.99, top=0.94, bottom=0.06, wspace=0.08, hspace=0.05)
    # Set a common x and y label for the figure
    fig.text(0.01, 0.5, r'$\Delta E_{Ns} \sqrt{Ns}$', va='center', rotation='vertical', fontsize=12)
    fig.text(0.45, 0.01, 'U/t', va='center', rotation='horizontal', fontsize=12)
    # Add one legend to the figure
    legend_reps = fig.legend(handles = list(legend_handles_reps.values()), loc='upper right', ncol = 7)    
    plt.gca().add_artist(legend_reps)  # This keeps the first legend visible when adding the second
    plt.legend(handles=list(legend_handles_Ns.values()), loc='lower center')  # bbox_to_anchor=(0.1, 0.2)
    fig.suptitle(f'Nmax = {Nmax}, t = {t}', fontsize= 16, x=0.15, y=0.98, ha='left')
    
    # Decide to save the plot as a file or display it
    directory = 'plots'  # Directory path
    resolution = 400
    if not os.path.exists(directory): os.makedirs(directory)  # Create the directory if it doesn't exist
    
    if save_plot:
        filename = f'{directory}/{lattice}_spectMulti_Nmax{Nmax}_t{t}_{"_Rep_GK" if Rep_cutoff else "_Rep_All"}.png'
        plt.savefig(filename, dpi=resolution)
    else:
        plt.show()


def analyze_Nmax(Nmax_vals, t, lattice, U, Ns, energy_cutoff=20):

    results_df = pd.DataFrame()
    for Nmax in Nmax_vals:
        df = extract_energies(Nmax, t, lattice)
        df = chemPot_correction(df)
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