import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os
import re
import glob
import math



def extract_correlations(Nmax, Ns, t, U_max, lattice, operator, tot_time, Dt):
    U_max_format = "{:.2f}".format(U_max)
    Dt= float(Dt)
    # Calculate the number of decimal places to represent Dt without scientific notation
    if Dt != 0:
        decimal_places = max(-int(math.floor(math.log10(abs(Dt)))), 1)
    else:
        decimal_places = 1
    # Format Dt with dynamic decimal places
    Dt_format = f'{{:.{decimal_places}f}}'
    Dt_str = Dt_format.format(Dt)
    tot_time= float(tot_time)
    directory = f"/Users/pietro/Runtime/OutTime"
    path_pattern = f"TE.Bosons_{Nmax}.{lattice}.{Ns}*.t.{t}.U.{U_max_format}.T.{tot_time}.DT.{Dt_str}.dat"
    files = glob.glob(os.path.join(directory, path_pattern))
    if not files:
        print(f"No files found for Nmax={Nmax}, Ns={Ns}, t={t}, U={U_max_format}, totTime={tot_time}, Dt={Dt}, lattice={lattice}")
        return pd.DataFrame()
    print("Reading file:\n", files, "\n")
    all_data = []
    DU = (1.0) / int(tot_time/Dt)
    for file_path in files:
        with open(file_path, 'r') as file:
                for line in file:

                    if line.startswith(f'Corr1Body_{operator}'):
                        # Extract the time value
                        time_match = re.search(f'Corr1Body_{operator}\((-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\)', line)
                        time = float(time_match.group(1)) if time_match else None
                        U = round(float(U_max) *(1.0 - int(time/Dt) * DU),2) if time_match else None
                        # Extract the site value
                        site_match = re.search(f'\[(\d+)\]', line)
                        site = int(site_match.group(1)) if site_match else None
                        # Extract the real and imaginary parts of the correlation value
                        corr_match = re.search(f'=\s*\((-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?),\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\)', line)
                        corr2 = float(corr_match.group(1)) if corr_match else None
                        imaginarypart = float(corr_match.group(2)) if corr_match else None
                        all_data.append({'time': time, 'U':U, 'site1': site, 'site2': site, 'Corr2Body': corr2, 'ImagPart': imaginarypart})

                    if line.startswith(f'Corr2Body_{operator}'):
                        # Extract the time value
                        time_match = re.search(f'Corr2Body_{operator}\((-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\)', line)
                        time = float(time_match.group(1)) if time_match else None
                        U = round(float(U_max) *(1.0 - int(time/Dt) * DU),2) if time_match else None
                        # Extract the site value
                        site_match = re.search(f'\[(\d+),(\d+)\]', line)
                        site1 = int(site_match.group(1)) if site_match else None
                        site2 = int(site_match.group(2)) if site_match else None
                        # Extract the real and imaginary parts of the correlation value
                        corr_match = re.search(f'=\s*\((-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?),\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\)', line)
                        corr2 = float(corr_match.group(1)) if corr_match else None
                        imaginarypart = float(corr_match.group(2)) if corr_match else None
                        all_data.append({'time': time, 'U':U, 'site1': site1, 'site2': site2, 'Corr2Body': corr2, 'ImagPart': imaginarypart})

    df = pd.DataFrame(all_data)

    directory = "/Users/pietro/Runtime/LatticeFiles"
    path_pattern = f"{lattice}.BoseHubbard.{Ns}*.lat"
    files = glob.glob(os.path.join(directory, path_pattern))
    # print("Reading lattice file:\n", files, "\n")
    coords = []
    start_reading_coords = False
    for file_path in files:
        with open(file_path, 'r') as file:
            for line in file:
                if '[Sites]' in line:
                    start_reading_coords = True
                    continue
                if start_reading_coords and '[Interactions]' in line:
                    start_reading_coords = False
                    continue
                if start_reading_coords:
                    if line.strip() and not line.startswith('['):
                        coords.append([float(num) for num in line.split()])

    coords = np.array(coords)
    for idx, row in df.iterrows():
        site1 = int(row['site1'])
        site2 = int(row['site2'])
        df.at[idx, 'x1'] = coords[site1,0]
        df.at[idx, 'y1'] = coords[site1,1]
        df.at[idx, 'x2'] = coords[site2,0]
        df.at[idx, 'y2'] = coords[site2,1]
        # print(df)
    
    return df

def read_lattice_file(lattice, Ns):
    directory = "/Users/pietro/Runtime/LatticeFiles"
    path_pattern = f"{lattice}.BoseHubbard.{Ns}*.lat"
    if lattice == "kagome":
        Nt = Ns + (Ns/3)
        path_pattern = f"triangular.BoseHubbard.{round(Nt)}*.lat"
    files = glob.glob(os.path.join(directory, path_pattern))
    # print("Reading lattice file:\n", files, "\n")
    
    k_points = []
    marked_with_star = []
    start_reading_k = False
    
    for file_path in files:
        with open(file_path, 'r') as file:
            for line in file:
                if '# K points (K wedge marked with *):' in line:
                    start_reading_k = True
                    continue
                if start_reading_k and '# High Symmetry Points:' in line:
                    start_reading_k = False
                    continue
                if start_reading_k:
                    if line.startswith('# ['):
                        numbers = line[line.index('[') + 1:line.index(']')]
                        is_marked = '*' in line
                        k_points.append(np.array([float(num) for num in numbers.split()]))
                        marked_with_star.append(is_marked)

    df_k = pd.DataFrame(k_points, columns=['kx', 'ky'])
    df_k['marked_with_star'] = marked_with_star

    if lattice == "kagome":
        df_k['kx'] *= 2
        df_k['ky'] *= 2

    # Reference values (could be extended or modified)
    reference_values = {
        'Gamma': (0.0, 0.0),
        'K': (4.1887902047863905, 0.0),
        'M': (3.141592653589793, 1.8137993642342178),
        'X': (2.0943951023931953, 0.0),
        'Y': (2.0943951023931953, 1.2091995761561452),
        'C1.1': (1.495996501709425, 0.518228389781205),
        'C1.2': (2.692793703076965, -0.5182283897812052),
        'C1.3': (2.99199300341885, 1.03645677956241),
        'M*': (2 * 3.141592653589793, 0.0),
        'Gamma*': (2 * 3.141592653589793, 2 * 1.8137993642342178),
        'K*': (2 * 4.1887902047863905, 0.0),
        'C1.1*': (2* 1.346396851538483, 2* 0.2591141948906026),
        'C1.2*': (2* 2.2439947525641384,2* -0.7773425846718076),
        'C1.3*': (2* 2.692793703076966, 2* 0.5182283897812052),
        'C1.4*': (2* 3.590391604102621, 2* -0.5182283897812048),
    }
    # print(reference_values)
    
    # Calculate the distance to each reference point and assign labels
    def assign_label(row):
        if row['marked_with_star']:
            min_distance = float('inf')
            label = None
            for key, value in reference_values.items():
                distance = np.sqrt((row['kx'] - value[0])**2 + (row['ky'] - value[1])**2)
                if distance < min_distance:
                    min_distance = distance
                    label = key
            return label
        return None
    
    df_k['label'] = df_k.apply(assign_label, axis=1)

    b1 = (2 * np.pi / (np.sqrt(3))) * np.array([np.sqrt(3), -1])
    b2 = (4 * np.pi / (np.sqrt(3))) * np.array([0, 1])

    return df_k


def compute_correlation_length(df_k, reference_label):
    """
    Compute the correlation length for marked points in a dataframe with respect to a specified reference point.
    The DataFrame should contain the momentum coordinates ('kx', 'ky'), a column for S_k values, and a marking
    column to indicate special points.

    Parameters:
    - df_k: DataFrame containing momentum space data, 'kx', 'ky', 'Sk', 'label', and 'marked_with_star'.
    - reference_label: The label of the reference point to compute distances from (e.g., 'Gamma').

    Returns:
    - df_k: Updated DataFrame with distance from reference point and correlation lengths.
    """
    # Step 1: Identify the reference point
    ref_idx = df_k[df_k['label'] == reference_label].index
    if not ref_idx.empty:
        # print("Reference k-point:  ", reference_label)
        ref_point = np.array([df_k.at[ref_idx[0], 'kx'], df_k.at[ref_idx[0], 'ky']])
        S_k_ref = df_k.at[ref_idx[0], 'Sk']

        # Step 2: Compute distances in momentum space from the reference point to all marked points
        df_k['distance'] = np.nan  # Initialize column with NaNs
        for idx, row in df_k.iterrows():
            if row['marked_with_star']:  # Compute distance only for marked points
                kx, ky = row['kx'], row['ky']
                distance = np.linalg.norm([kx - ref_point[0], ky - ref_point[1]])
                df_k.at[idx, 'distance'] = distance

        # Step 3: Find the minimal distance point(s) from the reference point
        # print("Momentum DataFrame:\n", df_k)
        min_distance = df_k['distance'].min()
        min_distance = df_k[df_k.index != ref_idx[0]]['distance'].min()
        min_distance_points = df_k[(df_k['distance'] == min_distance) & (df_k.index != ref_idx[0])].copy()
        # print("Minimal distance k-points:\n", min_distance_points)

        # Step 4: Compute the correlation length for minimal distance point(s)
        df_k['corr_length'] = np.nan  # Initialize column
        for idx in min_distance_points.index:
            S_k_idx = df_k.at[idx, 'Sk']
            if S_k_idx != 0 and min_distance != 0:
                corr_length = np.sqrt((S_k_ref / S_k_idx - 1)) / min_distance
                df_k.at[idx, 'corr_length'] = corr_length
                min_distance_points.at[idx, 'corr_length']= corr_length

        # print("Dataframe with distances and correlation lengths:\n", df_k)
        # print("Minimal distance points:\n", min_distance_points)
    else:
        print(f"No '{reference_label}' point found in the 'label' column.")
    return df_k


def compute_fourier_transform(df_x, df_k):

    if df_x.empty or df_k.empty:
        print("Empty DataFrame provided to compute_fourier_transform.")
        return df_k

    # if (lattice == "kagome"):
    Ns = len(df_x['site1'].unique())
    # print("Ns= ", Ns)
    r1 = np.stack((df_x['x1'], df_x['y1']), axis=1)
    r2 = np.stack((df_x['x2'], df_x['y2']), axis=1)
    R12 = r1 - r2
    C12 = np.array(df_x['Corr2Body'])
    momentum_coords = np.stack((df_k['kx'], df_k['ky']), axis=1)
    S_k = np.zeros(len(momentum_coords), dtype=complex)
    for idx, k in enumerate(momentum_coords):
        # print("k= ", k, "\n")
        for r, c in zip(R12, C12):
            phase = 1j * np.dot(k, r)
            # print("Coords= ", r, "Correlation= ", c, " --> term in the sum: ", c * np.exp(phase), "\n")
            S_k[idx] += c * np.exp(phase)
        # print("S_k = ", S_k[idx], "\n")
        S_k[idx] = S_k[idx]/Ns

    # Check that imaginary part of S_k is zero to machine precision
    if not np.allclose(np.imag(S_k), 0, atol=1e-10):
        print("Error: Non-zero imaginary part detected in Fourier transform calculations.\n")
        print("Structure factors S(k):\n", S_k)
        return  # Stop function execution here

    if np.any(np.real(S_k) < 0):
        print("Error: Found negative Structure Factors!\n")
        print("Structure Factors:\n", S_k)
        return
    
    df_k['Sk'] = np.real(S_k)
    # print("Momentum Correlations:\n", S_k)
    # print("Momentum coordinates:\n", momentum_coords)
    # print("Momentum DataFrame:\n", df_k)
    return df_k


def analyze_structure_factors(lattice, Nmax, Ns, t, U, operator):

    df_x = extract_correlations(Nmax, Ns, t, U, lattice, operator)
    df_k = read_lattice_file(lattice, Ns)
    df_k = compute_fourier_transform(df_x, df_k, lattice)
    
    # Construct the Brillouin zone 
    K1 = np.array([4.1887902047863905, 0.0])
    K2 = np.array([2.0943951023931953, 3.6275987284684357])
    b1 = (2 * np.pi / (np.sqrt(3))) * np.array([np.sqrt(3), -1])
    b2 = (4 * np.pi / (np.sqrt(3))) * np.array([0, 1])
    vertices = np.array([K1, K2, [-K2[0],K2[1]], -K1, -K2, [K2[0],-K2[1]]])
    vertices = np.vstack([vertices, vertices[0]])
    if lattice=="kagome":
        vertices2 = 2 * vertices
        vertices2 = np.vstack([vertices2, vertices2[0]])

    # Create an array for edge colors based on the 'marked_with_star' column
    edge_colors = ['blue' if star else 'black' for star in df_k['marked_with_star']]

    # Plotting
    plt.figure(figsize=(10, 8))
    S_k = np.array(df_k['Sk'])
    kx = np.array(df_k['kx'])
    ky = np.array(df_k['ky'])
    totSum = np.sum(S_k)
    maximum = np.max(S_k)
    magnitude = np.real(S_k)  # Ensure magnitude is real since S_k is complex
    plt.scatter(kx, ky, s=1100,
                c=magnitude, cmap='viridis', edgecolors=edge_colors, 
                linewidths=1.5, vmin=0, vmax=maximum)  # Set colormap scaling here    , vmin=0
    plt.plot(vertices[:, 0], vertices[:, 1], 'r-')
    if lattice=="kagome":
        # plt.plot(vertices2[:, 0], vertices2[:, 1], 'r-')
        # plt.plot(vertices3[:, 0], vertices3[:, 1], 'r-')
        plt.plot(vertices2[:, 0], vertices2[:, 1], 'r-')
    plt.colorbar(label='Magnitude of Fourier Components')
    plt.grid(True)
    plt.title(f'Structure factors of C(0,i) at U={U}, Ns={Ns}')
    plt.xlabel('kx')
    plt.ylabel('ky')
    plt.axis('equal')
    directory = 'plots'   # Directory path
    resolution = 300
    plt.savefig(f'{directory}/{lattice}_StructureFactors_Nmax{Nmax}_Ns{Ns}_t{t}_U{U}.png', dpi = resolution)
    # plt.show()


def analyze_structure_factors_U(lattice, Nmax, Ns_values, t, U_values, operator):
    results = []
    for Ns in Ns_values:
        for U in U_values:
            try:
                df_x = extract_correlations(Nmax, Ns, t, U, lattice, operator)
                if df_x.empty:
                    continue  # Skip this U value if no data is extracted
                df_k = read_lattice_file(lattice, Ns)
                if df_k.empty:
                    continue  # Skip this U value if no lattice file data is available

                df_k = compute_fourier_transform(df_x, df_k, lattice)
                for idx, row in df_k.dropna(subset=['label']).iterrows():
                    results.append({
                        'U': U,
                        'Ns': Ns,
                        'Sk': row['Sk'],
                        'label': row['label']
                    })
            except Exception as e:
                print(f"Failed to process U={U} for Ns={Ns}: {e}")

    results_df = pd.DataFrame(results)
    if results_df.empty:
        print("No data was processed.")
        return

    plt.figure(figsize=(12, 8))
    markers = ['o', 's', '^', 'v', '>', '<', 'p', '*', 'h', 'x']  # List of markers for different labels
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # Basic color list

    # Ensure there are enough markers and colors
    if len(results_df['label'].unique()) > len(markers):
        raise ValueError("Not enough markers specified for the number of labels.")
    if len(results_df['Ns'].unique()) > len(colors):
        raise ValueError("Not enough colors specified for the number of Ns values.")
    
    ns_handles = {}  # For storing Ns legend handles
    label_handles = {}  # For storing label legend handles
    # Create dictionaries for markers and colors
    marker_dict = {rep: markers[i % len(markers)] for i, rep in enumerate(sorted(results_df['label'].unique()))}
    color_dict = {ns: colors[i % len(colors)] for i, ns in enumerate(sorted(results_df['Ns'].unique()))}

    for i, ns in enumerate(sorted(results_df['Ns'].unique())):
        ns_subset = results_df[results_df['Ns'] == ns]
        color = color_dict[ns]
        for j, label in enumerate(sorted(ns_subset['label'].unique())):
            marker = marker_dict[label]
            label_subset = ns_subset[ns_subset['label'] == label]
            plt.plot(label_subset['U'], label_subset['Sk'], marker=marker, markeredgecolor=color, markerfacecolor='none', markeredgewidth=2, markersize=10,
                     linestyle='', color=color, label=None)
            # Create legend handles for reps if not already created
            if label not in label_handles:
                handle_rep = mlines.Line2D([], [], marker=marker, linestyle='None', markeredgecolor='black',
                                        markerfacecolor='none', markersize=10, label=label)
                label_handles[label] = handle_rep

        # Create legend handles for Nb if not already created
        if ns not in ns_handles:
            handle_ns = mlines.Line2D([], [], color=color, marker='s', linestyle='None',
                                      markersize=10, markeredgewidth=2, label=f'Ns = {ns}')
            ns_handles[ns] = handle_ns

    plt.title(f'Structure Factors for Nmax={Nmax}, t={t} and Nb=Ns')
    plt.xlabel('U/t')
    plt.ylabel('S(k)')
    # Position the legends on the plot
    legend_reps = plt.legend(handles=list(label_handles.values()), title='k-point:', loc='upper left', bbox_to_anchor=(1.005, 1))
    plt.gca().add_artist(legend_reps)  # This keeps the first legend visible when adding the second
    plt.legend(handles=list(ns_handles.values()), title='Lattice Size:', loc='upper left', bbox_to_anchor=(1.005, 0.4))
    plt.grid(True)
    plt.subplots_adjust(left=0.08, right=0.85, bottom=0.08, top=0.95)
    directory = 'plots'   # Directory path
    resolution = 400
    plt.savefig(f'{directory}/{lattice}_StructureFactorsVsU_Nmax{Nmax}Ns{Ns_values}_t{t}.png', dpi = resolution)
    plt.show()



def analyze_structure_factors_Nmax(lattice, Ns_values, t, U, operator):
    results = []
    for Ns in Ns_values:
        for Nmax in np.arange(2, Ns+1):
            try:
                df_x = extract_correlations(Nmax, Ns, t, U, lattice, operator)
                if df_x.empty:
                    continue  # Skip this Nmax value if no data is extracted
                df_k = read_lattice_file(lattice, Ns)
                if df_k.empty:
                    continue  # Skip this Ns value if no lattice file data is available

                df_k = compute_fourier_transform(df_x, df_k, lattice)
                for idx, row in df_k.dropna(subset=['label']).iterrows():
                    results.append({
                        'Nmax': Nmax,
                        'Ns': Ns,
                        'Sk': row['Sk'],
                        'label': row['label']
                    })
            except Exception as e:
                print(f"Failed to process U={U} for Ns={Ns}: {e}")

    results_df = pd.DataFrame(results)
    if results_df.empty:
        print("No data was processed.")
        return

    plt.figure(figsize=(12, 8))
    markers = ['o', 's', '^', 'v', '>', '<', 'p', '*', 'h', 'x']  # List of markers for different labels
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # Basic color list

    # Ensure there are enough markers and colors
    if len(results_df['label'].unique()) > len(markers):
        raise ValueError("Not enough markers specified for the number of labels.")
    if len(results_df['Ns'].unique()) > len(colors):
        raise ValueError("Not enough colors specified for the number of Ns values.")
    
    ns_handles = {}  # For storing Ns legend handles
    label_handles = {}  # For storing label legend handles
    # Create dictionaries for markers and colors
    marker_dict = {rep: markers[i % len(markers)] for i, rep in enumerate(sorted(results_df['label'].unique()))}
    color_dict = {ns: colors[i % len(colors)] for i, ns in enumerate(sorted(results_df['Ns'].unique()))}

    for i, ns in enumerate(sorted(results_df['Ns'].unique())):
        ns_subset = results_df[results_df['Ns'] == ns]
        color = color_dict[ns]
        for j, label in enumerate(sorted(ns_subset['label'].unique())):
            marker = marker_dict[label]
            label_subset = ns_subset[ns_subset['label'] == label]
            plt.plot(label_subset['Nmax'], label_subset['Sk'], marker=marker, markeredgecolor=color, markerfacecolor='none', markeredgewidth=2, markersize=10,
                     linestyle='', color=color, label=None)
            # Create legend handles for reps if not already created
            if label not in label_handles:
                handle_rep = mlines.Line2D([], [], marker=marker, linestyle='None', markeredgecolor='black',
                                        markerfacecolor='none', markersize=10, label=label)
                label_handles[label] = handle_rep

        # Create legend handles for Ns if not already created
        if ns not in ns_handles:
            handle_ns = mlines.Line2D([], [], color=color, marker='s', linestyle='None',
                                      markersize=10, markeredgewidth=2, label=f'Ns = {ns}')
            ns_handles[ns] = handle_ns

    plt.title(f'Structure Factors for U={U}, t={t} and Nb=Ns')
    plt.xlabel('Nmax')
    plt.ylabel('S(k)')
    # Position the legends on the plot
    legend_reps = plt.legend(handles=list(label_handles.values()), title='k-point:', loc='upper left', bbox_to_anchor=(1.005, 1))
    plt.gca().add_artist(legend_reps)  # This keeps the first legend visible when adding the second
    plt.legend(handles=list(ns_handles.values()), title='Lattice Size:', loc='upper left', bbox_to_anchor=(1.005, 0.4))
    plt.grid(True)
    plt.subplots_adjust(left=0.08, right=0.85, bottom=0.08, top=0.95)
    directory = 'plots'   # Directory path
    resolution = 400
    plt.savefig(f'{directory}/{lattice}_StructureFactorsVsNmax_U{U}Ns{Ns_values}_t{t}.png', dpi = resolution)
    plt.show()




def analyze_corr_length(lattice, Nmax, Ns_values, t, U_max, tot_times, Dts, operator, reference_momentum):
    results = []
    for Ns in Ns_values:
        for tot_time in tot_times:
            for Dt in Dts:
                try:
                    df_x = extract_correlations(Nmax, Ns, t, U_max, lattice, operator, tot_time, Dt)
                    if df_x.empty:
                        continue
                    df_k = read_lattice_file(lattice, Ns)
                    # Iterate over unique times within df_x
                    unique_times = df_x['time'].unique()
                    for time in unique_times:
                        # Create a fresh copy of df_k for each time iteration
                        df_k_temp = df_k.copy()
                        # Create a subset of df_x for the specific time
                        df_x_subset = df_x[df_x['time'] == time]
                        # Apply the Fourier transform computation to this subset
                        df_k_temp = compute_fourier_transform(df_x_subset, df_k_temp)
                        # Compute correlation length for the transformed data
                        df_k_temp = compute_correlation_length(df_k_temp, reference_momentum)
                        # Collect results with additional metadata about the computation
                        for idx, row in df_k_temp.dropna(subset=['corr_length']).iterrows():
                            results.append({
                                'time': time,
                                'tot_time': tot_time,
                                'Dt': Dt,
                                'Ns': Ns,
                                'Correlation Length': row['corr_length'],
                                'label': row['label']
                            })
                except Exception as e:
                    print(f"Failed to process totTime={tot_time} for Dt={Dt}: {e}")
                    
    results_df = pd.DataFrame(results)
    return results_df


def plot_corr_length(results_df, lattice, Nmax, t, U_max, reference_momentum, logscale=False):
    plt.figure(figsize=(12, 8))
    markers = ['o', 's', '^', 'v', '>', '<', 'p', '*', 'h', 'x']  # List of markers for different labels
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # Basic color list

    # Ensure there are enough markers and colors
    if len(results_df['label'].unique()) > len(markers):
        raise ValueError("Not enough markers specified for the number of labels.")
    if len(results_df['Ns'].unique()) > len(colors):
        raise ValueError("Not enough colors specified for the number of Ns values.")
    
    ns_handles = {}  # For storing Ns legend handles
    label_handles = {}  # For storing label legend handles
    # Create dictionaries for markers and colors
    marker_dict = {rep: markers[i % len(markers)] for i, rep in enumerate(sorted(results_df['label'].unique()))}
    color_dict = {ns: colors[i % len(colors)] for i, ns in enumerate(sorted(results_df['Ns'].unique()))}

    for i, ns in enumerate(sorted(results_df['Ns'].unique())):
        ns_subset = results_df[(results_df['Ns'] == ns) & (results_df['tot_time']-results_df['time']<1e-6)]
        color = color_dict[ns]
        min_dt_per_time = ns_subset.groupby('tot_time')['Dt'].min().reset_index()
        ns_subset = pd.merge(ns_subset, min_dt_per_time, on=['tot_time', 'Dt'])

        for j, label in enumerate(sorted(ns_subset['label'].unique())):
            marker = marker_dict[label]
            label_subset = ns_subset[(ns_subset['label'] == label)]
            plt.plot(label_subset['tot_time'], label_subset['Correlation Length'], marker=marker, markeredgecolor=color, markerfacecolor='none', markeredgewidth=2, markersize=10,
                    linestyle='--', color=color, label=None)
            # Create legend handles for reps if not already created
            if label not in label_handles:
                handle_rep = mlines.Line2D([], [], marker=marker, linestyle='None', markeredgecolor='black',
                                        markerfacecolor='none', markersize=10, label=label)
                label_handles[label] = handle_rep
        # Create legend handles for Nb if not already created
        if ns not in ns_handles:
            handle_ns = mlines.Line2D([], [], color=color, marker='s', linestyle='None',
                                      markersize=10, markeredgewidth=2, label=f'Ns = {ns}')
            ns_handles[ns] = handle_ns

    
    plt.title(f'Correlation Length for Nmax={Nmax}, t={t} and Nb=Ns (reference momentum {reference_momentum})')
    plt.xlabel('ramp time')
    plt.ylabel('Correlation Length')
    if logscale: plt.xscale('log')
    # Position the legends on the plot
    legend_reps = plt.legend(handles=list(label_handles.values()), title='Nearest k-Points', loc='upper left', bbox_to_anchor=(1.005, 1))
    plt.gca().add_artist(legend_reps)  # This keeps the first legend visible when adding the second
    plt.legend(handles=list(ns_handles.values()), title='Lattice Size', loc='upper left', bbox_to_anchor=(1.005, 0.4))
    plt.grid(True)
    plt.subplots_adjust(left=0.08, right=0.85, bottom=0.08, top=0.95)
    directory = 'plots'   # Directory path
    resolution = 400
    Ns_vals = sorted(results_df['Ns'].unique())
    Ns_vals_str = '_'.join(map(str, Ns_vals))
    plt.savefig(f'{directory}/{lattice}_CorrelationLenghtVsTramp_U{U_max}_Nmax{Nmax}_Ns{Ns_vals_str}_t{t}.png', dpi = resolution)
    plt.show()



def analyze_Dt_convergence(lattice, Nmax, Ns_vals, t, U_max, tot_times, Dts, operator, reference_momentum):

    results = []
    for Ns in Ns_vals:
        for tot_time in tot_times:
            for Dt in Dts:
                try:
                    df_x = extract_correlations(Nmax, Ns, t, U_max, lattice, operator, tot_time, Dt)
                    if df_x.empty:
                        continue
                    df_k = read_lattice_file(lattice, Ns)
                    # Iterate over unique times within df_x
                    unique_times = df_x['time'].unique()
                    for time in unique_times:
                        # Create a fresh copy of df_k for each time iteration
                        df_k_temp = df_k.copy()
                        # Create a subset of df_x for the specific time
                        df_x_subset = df_x[df_x['time'] == time]
                        # Apply the Fourier transform computation to this subset
                        df_k_temp = compute_fourier_transform(df_x_subset, df_k_temp)
                        # Compute correlation length for the transformed data
                        df_k_temp = compute_correlation_length(df_k_temp, reference_momentum)
                        # Collect results with additional metadata about the computation
                        for idx, row in df_k_temp.dropna(subset=['corr_length']).iterrows():
                            results.append({
                                'Ns': Ns,
                                'time': time,
                                'Dt': Dt,
                                'tot_time': tot_time,
                                'Correlation Length': row['corr_length'],
                                'label': row['label']
                            })
                except Exception as e:
                    print(f"Failed to process totTime={tot_time} for Dt={Dt}: {e}")
                    
    results_df = pd.DataFrame(results)
    return results_df


def plot_Dt_convergence(results_df, lattice, Nmax, t, U_max, reference_momentum, logscale=False):
    plt.figure(figsize=(12, 8))
    markers = ['o', 's', '^', 'v', '>', '<', 'p', '*', 'h', 'x']  # List of markers for different labels
    colors = ['k', 'b', 'g', 'r', 'c', 'm', 'y', 'lime', 'deepskyblue', 'sienna', 'blueviolet']  # Basic color list
    lines = ['--', ':', '-.']

    # Ensure there are enough markers and colors
    if len(results_df['Ns'].unique()) > len(markers):
        raise ValueError("Not enough markers specified for the number of labels.")
    if len(results_df['Dt'].unique()) > len(colors):
        raise ValueError("Not enough colors specified for the number of Ns values.")
    
    ns_handles = {}  # For storing Ns legend handles
    label_handles = {}  # For storing label legend handles
    Dt_handles = {}

    # Create dictionaries for markers and colors
    marker_dict = {ns: markers[i % len(markers)] for i, ns in enumerate(sorted(results_df['Ns'].unique()))}
    color_dict = {dt: colors[i % len(colors)] for i, dt in enumerate(sorted(results_df['Dt'].unique()))}
    lines_dict = {label: lines[i % len(lines)] for i, label in enumerate(sorted(results_df['label'].unique()))}
    
    for i, Dt in enumerate(sorted(results_df['Dt'].unique())):
        Dt_subset = results_df[results_df['Dt'] == Dt]
        color = color_dict[Dt]
        for j, Ns in enumerate(sorted(Dt_subset['Ns'].unique())):
            marker = marker_dict[Ns]
            Ns_subset = Dt_subset[Dt_subset['Ns'] == Ns]
            for l, label in enumerate(sorted(Ns_subset['label'].unique())):
                label_subset = Ns_subset[Ns_subset['label'] == label] 
                line = lines_dict[label]
                for m, tot_time in enumerate(sorted(label_subset['tot_time'].unique())):
                    tot_time_subset = label_subset[label_subset['tot_time'] == tot_time]
                    plt.plot(tot_time_subset['time'], tot_time_subset['Correlation Length'], marker=marker, markeredgecolor=color, markerfacecolor='none', markeredgewidth=1, markersize=15,
                        linestyle=line, color=color, label=None)
                if label not in label_handles:
                    handle_label = mlines.Line2D([], [], marker= 'None', linestyle=line, color='black', label=label)
                    label_handles[label] = handle_label

            # Create legend handles for reps if not already created
            if Ns not in ns_handles:
                handle_ns = mlines.Line2D([], [], marker=marker, linestyle='None', markeredgecolor='black',
                                        markerfacecolor='none', markersize=10, label=Ns)
                ns_handles[Ns] = handle_ns

        # Create legend handles for Nb if not already created
        if Dt not in Dt_handles:
            handle_Dt = mlines.Line2D([], [], color=color, marker='s', linestyle='None',
                                      markersize=10, markeredgewidth=2, label=f'Dt = {Dt}')
            Dt_handles[Dt] = handle_Dt

    
    plt.title(f'Correlation Length for Nmax={Nmax}, t={t} (reference momentum {reference_momentum})')  # and Nb=Ns={Ns}
    plt.xlabel('time')
    plt.ylabel('Correlation Length')
    if logscale: plt.xscale('log')
    # plt.yscale('log')
    # Position the legends on the plot
    legend_reps = plt.legend(handles=list(label_handles.values()), title='Nearest k-Points', loc='upper left', bbox_to_anchor=(1.005, 1))
    plt.gca().add_artist(legend_reps)  # This keeps the first legend visible when adding the next
    legend_Ns = plt.legend(handles=list(ns_handles.values()), title='Ns values', loc='upper left', bbox_to_anchor=(1.005, 0.6))
    plt.gca().add_artist(legend_Ns)  # This keeps the first legend visible when adding the next
    plt.legend(handles=list(Dt_handles.values()), title='Dt', loc='upper left', bbox_to_anchor=(1.005, 0.3))
    plt.grid(True)
    plt.subplots_adjust(left=0.08, right=0.85, bottom=0.08, top=0.95)
    directory = 'plots'   # Directory path
    resolution = 400
    tot_times = sorted(results_df['tot_time'].unique())
    tot_times_str = '_'.join(map(str, tot_times))
    Ns_vals = sorted(results_df['Ns'].unique())
    Ns_vals_str = '_'.join(map(str, Ns_vals))
    plt.savefig(f'{directory}/{lattice}_DtConvergence_rampU{U_max}_totTime[{tot_times_str}]_Nmax{Nmax}_Ns[{Ns_vals_str}]_t{t}.png', dpi = resolution)
    plt.show()










