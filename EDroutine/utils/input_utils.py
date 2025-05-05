import numpy as np
import warnings

def parse_input_file(filename):
    """
    Parses the input file and returns a dictionary containing key-value pairs.
    
    Recognizes:
      - String values (e.g., "Matrix/Spin.X")
      - Numerical values (int, float)
      - Ignores comments and empty lines.
      
    Example input:
      LatticeDescriptionFile="LatticeFiles/TC.mixedfields.20.lat";
      hxFilename="Matrix/Spin.X";
      hxFactor=-0.90;
      Nevals=1;
      OutFile="OutToric/Energies.dat";
    """
    params = {}

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Ignore empty lines and comments
            if not line or line.startswith("#"):
                continue

            # Split key and value
            if "=" not in line:
                raise ValueError(f"Invalid input line (missing '='): {line}")

            key, value = line.split("=", 1)  # Split at first '='
            key = key.strip()

            # Remove trailing semicolon
            value = value.strip().strip(";")

            # Detect numeric values
            if value.replace(".", "", 1).lstrip("-").isdigit():
                if "." in value:
                    value = np.float64(value)  # Convert to float
                else:
                    value = np.int64(value)    # Convert to int
            else:
                # Strip quotes for string values
                value = value.strip('"')

            params[key] = value  # Store key-value pair

    return params

def parse_input_string(input_text):
    """
    Parses an input string and returns a dictionary of parameters.
    Assumes the format: key="value";
    """
    params = {}
    lines = input_text.strip().split("\n")
    for line in lines:
        if line.startswith("#"): 
            continue
        if "=" in line:
            key, value = line.split("=")
            key = key.strip()
            value = value.strip().strip(";").strip('"')  # Remove trailing `;` and quotes
            try:
                # Convert numeric values if possible
                if "." in value or "e" in value.lower():
                    value = float(value)
                else:
                    value = int(value)
            except ValueError:
                pass  # Keep as string if conversion fails
            params[key] = value
    return params

def parse_tags_factors_matrices(params):
    """
    Parses the input parameters to extract:
    - Valid interaction tags
    - Interaction matrices (filenames)
    - Prefactors for each interaction tag
    
    :input params: Dictionary of input parameters.
    :return valid_tags, matrices: Dictionaries of valid interaction tags and matrices.
    """
    # Parsing of valid tags, matrices and prefactors
    valid_tags = {}  # Stores tags and their corresponding filenames and prefactors
    matrices = {}    # Stores interaction matrices
    for key, value in params.items():
        interaction_tag = key.replace("Filename", "").replace("Factor", "")  # Extract tag (e.g., "hx", "Jz")
        if "Filename" in key:
            valid_tags[interaction_tag] = {"filename": value}  # Store filename
            matrices[interaction_tag] = parse_matrix_file(value)  # Parse matrix immediately
        elif "Factor" in key:
            if interaction_tag not in valid_tags:
                valid_tags[interaction_tag] = {}  # Ensure dictionary entry exists
            valid_tags[interaction_tag]["factor"] = value  # Store prefactor

    return valid_tags, matrices

def parse_lattice_file(filename, valid_tags):
    """
    Parses the lattice file to extract:
    - Sites and their coordinates
    - Interactions (filtered by valid tags)
    
    Issues warnings if:
    - The lattice file contains interactions that are NOT in valid_tags.
    - Some valid_tags are NOT found in the lattice file.

    :param filename: Path to the lattice file.
    :param valid_tags: Set of interaction tags expected from input parameters.
    :return: sites (dict), interactions (dict)
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    sites = {}
    interactions = {}
    num_interactions = None
    reading_sites = False
    reading_interactions = False
    found_tags = set()  # Track interaction tags found in the lattice file

    for line in lines:
        line = line.strip()

        # Ignore comments and empty lines
        if line.startswith("#") or line == "":
            continue

        # Read dimension
        if line.startswith("[Dimension]"):
            dimension = int(line.split('=')[1])
            continue

        # Read number of sites
        if line.startswith("[Sites]"):
            num_sites = int(line.split('=')[1])
            reading_sites = True
            continue

        # Read number of interactions
        if line.startswith("[Interactions]"):
            num_interactions = int(line.split('=')[1])
            reading_sites = False
            reading_interactions = True
            continue

        # Read site coordinates
        if reading_sites:
            coords = tuple(map(float, line.split()))  # Convert coordinates to a tuple of floats
            sites[len(sites)] = coords  # Use the current length of the dictionary as the index
            continue

        # Read interactions
        if reading_interactions:
            parts = line.split()
            interaction_type = parts[0]  # Example: "CORE1BODY" or "CORE2BODY"
            tag = parts[1]  # Example: "hx", "Spin.Ising", etc.
            site_indices = tuple(map(int, parts[2:]))  # Extract site indices

            # Track found interaction tags
            found_tags.add(tag)

            # Only keep interactions that match valid tags
            if tag in valid_tags:
                interaction_index = len(interactions)
                interactions[interaction_index] = {
                    "type": interaction_type,
                    "sites": site_indices,
                    "tag": tag
                }
            else:
                warnings.warn(f"Warning: Interaction tag '{tag}' in lattice file is NOT specified in input parameters.")

            # Stop after reading all interactions
            if len(interactions) == num_interactions:
                break

    # Check for missing tags
    print(f"Found interaction tags in lattice file: {found_tags}")
    missing_tags = set(valid_tags.keys()) - found_tags
    if missing_tags:
        warnings.warn(f"Warning: The following expected interaction tags were NOT found in the lattice file: {missing_tags}")

    return sites, interactions

def parse_matrix_file(filename):

    """Parses an interaction matrix file and returns a dictionary mapping input states to (output state, value)."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    num_elements = int(lines[0].strip())  # First line: Number of nonzero elements
    matrix_elements = {
        (int(line.split()[0]), int(line.split()[1])): np.complex128(line.split()[2])
        for line in lines[1:num_elements+1]
    }
    
    return matrix_elements

def parse_movelist(n_sites):
    """Parses the movelist file for n_sites lattice and returns a list of integers."""
    with open(f"Matrix/movelist.{n_sites}", 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("#") or not line:
                continue  # Skip comments and empty lines
            movelist = eval(line)  # Convert the string representation of the list to an actual list
            if len(movelist) != n_sites:
                raise ValueError(f"Movelist length {len(movelist)} does not match the number of sites {n_sites}.")
            return movelist
    raise ValueError("No valid movelist found in the file.")


def write_to_file(output_file:str , mode:str , observables, observable_name:str):
    """
    Writes observables to a file.
    :param output_file: Path to the output file.
    :param mode: Mode to open the file ('w' for write, 'a' for append).
    :param observables: Array-like list of observables to write.
    :param observable_name: Name of the observable for labeling.
    """

    with open(output_file, mode) as file:
        if observable_name == "GroundState" or observable_name=="groundstate" or observable_name=="ground_state" or observable_name=="Ground_state":
            for i,value in enumerate(observables):
                file.write(f"{value} ")
            print(f"Ground state written to {output_file}")
        else :
            for i,value in enumerate(observables):
                file.write(f"{observable_name}[{i}] = {value:.16f}\n")
                # file.write(f"{observable_name}[{i}] = Re {value.real:16f} , Im {value.imag:16f}\n")
        print(f"{observable_name} list written to {output_file}")