import numpy as np
import numba
from scipy.linalg import eig

@numba.njit(cache=True)
def readsite(state, s):
    """
    Read the value of the bit at position 's' from the integer representation 'state'.

    - 'state': integer representation of a computational basis state;
    - 's': index of the bit to read (0-indexed).
    """
    # Extract the bit at position 's' from 'state' by masking and shifting
    return (state & (1 << s)) >> s


# Implement the action of the Tranlation operator on a state
@numba.njit(cache=True)
def permutation(state, sites, movelist):
    """
    Apply a lattice sites translation operation on a computational basis state.

    - 'state': integer representation of a computational basis state;
    - 'sites': number of sites in the system;
    - 'movelist': list of site indexes representing the translation operation on the lattice.
    """
    tstate = 0
    for s in range(sites):
        tstate |= readsite(state, s) << movelist[s]
    return tstate

@numba.njit(cache=True)
def translation_op(psi, sites, movelist):
    out_psi = np.zeros_like(psi)
    for i in range(len(psi)):
        out_psi[permutation(i, sites, movelist)] = psi[i]
    return out_psi


# --- Optimized hadamard_op using FWHT and Numba ---
@numba.njit(cache=True)
def hadamard_op(psi, sites):
    """
    Apply the global Hadamard gate (H_0 ⊗ H_1 ⊗ ... ⊗ H_{sites-1})
    to the state vector psi using an efficient in-place
    Fast Walsh-Hadamard Transform (FWHT) algorithm. Optimized with Numba.

    Args:
        psi (np.ndarray): Input state vector (float64 or complex128).
                          The operation is performed in-place.
        sites (int): Number of sites (qubits). len(psi) must be 2**sites.

    Returns:
        np.ndarray: The transformed state vector (modified in-place).
    """
    n_states = 1 << sites  # Calculate 2**sites efficiently

    if len(psi) != n_states:
        # Cannot proceed reliably if sizes don't match
        # Returning the original psi might be misleading. Consider raising ValueError.
        raise ValueError("Length of state vector psi must be 2**sites.")


    # --- FWHT Algorithm (Iterative, In-Place) ---
    # Iterate through each qubit/site 's'
    for s in range(sites):
        # 'd' determines the stride/distance between paired elements for this stage
        d = 1 << s  # d = 2**s

        # Iterate through the "butterfly" groups
        # Increment by 2*d because each butterfly involves pairs separated by d
        for block_start in range(0, n_states, 2 * d):
            # Iterate through the pairs within the current block
            for i in range(block_start, block_start + d):
                # Indices of the pair to transform
                idx1 = i
                idx2 = i + d

                # Perform the Hadamard butterfly operation IN-PLACE
                # Use temporary variables to avoid overwriting needed values
                val1 = psi[idx1]
                val2 = psi[idx2]
                psi[idx1] = val1 + val2
                psi[idx2] = val1 - val2

    # --- Apply Global Normalization ---
    # Each single Hadamard has 1/sqrt(2). N Hadamards give (1/sqrt(2))**N = 1 / sqrt(2**N)
    norm_factor = np.sqrt(n_states) # Calculate sqrt(2**sites)
    psi /= norm_factor

    return psi # Return the modified psi array


def compute_duality_symmetry(eigenvalues, eigenvectors, n_sites, movelist, tolerance=1e-12):
    """
    Compute the duality symmetry eigenvalues for a given Hamiltonian.

    Args:
        eigenvalues (np.ndarray): Eigenvalues of the Hamiltonian.
        eigenvectors (np.ndarray): Eigenvectors of the Hamiltonian.
        n_sites (int): Number of sites in the system.
        movelist (list): List of site indexes representing the translation operation on the lattice.
        tolerance (float): Tolerance for grouping eigenvalues by (near-)degeneracy.

    Returns:
        None: Prints the duality symmetry eigenvalues.
    """

    # --- Group eigenstates by (near-)degeneracy ---
    groups = []
    current_group = [0]
    for i in range(1, len(eigenvalues)):
        if abs(eigenvalues[i] - eigenvalues[i-1]) < tolerance:
            current_group.append(i)
        else:
            groups.append(current_group)
            current_group = [i]
    if current_group:
        groups.append(current_group)
    print(groups)

    # --- Analyze each group ---
    duality_expectations = []
    for group in groups:
        block_size = len(group)
        if block_size == 1:
            # Non-degenerate case
            i = group[0]
            psi = translation_op(eigenvectors[:, i], n_sites, movelist)
            psi = hadamard_op(psi, n_sites)
            duality_expectations.append(np.vdot(eigenvectors[:, i], psi))
        else:
            # Degenerate block: construct D matrix
            O_matrix = np.zeros((block_size, block_size), dtype=np.complex128)
            D_matrix = np.zeros((block_size, block_size), dtype=np.complex128)
            for idx_i, i in enumerate(group):
                for idx_j, j in enumerate(group):
                    O_matrix[idx_i, idx_j] = np.vdot(eigenvectors[:, i], eigenvectors[:, j])
                    Dpsi_j = translation_op(eigenvectors[:, j], n_sites, movelist)
                    Dpsi_j = hadamard_op(Dpsi_j, n_sites)
                    D_matrix[idx_i, idx_j] = np.vdot(eigenvectors[:, i], Dpsi_j)

            # Diagonalize the D matrix in this subspace
            duality_eigs, _ = eig(a = D_matrix, b = O_matrix)
            duality_expectations.extend(np.complex128(duality_eigs))

    return duality_expectations