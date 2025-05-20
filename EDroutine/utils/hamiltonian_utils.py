import numpy as np
import numba
from typing import Dict, List, Tuple, Any, Union
from collections import defaultdict
from numba import njit, prange

import line_profiler
import atexit
from numba.typed import List
from scipy.sparse import csr_matrix

# --- Set up line profiler for performance analysis ---
profile = line_profiler.LineProfiler()
atexit.register(profile.print_stats)



# ---- Numba-compiled helper functions ----

# Non-padded version for the SERIAL approach (used when the sites are not padded)
@numba.njit(cache=True)
def extract_local_state(global_index: int, sites: Tuple[int, ...]) -> int:
    local_state = 0
    for i in range(len(sites)):
        site = sites[i]
        if (global_index >> site) & 1:
            local_state |= (1 << i)
    return local_state

@numba.njit(cache=True)
def modify_global_state(global_index: int, sites: Tuple[int, ...], new_local_state: int) -> int:
    new_global_index = global_index
    for i in range(len(sites)):
        site = sites[i]
        bit_value = (new_local_state >> i) & 1
        if bit_value:
            new_global_index |= (1 << site)
        else:
            new_global_index &= ~(1 << site)
    return new_global_index

# Padded version for the PARALLEL approach (used when the sites are padded to a fixed size)
@numba.njit(cache=True)
def extract_local_state_padded(global_index: int, sites_padded: np.ndarray, actual_arity: int) -> int:
    local_state = 0
    for i in range(actual_arity): # Loop only over actual sites
        site = sites_padded[i]
        if (global_index >> site) & 1:
            local_state |= (1 << i)
    return local_state

@numba.njit(cache=True)
def modify_global_state_padded(global_index: int, sites_padded: np.ndarray, int_type: int, new_local_state: int) -> int:
    new_global_index = np.int64(global_index) # Ensure integer type
    for i in range(int_type): # Loop only over actual sites
        site = sites_padded[i]
        bit_value = (new_local_state >> i) & 1
        mask = np.int64(1 << site)
        if bit_value:
            new_global_index |= mask
        else:
            new_global_index &= ~mask
    return new_global_index



# ---- Builder function for the Numba dict (for the SERIAL approach) ----

# Generic type hint for use in the Numba function signature if preferred,
# though less specific. Numba primarily cares about the type of the *list* passed.
PreprocessedInteractionGeneric = Tuple[
    Tuple[int, ...],         # sites_tuple: Variable length integer tuple
    np.ndarray,              # matrix_elements_array: Nx2 complex array [in, out]
    np.ndarray               # matrix_values_array (complex128, shape=(N,))
]

def build_numba_preprocessed_dict(
                                  interactions: Dict[int, Dict[str, Any]],
                                  matrices: Dict[str, Dict[Tuple[int, int], complex]],
                                  valid_tags: Dict[str, Dict[str, Any]]
                                 ) ->  Dict[int, List[PreprocessedInteractionGeneric]]:

    """
    Converts intermediate dictionaries into a dictionary mapping arity
    to lists of Numba-ready interaction terms.
    """
    preprocessed_dict = defaultdict(list)

    for interaction_details in interactions.values():

        tag = interaction_details["tag"]
        factor = np.float64(valid_tags[tag]["factor"])
        # Check if zero factor 
        if factor == (0.0): continue
        sites_tuple = interaction_details["sites"]
        interaction_type = len(sites_tuple)
        matrix_dict = matrices[tag]
        num_elements = len(matrix_dict)

        matrix_indices_array = np.zeros((num_elements, 2), dtype=np.int64)
        matrix_values_array = np.zeros(num_elements, dtype=np.complex128)

        idx = 0
        for (input_local, output_local), value in matrix_dict.items():

            matrix_indices_array[idx, 0] = input_local
            matrix_indices_array[idx, 1] = output_local
            matrix_values_array[idx] = factor * np.complex128(value)
            idx += 1

        # Create the 4-element tuple
        interaction_term_data: PreprocessedInteractionGeneric = (
            sites_tuple, matrix_indices_array, matrix_values_array
        )
        preprocessed_dict[interaction_type].append(interaction_term_data)

    return dict(preprocessed_dict)


# ------ Builder Function for the Numba-friendly CSR-like List of interaction terms (for the PARALLEL gather approach) ------

# --- Define Type Hint for the FLAT CSR Gather list element ---
# We now store offsets and flattened target/value arrays per term
PreprocessedInteractionTermCSR = Tuple[
                                        np.ndarray,        # sites_padded_array (int64, shape=(MAX_ARITY,))
                                        int,               # actual interaction type int_type
                                        # --- CSR-like structure for off-diagonal elements ---
                                        np.ndarray,        # offdiag_offsets (int64, shape=(max_local_state + 2,))
                                        np.ndarray,        # offdiag_target_locals (int64, shape=(total_nnz,))
                                        np.ndarray,        # offdiag_values (complex128, shape=(total_nnz,))
                                        # --- Optional: Store diagonal elements separately ---
                                        np.ndarray         # diagonal_values (complex128, shape=(max_local_state + 1,))
                                    ]

def build_numba_preprocessed_list(
        interactions: Dict[int, Dict[str, Any]],
        matrices: Dict[str, Dict[Tuple[int, int], float]],
        valid_tags: Dict[str, Dict[str, Any]],
        max_int_type: int,
        sentinel: int = -1
    ) -> List:
    """
    Builds a single flat typed list of Numba-ready interaction terms,
    using fixed-length padding for sites and a CSR-like format
    for efficient gather lookups based on source local state.
    """
    preprocessed_list_py = [] # Temp python list

    for interaction_details in interactions.values():
        tag = interaction_details["tag"]
        factor = float(valid_tags[tag]["factor"])
        if factor == 0.0: continue

        sites_tuple = interaction_details["sites"]
        int_type = len(sites_tuple) # or extract from interaction_details["type"]
        if int_type > max_int_type: raise ValueError(f"Interaction type {int_type}-BODY not supported. Max interaction type: {max_int_type}-BODY")

        sites_padded_array = np.full(max_int_type, sentinel, dtype=np.int64)
        sites_padded_array[:int_type] = sites_tuple

        matrix_dict = matrices[tag]
        max_local_state = (1 << int_type) - 1

        # --- Separate diagonal and off-diagonal elements ---
        diag_values_dict = {}
        offdiag_pairs_dict = defaultdict(list) # Key: source_local, Value: List[(target_local, value)]
        n_offdiag = 0
        for (input_local, output_local), value in matrix_dict.items():
            if input_local == output_local:
                diag_values_dict[input_local] = np.complex128(factor * value)
            else:
                offdiag_pairs_dict[output_local].append((input_local, np.complex128(factor * value)))
                n_offdiag += 1

        # --- Build CSR structure for off-diagonal to enable fast lookup ---
        # Size = max_local_state + 2 for easy range calculation
        offdiag_source_pointers = np.zeros((max_local_state+1)+1, dtype=np.int64)
        offdiag_target_locals = np.zeros(n_offdiag, dtype=np.int64)
        offdiag_values = np.zeros(n_offdiag, dtype=np.complex128)

        current_pointer = 0
        for source_local in range(max_local_state + 1):
            offdiag_source_pointers[source_local] = current_pointer
            if source_local in offdiag_pairs_dict:
                pairs = offdiag_pairs_dict[source_local]
                for target_local, value in pairs:
                    offdiag_target_locals[current_pointer] = target_local
                    offdiag_values[current_pointer] = value
                    current_pointer += 1
        offdiag_source_pointers[max_local_state + 1] = current_pointer # Store end offset

        # --- Build diagonal array ---
        diagonal_values = np.zeros(max_local_state + 1, dtype=np.complex128)
        for local_state, value in diag_values_dict.items():
             diagonal_values[local_state] = value

        # --- Append data for this term to the temp list ---
        term_data: Tuple = ( # Use basic tuple for appending to Python list
            sites_padded_array,
            int_type,
            offdiag_source_pointers, 
            offdiag_target_locals, 
            offdiag_values,
            diagonal_values
        )
        preprocessed_list_py.append(term_data)

    # Convert the Python list of tuples (containing NumPy arrays) to a Numba Typed List
    # Numba handles the nested NumPy arrays within the tuples correctly here.
    typed_preprocessed_list = List()
    for term in preprocessed_list_py:
        typed_preprocessed_list.append(term)

    return typed_preprocessed_list





# -------- Numba Kernels --------

# --- Numba kernel for the SERIAL approach ---
# It applies a homogeneous list of terms for a single interaction type.
@numba.njit(cache=True, parallel=False)
def _apply_terms_numba_serial(state: np.ndarray,
                       terms_list: List[PreprocessedInteractionGeneric],  # Use generic signature, Numba specializes implicitly ?
                       output_state: np.ndarray):
    """
    Numba kernel to apply a list of interaction terms of the same *interaction type*.
    Adds results into output_state (modified in-place).
    """
    dim = len(state)
    num_terms_in_list = len(terms_list)

    for i in range(num_terms_in_list):

        sites_tuple, matrix_indices_array, matrix_values_array = terms_list[i]
        num_matrix_elements = len(matrix_values_array)

        # --- Non-parallel loop over basis states ---
        for basis_index in range(dim):

            if state[basis_index] == 0: continue
            local_state = extract_local_state(basis_index, sites_tuple)

            for j in range(num_matrix_elements):

                input_local_state = matrix_indices_array[j, 0] 
                output_local_state = matrix_indices_array[j, 1] 
                matrix_value = matrix_values_array[j]
                
                if local_state == input_local_state:
                    new_basis_index = modify_global_state(basis_index, sites_tuple, output_local_state)
                    # Calculate contribution
                    contribution = matrix_value * state[basis_index]
                    output_state[new_basis_index] += contribution
                    

# --- Numba kernel for the PARALLEL approach ---
# It applies all interaction terms in a flat list of preprocessed interaction terms, using padded sites.
@numba.njit(cache=True, parallel=True)
def _apply_hamiltonian_numba_parallel(
    state: np.ndarray, # Input state w
    preprocessed_flat_list_csr: List, # Typed List of CSR term tuples
    output_state: np.ndarray # Output state u
    ):
    """
    Numba kernel implementing gather using flat list with CSR off-diagonal lookups.
    Parallelized over output basis index j.
    """
    dim = len(state)
    num_total_terms = len(preprocessed_flat_list_csr)

    # Parallel loop over OUTPUT indices j
    for j in numba.prange(dim):
        total_contribution_j = np.complex128(0.0)
        j_raw = np.int64(j) # Ensure j is usable as index

        # --- Loop over ALL terms ---
        for int_term_idx in range(num_total_terms):
            int_term_data = preprocessed_flat_list_csr[int_term_idx]
            
            (sites_padded_array, 
             int_type, 
             offdiag_offsets, 
             offdiag_target_locals, 
             offdiag_values,
             diagonal_values ) = int_term_data

            # --- Calculate local state for output j ---
            local_state_j = extract_local_state_padded(j_raw, sites_padded_array, int_type)

            # --- 1. Add diagonal contribution ---
            # u[j] = M_jj * w[j] contribution
            diag_val = diagonal_values[local_state_j]
            if diag_val != 0j: # Check if diagonal element exists (OPTIONAL?)
                total_contribution_j += diag_val * state[j_raw] 

            # --- 2. Add off-diagonal contributions ---
            # u[j] = M_jk * w[k] contributions
            # Find range in flattened arrays for source_local == local_state_j
            start = offdiag_offsets[local_state_j]
            end = offdiag_offsets[local_state_j + 1]
            # Loop over off-diagonal elements where local_state_j is the SOURCE
            for idx in range(start, end):
                local_state_k = offdiag_target_locals[idx] # Target local state k
                offdiag_val = offdiag_values[idx]     # M_jk value (local)
                # Construct the input global state k that is mapped to j via this interaction term
                k = modify_global_state_padded(j_raw, sites_padded_array, int_type, local_state_k)
                # Add M_jk * w[k] contribution
                if 0 <= k < dim:
                    total_contribution_j += offdiag_val * state[k]

        # Assign the FINAL accumulated value for output_state[j] (Numba-SAFE)
        output_state[j] = total_contribution_j




# -------- apply_hamiltonian Wrapper --------

# Define type hint for the dictionary structure
PreprocessedDictType = Dict[int, List[PreprocessedInteractionGeneric]]

def apply_hamiltonian_serial(state: np.ndarray,
                      preprocessed_dict: PreprocessedDictType):
    """
    Thin wrapper to apply the Hamiltonian using the pre-calculated dictionary
    mapping arity to lists of Numba-ready terms.
    """
    # Ensure input state vector is complex128
    if state.dtype != np.complex128:
       state = state.astype(np.complex128)

    # Initialize the output state vector
    new_state = np.zeros_like(state)

    # Iterate through the dictionary of preprocessed terms (grouped by arity)
    for interaction_type, interaction_term_list in preprocessed_dict.items():
        if interaction_term_list: # Only call if there are terms of this interaction type
            # Call the Numba function, passing the list for this specific intercation type
            # _apply_terms_numba(state, interaction_term_list, new_state)
            _apply_terms_numba_serial(state, interaction_term_list, new_state)

    return new_state


def apply_hamiltonian_parallel(state: np.ndarray,
                                  preprocessed_flat_list: List):
    """ Thin wrapper for the flat_gather kernel. """
    if state.dtype != np.complex128:
       state = state.astype(np.complex128)
    output_state = np.zeros_like(state)
    # Call the kernel ONCE
    _apply_hamiltonian_numba_parallel(state, preprocessed_flat_list, output_state)
    return output_state






# -------- Store Hamiltonian matrix in sparse format approach --------

@numba.njit(parallel=True)
def count_nnz_per_row(hilb_size, interaction_terms):
    nnz_per_row = np.zeros(hilb_size, dtype=np.int32)
    for j in prange(hilb_size):
        count = 0
        for term in interaction_terms:
            sites, arity, rowptr, colinds, values, diag = term
            local_state = extract_local_state_padded(j, sites, arity)
            # Diagonal contribution
            if diag[local_state] != 0.0:
                count += 1
            # Off-diagonal
            start = rowptr[local_state]
            end = rowptr[local_state + 1]
            count += (end - start)
        nnz_per_row[j] = count
    return nnz_per_row

@numba.njit()
def prefix_sum(arr):
    n = arr.shape[0]
    result = np.zeros(n + 1, dtype=np.int32)
    for i in range(n):
        result[i + 1] = result[i] + arr[i]
    return result

@numba.njit(parallel=True)
def fill_sparse_data(hilb_size, interaction_terms, rowptr, colind, data):
    for j in prange(hilb_size):
        index = rowptr[j]  # Start index in colind/data for row j
        for term in interaction_terms:
            sites, arity, offptrs, colinds_local, vals_local, diag = term
            local_state = extract_local_state_padded(j, sites, arity)

            # Diagonal
            val = diag[local_state]
            if val != 0.0:
                colind[index] = j
                data[index] = val
                index += 1

            # Off-diagonal
            start = offptrs[local_state]
            end = offptrs[local_state + 1]
            for k in range(start, end):
                local_in = colinds_local[k]
                global_in = modify_global_state_padded(j, sites, arity, local_in)
                colind[index] = global_in
                data[index] = vals_local[k]
                index += 1

# Wrapper function

def build_sparse_hamiltonian(hilb_size: int, interaction_terms: List):
    """
    Constructs a CSR-format sparse Hamiltonian using Numba-parallelized assembly.

    Parameters:
        hilb_size : int
            Hilbert space dimension (2^N).
        interaction_terms : numba.typed.List
            Preprocessed flat interaction terms (same structure used in on-the-fly parallel).

    Returns:
        scipy.sparse.csr_matrix
            The sparse Hamiltonian.
    """
    # Phase 1: Count NNZ per row
    nnz_per_row = count_nnz_per_row(hilb_size, interaction_terms)

    # Phase 2: Build row pointer
    rowptr = prefix_sum(nnz_per_row)
    total_nnz = rowptr[-1]

    # Phase 3: Allocate and fill
    colind = np.empty(total_nnz, dtype=np.int32)
    data = np.empty(total_nnz, dtype=np.complex128)
    fill_sparse_data(hilb_size, interaction_terms, rowptr, colind, data)

    # Wrap in scipy CSR
    return csr_matrix((data, colind, rowptr), shape=(hilb_size, hilb_size))





















# --- Define Explicit Type Hints for Preprocessed Interactions by Arity ---
# For 1-body interactions (e.g., hx, hy, hz acting on a single site)
PreprocessedInteraction1Body = Tuple[
    Tuple[int],              # sites_tuple: Exactly one site index
    np.complex128,           # factor: The prefactor (always complex)
    np.ndarray,              # matrix_elements_array: Nx3 complex array [(in, out), val]
    np.ndarray               # matrix_values_array (complex128, shape=(N,))
]

# For 2-body interactions (e.g., Jzz * Sz_i * Sz_j)
PreprocessedInteraction2Body = Tuple[
    Tuple[int, int],         # sites_tuple: Exactly two site indices
    np.complex128,           # factor: The prefactor (always complex)
    np.ndarray,              # matrix_elements_array: Nx3 complex array [(in, out), val]
    np.ndarray               # matrix_values_array (complex128, shape=(N,))
]

# For 3-body interactions
PreprocessedInteraction3Body = Tuple[
    Tuple[int, int, int],    # sites_tuple: Exactly three site indices
    np.complex128,           # factor: The prefactor (always complex)
    np.ndarray               # matrix_elements_array: Nx3 complex array [(in, out), val]
]

# For 4-body interactions (like your plaquette 'PB' or star 'PA' terms)
PreprocessedInteraction4Body = Tuple[
    Tuple[int, int, int, int],# sites_tuple: Exactly four site indices
    np.complex128,           # factor: The prefactor (always complex)
    np.ndarray,              # matrix_elements_array: Nx3 complex array [(in, out), val]
    np.ndarray               # matrix_values_array (complex128, shape=(N,))
]


# --- Define Type Hint for the FLAT list element to use in the gather approach ---
PreprocessedTermFlat = Tuple[
    np.ndarray,        # sites_padded_array (int64, shape=(MAX_ARITY,))
    int,               # actual_arity
    np.complex128,     # factor
    # Using separate index/value arrays for matrix data:
    np.ndarray,        # matrix_indices_array (int64, shape=(N_elements, 2))
    np.ndarray         # matrix_values_array (complex128, shape=(N_elements,))
]