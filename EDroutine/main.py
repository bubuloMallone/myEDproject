import numpy as np
import sys
import line_profiler
import atexit
import os
import numba
import time

# Import the builder function
from utils.hamiltonian_utils import (
    build_numba_preprocessed_dict,
    build_numba_preprocessed_list,
    build_sparse_hamiltonian  # <-- newly added parallel builder
)
from utils.input_utils import (
    parse_lattice_file, parse_matrix_file, parse_input_string,
    parse_tags_factors_matrices, parse_movelist, write_to_file
)
from utils.lanczos import (
    lanczos_scipy_serial, lanczos_scipy_parallel,
    sort_eigenpairs, lanczos_scipy_sparse, lanczos_primme_sparse
)
from utils.observable_utils import (
    translation_op, hadamard_op, compute_duality_symmetry
)
from utils.debug_utils import (
    check_duality_symmetry_fly, check_duality_symmetry_sparse,
    print_hamiltonian_terms
)

# # --- Set up line profiler for performance analysis ---
# profile = line_profiler.LineProfiler()
# atexit.register(profile.print_stats)

MAX_INT_TYPE = 4

# @profile
def main():

    # ---- Set maximum number of available cores ----
    default_threads = 1
    max_threads = os.cpu_count()
    if max_threads is None: # cpu_count can return None
        print(f"Warning: Could not determine number of CPU cores. Defaulting to {default_threads} threads.")
        max_threads = default_threads
    # --- End maximum number of available cores ----
        
    # ---- Read input from stdin ----
    input_text = sys.stdin.read()
    params = parse_input_string(input_text)
    print("Parsed input parameters")
    # ---- End stdin input parsing ----

    # ---- Extract input parameters ----
    threads = params["threads"]
    if threads > max_threads:
        print(f"Warning: Requested {threads} threads exceeds available {max_threads} cores. Defaulting to {max_threads} threads.")
        threads = max_threads
    elif threads < 1:
        print(f"Warning: Requested {threads} threads is less than 1. Defaulting to {default_threads} threads.")
        threads = default_threads
    sparse_method = params["sparseMethod"]
    lattice_file = params["LatticeDescriptionFile"]
    output_file = params["OutFile"]
    num_eigenvalues = params["Nevals"]
    if "dumpgs" in params:
        dumpgs = params["dumpgs"]
        gs_file = params["dumpgsFile"]
    else: dumpgs = False
    print("Lattice file:", lattice_file)
    print("Output file:", output_file)
    print("Number of eigenvalues to compute:", num_eigenvalues)
    valid_tags, matrices = parse_tags_factors_matrices(params)
    print("Parsed interaction types and corresponding matrix files")
    sites, interactions = parse_lattice_file(lattice_file, valid_tags)
    n_sites = len(sites)
    hilbert_size = 2 ** n_sites
    print("Hilbert space dim:", hilbert_size)
    print("Parsed lattice interactions")
    # ---- End input parameters parsing ----

    # --- Set Numba Threads ---
    print(f"Requested {threads} threads for Numba parallelization")
    try:
        numba.set_num_threads(threads)
        print(f"Numba configured to use {numba.get_num_threads()} threads.") # Verify setting
    except ValueError as e:
         print(f"Error setting Numba threads: {e}. Using Numba default.")
    # --- End Numba Threads ---

    # --- Report OMP / BLAS threading setup ---
    print("OPENBLAS_NUM_THREADS =", os.environ.get("OPENBLAS_NUM_THREADS", "not set"))
    # OpenMP (used by Numba, PRIMME, etc.)
    print("OMP_NUM_THREADS =", os.environ.get("OMP_NUM_THREADS", "not set"))
    # --- End NumPy / BLAS threading setup ---

    if sparse_method:
        # ---- Perform Numba preparation of data ----
        preprocessed_hamiltonian_list = build_numba_preprocessed_list(interactions, matrices, valid_tags, MAX_INT_TYPE)
        num_processed_terms = len(preprocessed_hamiltonian_list)
        print(f"Built Numba-friendly interaction list containing {num_processed_terms} terms.")
        # print_hamiltonian_terms(preprocessed_hamiltonian_dict)
        # --- End Numba preparation ---

        print("Constructing full sparse Hamiltonian")
        build_sparse_hamiltonian_start_time = time.time()
        H_sparse = build_sparse_hamiltonian(hilbert_size, interaction_terms=preprocessed_hamiltonian_list)
        build_sparse_hamiltonian_end_time = time.time()
        print("Sparse Hamiltonian construction finished")
        print(f"Sparse Hamiltonian construction time: {build_sparse_hamiltonian_end_time - build_sparse_hamiltonian_start_time:.3f} seconds")
        print("Entering Lanczos computation")
        lanczos_start_time = time.time()
        # ---- Perform Lanczos computation ----
        eigenvalues, eigenvectors = lanczos_primme_sparse(H_sparse, hilbert_size, num_eigenvalues, tolerance=1e-12, maxiter=100000)
        lanczos_end_time = time.time()
        print("Lanczos computation finished")
        print(f"Lanczos computation time: {lanczos_end_time - lanczos_start_time:.3f} seconds")

    else:
        if threads == 1:
            # ---- Perform Numba preparation of data ----
            preprocessed_hamiltonian_dict = build_numba_preprocessed_dict(interactions, matrices, valid_tags)
            print(f"Built Numba-friendly interaction dictionary (keys are interaction types: {list(preprocessed_hamiltonian_dict.keys())})")
            # print_hamiltonian_terms(preprocessed_hamiltonian_dict)
            # --- End Numba preparation ---

            # ---- Lanczos computation ----
            print("Entering Lanczos computation")
            eigenvalues, eigenvectors = lanczos_scipy_serial(hilbert_size=hilbert_size,
                                                            preprocessed_hamiltonian_dict=preprocessed_hamiltonian_dict,
                                                            num_eigenvalues=num_eigenvalues )
            print("Lanczos computation finished")

        else:  # Parallel computation
            # ---- Perform Numba preparation of data ----
            preprocessed_hamiltonian_list = build_numba_preprocessed_list(interactions, matrices, valid_tags, MAX_INT_TYPE)
            num_processed_terms = len(preprocessed_hamiltonian_list)
            print(f"Built Numba-friendly interaction list containing {num_processed_terms} terms.")
            # print_hamiltonian_terms(preprocessed_hamiltonian_dict)
            # --- End Numba preparation ---

            # ---- Lanczos computation ----
            print("Entering Lanczos computation")
            eigenvalues, eigenvectors = lanczos_scipy_parallel( hilbert_size=hilbert_size,
                                                                preprocessed_hamiltonian_list=preprocessed_hamiltonian_list,
                                                                num_eigenvalues=num_eigenvalues )
            print("Lanczos computation finished")
            # --- End Lanczos computation ---


    eigenvalues, eigenvectors = sort_eigenpairs(eigenvalues, eigenvectors)

    # ---- Compute duality symmetry eigenvalues ----
    movelist = parse_movelist(n_sites)
    if sparse_method:
        print("Checking duality symmetry using sparse Hamiltonian")
        check_duality_symmetry_sparse(H_sparse, sites, movelist)
    else:
        print("Checking duality symmetry using Hamiltonian action")
        check_duality_symmetry_fly(preprocessed_hamiltonian_dict, sites, movelist)

    tolerance = 1e-10
    duality_expectations = compute_duality_symmetry(eigenvalues, eigenvectors, n_sites, movelist, tolerance)

    # ---- Write eigenvalues to file ----
    write_to_file(output_file, 'w', eigenvalues, "Eigenvalue")
    # ---- Write duality quantum numbers to file ----
    write_to_file(output_file, 'a', duality_expectations, "Duality")
    # ---- Write ground state to file ---- (optional)
    if dumpgs: write_to_file(gs_file, 'w', eigenvectors[:, 0], "ground_state")


if __name__ == "__main__":
    main()