import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh, ArpackNoConvergence
# Import the actual function that applies the Hamiltonian
from utils.hamiltonian_utils import apply_hamiltonian_serial, apply_hamiltonian_parallel, PreprocessedDictType
from numba.typed import List

import line_profiler
import atexit

# --- Set up line profiler for performance analysis ---
profile = line_profiler.LineProfiler()
atexit.register(profile.print_stats)

# Modify signature: remove old dicts, add preprocessed_hamiltonian_dict
# @profile
def lanczos_scipy_serial(hilbert_size: int,
                  preprocessed_hamiltonian_dict: PreprocessedDictType,
                  num_eigenvalues: int = 1):
     """
     Compute the lowest eigenvalues using SciPy's Lanczos method
     with a pre-processed Hamiltonian structure.
     """

     # Define the matvec function using the preprocessed dictionary
     # Option 1: Lambda (simpler change here)
     def matvec(psi):
          # apply_hamiltonian now takes the state and the preprocessed dict
          return apply_hamiltonian_serial(psi, preprocessed_hamiltonian_dict)

     # Use the chosen matvec implementation
     H = LinearOperator((hilbert_size, hilbert_size), matvec=matvec, dtype=np.complex128)

     # ARPACK requires k < N, handle small Hilbert spaces if necessary
     k_eff = min(num_eigenvalues, hilbert_size - 1)
     if k_eff <= 0:
          raise ValueError("Number of eigenvalues 'k' must be > 0 and < Hilbert space size.")
     if k_eff != num_eigenvalues:
          print(f"Warning: Requesting {num_eigenvalues} eigenvalues, but Hilbert space size {hilbert_size} only allows k={k_eff}.")

     try:
          eigenvalues, eigenvectors = eigsh(H, k=k_eff, which='SA', tol=1e-15, maxiter=100000)
     except ArpackNoConvergence as e:
          print("WARNING: ARPACK did not converge!")
          eigenvalues = e.eigenvalues
          eigenvectors = e.eigenvectors
     
     return eigenvalues, eigenvectors


def lanczos_scipy_parallel(hilbert_size: int,
                  preprocessed_hamiltonian_list: List,
                  num_eigenvalues: int = 1):
    """
    Compute the lowest eigenvalues using SciPy's Lanczos method
    with a pre-processed Hamiltonian structure.
    """

    # Define the matvec function using the preprocessed dictionary
    # Option 1: Lambda (simpler change here)
    def matvec(psi):
        # apply_hamiltonian now takes the state and the preprocessed dict
        return apply_hamiltonian_parallel(psi, preprocessed_hamiltonian_list)

    # Option 2: Class (more verbose but potentially cleaner)
    # class HamiltonianOperator(LinearOperator):
    #     def __init__(self, dtype, shape, hamiltonian_dict):
    #         super().__init__(dtype, shape)
    #         self.hamiltonian_dict = hamiltonian_dict
    #
    #     def _matvec(self, psi):
    #         return apply_hamiltonian(psi, self.hamiltonian_dict)
    # H = HamiltonianOperator(dtype=np.complex128, shape=(hilbert_size, hilbert_size), hamiltonian_dict=...)

    # Use the chosen matvec implementation
    H = LinearOperator((hilbert_size, hilbert_size), matvec=matvec, dtype=np.complex128)

    # ARPACK requires k < N, handle small Hilbert spaces if necessary
    k_eff = min(num_eigenvalues, hilbert_size - 1)
    if k_eff <= 0:
         raise ValueError("Number of eigenvalues 'k' must be > 0 and < Hilbert space size.")
    if k_eff != num_eigenvalues:
         print(f"Warning: Requesting {num_eigenvalues} eigenvalues, but Hilbert space size {hilbert_size} only allows k={k_eff}.")


    eigenvalues, eigenvectors = eigsh(H, k=k_eff, which='SA')

    return eigenvalues, eigenvectors



def sort_eigenpairs(eigenvalues, eigenvectors):
    """
    Sort eigenvalues and corresponding eigenvectors in ascending order.
    """

     # Sort eigenvalues and reorder eigenvectors accordingly
    idx_sorted = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[idx_sorted]
    eigenvectors = eigenvectors[:, idx_sorted]

    return eigenvalues, eigenvectors



def lanczos_scipy_sparse(matrix_sparse,
                         hilbert_size: int,
                         num_eigenvalues: int = 1,
                         tolerance: float = 1e-14,
                         maxiter: int = 100000):
     """
     Compute the lowest eigenvalues using SciPy's Lanczos method
     with a pre-processed Hamiltonian sparse matrix.
     """

     # ARPACK requires num_eigenvalues < hilbert_size
     k_eff = min(num_eigenvalues, hilbert_size - 1)
     if k_eff <= 0:
          raise ValueError("Number of eigenvalues 'k' must be > 0 and < Hilbert space size.")
     if k_eff != num_eigenvalues:
          print(f"Warning: Requesting {num_eigenvalues} eigenvalues, but Hilbert space size {hilbert_size} only allows k={k_eff}.")

     try:
          eigenvalues, eigenvectors = eigsh(matrix_sparse, k=k_eff, which='SA', tol=tolerance, maxiter=maxiter, ncv= (4 * k_eff + 1))
     except ArpackNoConvergence as e:
          print("WARNING: ARPACK did not converge!")
          eigenvalues = e.eigenvalues
          eigenvectors = e.eigenvectors
     
     return eigenvalues, eigenvectors