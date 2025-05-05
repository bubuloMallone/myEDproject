import numpy as np
from scipy.sparse.linalg import eigsh
from hamiltonian import get_hamiltonian_operator

def lanczos_ed(N, h, k=5):
    """Computes the lowest k eigenvalues of the TFIM Hamiltonian using Lanczos."""
    H = get_hamiltonian_operator(N, h)
    eigvals, eigvecs = eigsh(H, k=k, which='SA')  # 'SA' → Smallest Algebraic eigenvalues
    return eigvals, eigvecs