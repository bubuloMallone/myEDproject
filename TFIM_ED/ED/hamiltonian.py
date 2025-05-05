import numpy as np
from scipy.sparse.linalg import LinearOperator

def hamiltonian_action(psi, N, h):
    """Computes H|psi⟩ for the transverse field Ising model dynamically."""
    psi_out = np.zeros_like(psi)

    # Interaction term S^z S^z
    for basis in range(len(psi)):
        if psi[basis] == 0:
            continue  # Skip zero coefficients
        
        for i in range(N - 1):
            if ((basis >> i) & 1) == ((basis >> (i + 1)) & 1):  
                psi_out[basis] -= psi[basis]  # Aligned spins (↓↓ or ↑↑) → -1
            else:
                psi_out[basis] += psi[basis]  # Anti-aligned spins (↓↑ or ↑↓) → +1

        # Transverse field term S^x
        for i in range(N):
            flipped_basis = basis ^ (1 << i)  # Flip spin i using XOR
            psi_out[flipped_basis] -= h * psi[basis]

    return psi_out

def get_hamiltonian_operator(N, h):
    """Returns a SciPy LinearOperator representing the TFIM Hamiltonian."""
    dim = 2**N
    return LinearOperator((dim, dim), matvec=lambda psi: hamiltonian_action(psi, N, h))