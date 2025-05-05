import numpy as np
from ED.lanczos import lanczos_ed

N = 6  # Number of spins
h = 1.0  # Transverse field strength
k = 5   # Number of eigenvalues to compute

eigvals, eigvecs = lanczos_ed(N, h, k)

print("Lowest eigenvalues:", eigvals)