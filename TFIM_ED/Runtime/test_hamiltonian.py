import numpy as np
from ED.hamiltonian import hamiltonian_action

N = 4
h = 1.0
dim = 2**N

# Test on a random quantum state
psi = np.random.rand(dim)
Hpsi = hamiltonian_action(psi, N, h)

print("H|psi⟩ computed successfully!")