import numpy as np
import math

Ns = 12
Nmax = 4
Nb = Ns

# n = math. factorial(Nmax * Ns)
# d = math.factorial(Nb)
# d2 = math.factorial(Nmax * Ns - Nb)


# n = math. factorial(Nb + Ns -1)
# d = math.factorial(Nb)
# d2 = math.factorial(Ns - 1)

# maxdim = n/d/d2
#  rawdim (Nmax=4) = 187500119445
#  rawdim (Nmax=3) = 42478745781

rawdim = 1335698
SymmetryOps=144
print('Ns = ', Ns)
print('Nmax = ', Nmax)
print('rawdim: ', rawdim)
lookupTable = rawdim * 8 / (1024)**3
normTable = 1 * rawdim / (1024)**3
LanczosVecDim = 16 * rawdim/SymmetryOps / (1024)**3
SymmetrizedState = 8 * rawdim/SymmetryOps / (1024)**3

print('lookup table (GB): ', lookupTable)
print('normtable (GB): ', normTable)
print('Lanczos vector (GB): ', LanczosVecDim)
print('Symmetrized states table (GB): ', SymmetrizedState)

memory = lookupTable + normTable + 4 * LanczosVecDim + SymmetrizedState

print('memory (GB): ', memory)


# For N=21 triangular
# Found first rawdim= 1105350729
# Maxdim= 10460353203

# For N=24 kagome
# rawdim = 27114249960