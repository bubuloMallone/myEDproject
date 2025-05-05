import numpy as np
import math

Ns = 26
rawdim = 67108864
SymmetryOps=52
print('Ns = ', Ns)
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