from __future__ import print_function, division
import numpy as np  # generic math functions

np.object = object

import os
import math
#####################################################################
#                            example 0                              #
#    In this script we demonstrate how to use QuSpin's exact        #
#    diagonlization routines to solve for the eigenstates and       #
#    energies of the XXZ chain.                                     #
#####################################################################
from quspin.operators import hamiltonian  # Hamiltonians and operators
from quspin.basis import spin_basis_1d  # Hilbert space spin basis

from scipy.sparse.linalg import eigsh

##### define model parameters #####
L = 20  # system size
Jzz = 1
theta = math.pi / 4
alpha = 1
omega_0 = math.cos(theta)
omega_1 = math.sin(theta)

os.environ['OMP_NUM_THREADS'] = '4'  # set number of OpenMP threads to run in parallel
os.environ['MKL_NUM_THREADS'] = '4'  # set number of MKL threads to run in parallel
#

s = np.arange(L)  # sites [0,1,2,....]
Z = -(s + 1)  # spin inversion
P = s[::-1]

E_list = []

for k_int in range(L):
    # basis = spin_basis_general(L,pauli=0,Nup=L//2,zblock=(Z,0),pblock=(P,0)) # and positive parity sector
    basis = spin_basis_1d(L, pauli=1, kblock=k_int)
    dim = basis.Ns
    print("basis dimension : {}".format(dim))
    print("basis blocks : {}".format(basis.blocks))
    print("basis info : {}".format(basis.description))

    # define operators with OBC using site-coupling lists
    J_zz = [[Jzz / 2.0 / min((j - i) % L, (i - j) % L) ** alpha, i, j] for i in range(L) for j in range(i + 2, L)]
    Ja = [[0.5 * omega_0, i, (i + 1) % L] for i in range(L)]
    Jb = [[0.5 * omega_1, i] for i in range(L)]
    Jc = [[-0.5 * omega_1, i, (i + 1) % L, (i + 2) % L] for i in range(L)]

    # static and dynamic lists
    static = [["zx", Ja], ["xz", Ja], ["x", Jb], ["zxz", Jc], ["zz", J_zz]]
    dynamic = []
    # compute the time-dependent Heisenberg Hamiltonian
    H_XXZ = hamiltonian(static, dynamic, basis=basis, dtype=np.complex128)
    #
    ##### various exact diagonalisation routines #####
    # E,V=eigsh(H_XXZ.aslinearoperator())
    E = H_XXZ.eigsh(k=min(8, dim), which='SA', return_eigenvectors=False)
    E.sort()
    print(E)
    E_list.append(E.tolist())

np.savetxt("data.txt", E_list, fmt="%s")
