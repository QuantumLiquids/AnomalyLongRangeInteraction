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
from quspin.operators import hamiltonian,quantum_operator  # Hamiltonians and operators
from quspin.basis import spin_basis_1d  # Hilbert space spin basis

from scipy.sparse.linalg import eigsh

##### define model parameters #####
L = 16  # system size
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
BZZ_list = [] # boundary sigma_z * sigma_z, distinguish even/odd number of domain wall


for p_int in range(2):
    # basis = spin_basis_general(L,pauli=0,Nup=L//2,zblock=(Z,0),pblock=(P,0)) # and positive parity sector
    basis = spin_basis_1d(L, pauli=1, pblock=p_int * 2 - 1)
    dim = basis.Ns
    print("basis dimension : {}".format(dim))
    print("basis blocks : {}".format(basis.blocks))
    print("basis info : {}".format(basis.description))

    # define operators with OBC using site-coupling lists
    J_zz = [[Jzz / 2.0 / (j - i) ** alpha, i, j] for i in range(L) for j in range(i + 2, L)]
    Ja1 = [[0.5 * omega_0, i, i + 1] for i in range(L - 2)]
    Ja2 = [[0.5 * omega_0, i, i + 1] for i in range(1, L - 1)]
    Jb = [[0.5 * omega_1, i] for i in range(1, L - 1)]
    Jc = [[-0.5 * omega_1, i, (i + 1), (i + 2)] for i in range(L - 2)]

    # static and dynamic lists
    static = [["zx", Ja1], ["xz", Ja2], ["x", Jb], ["zxz", Jc], ["zz", J_zz]]
    dynamic = []
    # compute the time-dependent Heisenberg Hamiltonian
    H_XXZ = hamiltonian(static, dynamic, basis=basis, dtype=np.float64)
    #
    ##### various exact diagonalisation routines #####
    # E,V=eigsh(H_XXZ.aslinearoperator())
    E, V = H_XXZ.eigsh(k=min(10, dim), which='SA', return_eigenvectors=True)
    E.sort()
    print(E)
    E_list.append(E.tolist())

    bdry_sites = [[1.0, 0, L-1]]
    bdry_zz_dict={"a": [["zz", bdry_sites]]}
    bdry_zz_opt = quantum_operator(bdry_zz_dict, N = L, basis = basis)
    bzz = []
    for s in range(10):
        v = V[:, s]
        bzz.append(int(round(bdry_zz_opt.matrix_ele(v,v,{"a":1}))))
    BZZ_list.append(bzz)

theta_str = "{:.4f}".format(theta)
filename = f"EnergyOBCN{L}theta{theta_str}alpha{alpha}.txt"
np.savetxt(filename, E_list, fmt="%s")

filename = f"BZZOBCN{L}theta{theta_str}alpha{alpha}.txt"
np.savetxt(filename, BZZ_list, fmt="%s")