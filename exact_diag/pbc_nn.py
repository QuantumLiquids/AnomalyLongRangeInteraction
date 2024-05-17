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
from quspin.operators import hamiltonian, quantum_operator, quantum_LinearOperator  # Hamiltonians and operators
from quspin.basis import spin_basis_1d  # Hilbert space spin basis
from quspin.tools.evolution import expm_multiply_parallel

from scipy.sparse.linalg import eigsh

##### define model parameters #####
L = 24  # system size
Jzz = 1
theta = 0.6
omega_0 = math.cos(theta)
omega_1 = math.sin(theta)

os.environ['OMP_NUM_THREADS'] = '4'  # set number of OpenMP threads to run in parallel
os.environ['MKL_NUM_THREADS'] = '4'  # set number of MKL threads to run in parallel
#

s = np.arange(L)  # sites [0,1,2,....]
Z = -(s + 1)  # spin inversion
P = s[::-1]

E_list = []
UX_list = []

for k_int in range(L):
    # basis = spin_basis_general(L,pauli=0,Nup=L//2,zblock=(Z,0),pblock=(P,0)) # and positive parity sector
    basis = spin_basis_1d(L, pauli=1, kblock=k_int)
    dim = basis.Ns
    print("basis dimension : {}".format(dim))
    print("basis blocks : {}".format(basis.blocks))
    print("basis info : {}".format(basis.description))

    # define operators with OBC using site-coupling lists
    # J_zz = [[Jzz / 2.0 / min((j - i) % L, (i - j) % L) ** alpha, i, (i+2)%L] for i in range(L)]
    Ja = [[0.5 * omega_0, i, (i + 1) % L] for i in range(L)]
    Jb = [[0.5 * omega_1, i] for i in range(L)]
    Jc = [[-0.5 * omega_1, i, (i + 1) % L, (i + 2) % L] for i in range(L)]

    # static and dynamic lists
    static = [["zx", Ja], ["xz", Ja], ["x", Jb], ["zxz", Jc]]
    dynamic = []
    # compute the time-dependent Heisenberg Hamiltonian
    H_XXZ = hamiltonian(static, dynamic, basis=basis, dtype=np.complex128)
    #
    ##### various exact diagonalisation routines #####
    # E,V=eigsh(H_XXZ.aslinearoperator())
    E, V = H_XXZ.eigsh(k=min(10, dim), which='SA', return_eigenvectors=True)
    print(E)
    E_list.append(E.tolist())

    # x_prod_coup_sites = [[1.0] + range(L)]
    # x_prod_dict={"spin_flip_operator": [["x"*L, x_prod_coup_sites]]}
    # x_prod_opt = quantum_operator(x_prod_dict, N = L, basis = basis)

    x_prod_coup_sites = [[1.0] + list(range(L))]
    x_prod_static=[["x"*L, x_prod_coup_sites]]
    x_prod_opt = quantum_LinearOperator(x_prod_static, basis=basis, dtype=np.complex128)

    ising_zz_coup_sites = [[1.0, i, (i + 1) % L] for i in range(L)]
    ising_zz_static = [["zz",ising_zz_coup_sites ]]
    ising_zz_opt = hamiltonian(ising_zz_static, dynamic, basis=basis, dtype=np.complex128)
    
    uxs = []
    for s in range(10):
        v = V[:, s]
        v1 = x_prod_opt.dot(v)
        v2 = expm_multiply_parallel(ising_zz_opt.tocsr(), a = -1.0j*math.pi/4.0).dot(v1)
        ux = np.conjugate(v).dot(v2) * np.exp(1.0j*math.pi/4.0 * L)
        uxs.append(int(round(ux)))

    UX_list.append(uxs)

    
theta_str = "{:.4f}".format(theta)
filename = f"EnergyLocHamN{L}theta{theta_str}.txt"
np.savetxt(filename, E_list, fmt="%s")

filename = f"UXLocHamN{L}theta{theta_str}.txt"
np.savetxt(filename, UX_list, fmt="%s")
