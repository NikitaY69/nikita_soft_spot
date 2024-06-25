import numpy as np
import os
import scipy.sparse as sp
from scipy.sparse.linalg import eigs

import soft_spot_project.sim_lib as sim



def participation_ratio(z):
    if len(z.shape) == 1:
        z = np.reshape(z, (sim.system.npart, -1))

    s1 = 0.
    s2 = 0.
    for i in range(z.shape[0]):
        s1 += np.dot(z[i,:],z[i,:])
        s2 += np.power(np.dot(z[i,:],z[i,:]),2.)
    return s1*s1/(s2*z.shape[0])




def get_sparse_hessian():
    size_hess = sim.system.ndim*sim.system.ndim*sim.system.npart + sim.system.ndim*sim.system.ndim*sim.system.ndim*sim.pairwise.pair_count
    hess = np.zeros(shape=(3,size_hess),order='F')
    sim.pairwise.compute_hessian(hess)
    hessian = sp.csr_matrix((hess[2,:], (hess[0,:].astype(int), hess[1,:].astype(int))), dtype=np.float64)
    return hessian



