cimport numpy as np
import numpy as np
cimport cython

from libc.math cimport sqrt

@cython.boundscheck(False)
@cython.profile(False)
@cython.wraparound(False)
cpdef np.ndarray[np.int32_t, ndim=1] norm_coods(np.ndarray g, np.ndarray c, float R, float half_deltac, double L=500):

    assert (g.dtype == np.float and c.dtype == np.float)

    g -= c
    g[g > L/2] -= L
    g[g < -L/2] += L

    cdef np.ndarray z_filt
    z_filt = g[:,2] < 0

    g[z_filt, 2] = -g[z_filt, 2]

    return ((g[:,0]**2 + g[:,1]**2)**0.5 <= R) & ((g[:,2]) <= half_deltac)
