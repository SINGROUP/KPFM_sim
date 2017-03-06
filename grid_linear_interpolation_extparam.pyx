import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

def interpolate_2d(np.ndarray[DTYPE_t, ndim=2] V_lower, np.ndarray[DTYPE_t, ndim=2] V_higher,
                double s, double s_lower, double s_higher):
    assert V_lower.dtype == DTYPE and V_higher.dtype == DTYPE
    
    cdef double s_norm
    cdef int nx = V_lower.shape[0]
    cdef int ny = V_higher.shape[1]
    cdef int ix, iy
    cdef np.ndarray[DTYPE_t, ndim=2] V_s = np.zeros([nx, ny], dtype=DTYPE)
    
    s_norm = (s-s_lower)/(s_higher-s_lower)
    for ix in range(nx):
        for iy in range(ny):
            V_s[ix, iy] = (1-s_norm)*V_lower[ix, iy] + s_norm*V_higher[ix, iy]
            
    return V_s
