import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

def interpolate(np.ndarray[DTYPE_t, ndim=1] rs, np.ndarray[DTYPE_t, ndim=1] zs,\
                np.ndarray[DTYPE_t, ndim=2] pot, np.ndarray[DTYPE_t, ndim=1] xgrid,\
                np.ndarray[DTYPE_t, ndim=1] ygrid, np.ndarray[DTYPE_t, ndim=1] zgrid):
   assert rs.dtype == DTYPE and rs.dtype == DTYPE and pot.dtype == DTYPE and\
      xgrid.dtype == DTYPE and ygrid.dtype == DTYPE and zgrid.dtype == DTYPE
   cdef double dr, dz, zoffset, r, normR, normRRev, z, normZ, normZRev
   cdef int ix, iy, iz, lowRInd, highRInd, lowZInd, highZInd
   cdef int xgrid_len = len(xgrid)
   cdef int ygrid_len = len(ygrid)
   cdef int zgrid_len = len(zgrid)
   cdef np.ndarray[DTYPE_t, ndim=3] pot_on_grid = np.zeros((xgrid_len,ygrid_len,zgrid_len),dtype=DTYPE)
   
   dr = rs[1]-rs[0]
   dz = zs[1]-zs[0]
   zoffset = zs[0]
   for ix in range(xgrid_len):
      for iz in range(zgrid_len):
         r = np.sqrt(xgrid[ix]**2+zgrid[iz]**2)
         lowRInd = int(r/dr)
         highRInd = lowRInd+1
         normR = (r-rs[lowRInd])/(rs[highRInd]-rs[lowRInd])
         normRRev = 1-normR
         for iy in range(ygrid_len):
            z = ygrid[iy]
            lowZInd = int((z-zoffset)/dz)
            highZInd = lowZInd+1
            normZ = (z-zs[lowZInd])/(zs[highZInd]-zs[lowZInd])
            normZRev = 1-normZ
            pot_on_grid[ix,iy,iz] = pot[lowRInd,lowZInd]*normRRev*normZRev + \
                                    pot[highRInd,lowZInd]*normR*normZRev + \
                                    pot[lowRInd,highZInd]*normRRev*normZ + \
                                    pot[highRInd,highZInd]*normR*normZ

   return pot_on_grid
