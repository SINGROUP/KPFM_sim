# -*- coding: utf-8 -*-

def pot_rounding(xgrid, ygrid, zgrid, V, round_dist):
  nx = len(xgrid)
  ny = len(ygrid)
  nz = len(zgrid)
  dx = xgrid[1]-xgrid[0]
  dz = zgrid[1]-zgrid[0]
  xtrans = xgrid-xgrid[0]
  ztrans = zgrid-zgrid[0]
  
  for iy in range(ny):
    # round the edges along x-axis
    for ix in range(round_dist+1,nx-round_dist):
      dV = 0.5*(V[ix,iy,round_dist+1]-V[ix,iy,round_dist-1])/dz;
      a = 0.5*dV/ztrans[round_dist]
      b = V[ix,iy,round_dist]-0.5*ztrans[round_dist]*dV
      V[ix,iy,0] = b
      for iz in range(1,round_dist):
        V[ix,iy,iz] = a*ztrans[iz]**2+b
        V[ix,iy,nz-iz] = V[ix,iy,iz]
        
    # round the edges along z-axis
    for iz in range(round_dist+1,nz-round_dist):
      dV = 0.5*(V[round_dist+1,iy,iz]-V[round_dist-1,iy,iz])/dx;
      a = 0.5*dV/xtrans[round_dist]
      b = V[round_dist,iy,iz]-0.5*xtrans[round_dist]*dV
      V[0,iy,iz] = b
      for ix in range(1,round_dist):
        V[ix,iy,iz] = a*xtrans[ix]**2+b
        V[nx-ix,iy,iz] = V[ix,iy,iz]
        
    # round the corner areas
    dxV = 0.5*(V[round_dist+1,iy,round_dist]-V[round_dist-1,iy,round_dist])/dx
    dzV = 0.5*(V[round_dist,iy,round_dist+1]-V[round_dist,iy,round_dist-1])/dz
    a = 0.5*dxV/xtrans[round_dist]
    b = 0.5*dzV/ztrans[round_dist]
    c = V[round_dist,iy,round_dist]-0.5*xtrans[round_dist]*dV-0.5*ztrans[round_dist]*dV
    V[0,iy,0] = c
    for k in range(1,round_dist+1):
      V[k,iy,0] = a*xtrans[k]**2+c
      V[nx-k,iy,0] = V[k,iy,0]
      V[0,iy,k] = b*ztrans[k]**2+c
      V[0,iy,nz-k] = V[0,iy,k]
    for ix in range(1,round_dist+1):
      for iz in range(1,round_dist+1):
        V[ix,iy,iz] = a*xtrans[ix]**2+b*ztrans[iz]**2+c
        V[ix,iy,nz-iz] = V[ix,iy,iz]
        V[nx-ix,iy,iz] = V[ix,iy,iz]
        V[nx-ix,iy,nz-iz] = V[ix,iy,iz]
  
  return V
