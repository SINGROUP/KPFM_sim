# -*- coding: utf-8 -*-

import pyximport; pyximport.install()
import sys, math, time
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d

import cube_file_writer as cube_writer
import axisym_to_3d_grid_bilinear_interpolation as rs_grid_interp
import grid_linear_interpolation_extparam as s_interp
from cp2k_grid_size import get_grid_size
from pot_rounding import pot_rounding

eV_to_au = 0.03674932540261308
rydberg_to_au = 0.5

# Gap between the lower atoms of the sample and the surface of the macroscopic sample                                                                                                                            
# and between the upper atoms of the tip and the macroscopic tip                                                                                                                                                 
macro_atomic_gap = 6.0 # Å

# Standard deviation for the Gaussian blur
gaussian_blur_stdev = 0.5 # Å


# Function for loading the electrostatic potential from database first
# interpolating with respect to available macroscopic tip-sample distance
# values and then interpolating it to the 3d realspace grid of CP2k and
# writing it to a cube file.
def axisym_pot_in_db_to_cube(out_file_name, s, atoms, result_db):
    print("Fetching external potential data from database and interpolating " +\
            "with respect to macroscopic tip-sample distance s")
    s_macro = s + 2*macro_atomic_gap
    lower_s_point, higher_s_point = result_db.get_closest_pot_scan_points(s_macro)
    if (lower_s_point is None) and (higher_s_point is None):
        raise Exception("External potential scan points were not found")
    elif lower_s_point is None:
        lower_s_point = higher_s_point
    elif higher_s_point is None:
        higher_s_point = lower_s_point
        
    id_lower = lower_s_point[0]
    s_lower = lower_s_point[1]
    id_higher = higher_s_point[0]
    s_higher = higher_s_point[1]
    
    print("s = {}, s_lower = {}, s_higher = {}\n".format(s_macro, s_lower, s_higher))
    
    rs, zs, V_lower = result_db.get_external_potential(id_lower)
    rs, zs, V_higher = result_db.get_external_potential(id_higher)
    
    print("V_lower = ")
    print(V_lower)
    print("V_higher = ")
    print(V_higher)
    
    pot = s_interp.interpolate_2d(V_lower, V_higher, s_macro, s_lower, s_higher)
    
    print("Potential interpolated with respect to s:")
    print(pot)
    print("")
    
    cell = atoms.get_cell()
    cell_dim = [cell[0,0], cell[1,1], cell[2,2]]
    cutoff = float(result_db.get_sim_parameter("cutoff"))
    vacuum = float(result_db.get_sim_parameter("vacuum"))
    if cutoff is None:
        raise Exception("Cutoff parameter was not found in the results database.")
    if vacuum is None:
        raise Exception("Vacuum parameter was not found in the results database.")
    xgrid, ygrid, zgrid = init_pot_grid_(cell_dim, vacuum, macro_atomic_gap, cutoff)
    
    print("Interpolating axisymmetric potential to 3d grid...")
    pot_on_grid = rs_grid_interp.interpolate(rs, zs, pot, xgrid, ygrid, zgrid)
    print('Done\n')
    
    print("Gaussian blurring the potential to get rid of the discontinuities of the first derivative...")
    dx = xgrid[1] - xgrid[0]
    dy = ygrid[1] - ygrid[0]
    dz = zgrid[1] - zgrid[0]
    sigma_x = gaussian_blur_stdev / dx
    sigma_y = gaussian_blur_stdev / dy
    sigma_z = gaussian_blur_stdev / dz
    pot_on_grid = gaussian_filter1d(pot_on_grid, sigma_x, axis=0, mode='wrap', truncate=6.0)
    pot_on_grid = gaussian_filter1d(pot_on_grid, sigma_y, axis=1, mode='nearest', truncate=6.0)
    pot_on_grid = gaussian_filter1d(pot_on_grid, sigma_z, axis=2, mode='wrap', truncate=6.0)
    print('Done\n')
    
    '''
    round_dist = (len(xgrid)+len(zgrid))/100
    print 'Rounding the edges of the potential... (rounding distance = ' + \
        repr(round_dist) + ')'
    pot_on_grid = pot_rounding(xgrid, ygrid, zgrid, pot_on_grid, round_dist)
    print 'Done\n'
    '''

    print('Writing interpolated electrostatic potential on grid to file...')
    start_time = time.clock()
    pot_on_grid = pot_on_grid*eV_to_au
    cube_writer.write_cube_orthogonal(out_file_name, xgrid, ygrid, zgrid, pot_on_grid,
                                    comment_line='CP2k external electrostatic potential')
    end_time = time.clock()
    print('Done in {} seconds\n'.format(end_time-start_time))


def init_pot_grid_(cell_dim, vacuum, macro_atomic_gap, cutoff):
    n_grid = get_grid_size(cell_dim, cutoff)
    xgrid = np.linspace(-0.5*cell_dim[0], 0.5*cell_dim[0], n[0],
                        endpoint=False)
    ygrid = np.linspace(macro_atomic_gap-vacuum, cell_dim[1]+macro_atomic_gap-vacuum,
                        n[1], endpoint=False)
    zgrid = np.linspace(-0.5*cell_dim[2], 0.5*cell_dim[2], n[2],
                        endpoint=False)
    return xgrid, ygrid, zgrid
