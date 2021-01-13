# -*- coding: utf-8 -*-

import pyximport; pyximport.install()
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d

from cp2k_grid_size import get_grid_size
import cube_file_writer as cube_writer

eV_to_au = 0.03674932540261308

# Standard deviation for the Gaussian blur
gaussian_blur_stdev = 0.2 # Å

def create_linear_pot_cube(cube_file, linear_dim_upper_bottom, linear_dim_lower_top, V_bias, cell_vectors,
                            cutoff, linear_dim, extended_fft_lengths=False, pbc=False):
    # Initialize CP2k grid
    print('- Initializing grid:')
    n_grid, grid_vectors, pot_on_grid = init_grid_(cell_vectors, cutoff, extended_fft_lengths)
    
    # Calculate the values of the potential at each point (assuming that the cell vector in
    # linear potential dimension is parallel to one of the Cartesian axes)
    print('- Calculating potential values on grid')
    for i_linear_dim in range(n_grid[linear_dim]):
        pos_linear_dim = i_linear_dim*grid_vectors[linear_dim, linear_dim]
        
        if pos_linear_dim >= linear_dim_upper_bottom:
            potential_value = V_bias
        elif pos_linear_dim <= linear_dim_lower_top:
            potential_value = 0
        else:
            potential_value = (pos_linear_dim-linear_dim_lower_top)/(linear_dim_upper_bottom-linear_dim_lower_top) * V_bias
        
        if linear_dim == 0:
            pot_on_grid[i_linear_dim, :, :] = potential_value
        elif linear_dim == 1:
            pot_on_grid[:, i_linear_dim, :] = potential_value
        else:
            pot_on_grid[:, :, i_linear_dim] = potential_value
    
    #for ix in range(n_grid[0]):
        #for iy in range(n_grid[1]):
            #for iz in range(n_grid[2]):
                #position = ix*grid_vectors[0, :] + iy*grid_vectors[1, :] + iz*grid_vectors[2, :]
                #if position[linear_dim] >= linear_dim_upper_bottom:
                    #pot_on_grid[ix, iy, iz] = V_bias
                #elif position[linear_dim] <= linear_dim_lower_top:
                    #pot_on_grid[ix, iy, iz] = 0
                #else:
                    #pot_on_grid[ix, iy, iz] = (position[linear_dim]-linear_dim_lower_top)/(linear_dim_upper_bottom-linear_dim_lower_top) * V_bias
    
    print('- Gaussian blurring the potential')
    # Gaussian blur to round the derivative discontinuity of the piecewise linear potential
    dgrid = np.sqrt(grid_vectors[linear_dim, :].dot(grid_vectors[linear_dim, :]))
    normalized_gaussian_sigma = gaussian_blur_stdev/dgrid
    if pbc:
        blur_mode = 'wrap'
    else:
        blur_mode = 'nearest'
    pot_on_grid = gaussian_filter1d(pot_on_grid, normalized_gaussian_sigma, axis=linear_dim, mode=blur_mode, truncate=6.0)
    
    # Write to cube file
    print('- Writing potential to cube file: {}'.format(cube_file))
    pot_on_grid = pot_on_grid*eV_to_au
    cube_writer.write_cube_nonorthogonal(cube_file, n_grid, grid_vectors, pot_on_grid, np.array([0.0, 0.0, 0.0]),
                                        comment_line='CP2k external electrostatic potential')


def init_grid_(cell_vectors, cutoff, extended_fft_lengths):
    n_grid = get_grid_size(cell_vectors, cutoff, use_extended_fft_lengths=extended_fft_lengths)
    grid_vectors = np.zeros((3, 3))
    grid_vectors[0, :] = cell_vectors[0, :]/n_grid[0]
    grid_vectors[1, :] = cell_vectors[1, :]/n_grid[1]
    grid_vectors[2, :] = cell_vectors[2, :]/n_grid[2]
    grid_data = np.zeros(n_grid)
    return n_grid, grid_vectors, grid_data
