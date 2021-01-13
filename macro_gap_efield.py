# -*- coding: utf-8 -*-

import numpy as np
import grid_linear_interpolation_extparam as s_interp
#import matplotlib.pyplot as plt

def calc_macro_gap_efield(s, result_db):
    with result_db:
        lower_s_point, higher_s_point = result_db.get_closest_pot_scan_points(s)
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
    
    print("s = {}, s_lower = {}, s_higher = {}\n".format(s, s_lower, s_higher))
    
    with result_db:
        rs, zs, V_lower = result_db.get_external_potential(id_lower)
        rs, zs, V_higher = result_db.get_external_potential(id_higher)
    
    print("V_lower = ")
    print(V_lower)
    print("V_higher = ")
    print(V_higher)
    
    pot = s_interp.interpolate_2d(V_lower, V_higher, s, s_lower, s_higher)
    
    print("Potential interpolated with respect to s:")
    print(pot)
    print("")
    
    gap_z_values = zs[(zs > 0.0) & (zs < s)]
    gap_pot_values = pot[0, (zs > 0.0) & (zs < s)]
    pot_line_fit = np.polyfit(gap_z_values, gap_pot_values, 1, cov=False)
    efield = pot_line_fit[0]
    
    #pot_line = np.poly1d(pot_line_fit)
    #pot_line_values = pot_line(gap_z_values)
    
    #plt.plot(gap_z_values, gap_pot_values, 's', gap_z_values, pot_line_values, '-')
    #plt.show()
    
    return efield
