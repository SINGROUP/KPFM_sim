# -*- coding: utf-8 -*-
#!/usr/bin/python

import sys, os

from kpfm_sim_result_db import Result_db
from cp2k_init import CP2k_init
from cp2k_output_tools import get_output_from_file, get_forces_from_output
from axisym_pot_to_cube import axisym_pot_in_db_to_cube
from macro_gap_efield import calc_macro_gap_efield
import kpfm_sim_global_constants as global_const

eps = 1.0e-13

project_name = "kpfm_atomic_forces"
wfn_file_name = project_name + global_const.cp2k_wfn_suffix

if len(sys.argv) == 4:
    project_path = sys.argv[1]
    result_db_file = sys.argv[2]
    global_res_db_file = sys.argv[3]
else:
    sys.exit("Usage: python calc_atomic_forces.py <project_path> <result_db_file> <global_res_db_file>")

worker_dir = os.path.dirname(result_db_file)
worker_path = os.path.join(project_path, worker_dir)
result_db_path = os.path.join(project_path, result_db_file)
global_res_db_path = os.path.join(project_path, global_res_db_file)

results = Result_db(result_db_path)
global_results = Result_db(global_res_db_path)

with results:
    no_forces_scan_points = results.get_no_forces_scan_points()
    previous_s = None
    for scan_point in no_forces_scan_points:
        scan_point_id = scan_point[0]
        s = scan_point[1]
        V = scan_point[2]
        print "Calculating atomic forces for scan point {}".format(scan_point_id)
        atoms = results.extract_atoms_object(scan_point_id)
        results.extract_wf_data(scan_point_id, wfn_file_name, project_path)
        cp2k_initializer = CP2k_init(project_name, atoms)
        
        if abs(V) > eps:
            if (previous_s is None) or (s != previous_s):
                previous_s = s
                
                if global_const.use_uniform_efield:
                    E_per_V = calc_macro_gap_efield(s, global_results)
                else:
                    E_per_V = None
                    with global_results:
                        axisym_pot_in_db_to_cube(global_const.cp2k_extpot_file, s, atoms, global_results)
            
            cp2k_calc = cp2k_initializer.init_calc_forces(V, E_per_V=E_per_V)
            
        else:
            cp2k_calc = cp2k_initializer.init_calc_forces()
        
        cp2k_calc.run()
        output = get_output_from_file(cp2k_calc.get_output_path())
        forces = get_forces_from_output(output)
        results.write_atomic_forces(scan_point_id, forces)
        results.write_calc_forces_output(scan_point_id, output)
        os.remove(cp2k_calc.get_output_path())
        os.remove(wfn_file_name)
