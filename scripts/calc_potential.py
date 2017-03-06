# -*- coding: utf-8 -*-
#!/usr/bin/python

import sys, os

from kpfm_sim_result_db import Result_db
from cp2k_init import CP2k_init
from axisym_pot_to_cube import axisym_pot_in_db_to_cube
import kpfm_sim_global_constants as global_const

eps = 1.0e-13

if len(sys.argv) == 8:
    x = float(sys.argv[1])
    y = float(sys.argv[2])
    s = float(sys.argv[3])
    V = float(sys.argv[4])
    project_path = sys.argv[5]
    result_db_file = sys.argv[6]
    global_res_db_file = sys.argv[7]
else:
    sys.exit("Usage: python {} <x_tip> <y_tip> <s_tip> <V> <project_path> <result_db_file> <global_res_db_file>".format(sys.argv[0]))

project_name = "kpfm_x{}_y{}_s{}_V{}".format(x, y, s, V)
worker_dir = os.path.dirname(result_db_file)
worker_path = os.path.join(project_path, worker_dir)
wfn_file_name = project_name + global_const.cp2k_wfn_suffix
result_db_path = os.path.join(project_path, result_db_file)
global_res_db_path = os.path.join(project_path, global_res_db_file)

results = Result_db(result_db_path)
global_results = Result_db(global_res_db_path)
with results:
    scan_point_id = results.get_scan_point_id(x, y, s, V)
    atoms = results.extract_atoms_object(scan_point_id)
    results.extract_wf_data(scan_point_id, wfn_file_name, project_path)
    cp2k_initializer = CP2k_init(project_name, atoms)
    if abs(V) > eps:
        with global_results:
            axisym_pot_in_db_to_cube(global_const.cp2k_extpot_file, s, atoms, global_results)
        cp2k_calc = cp2k_initializer.init_calc_potential(V)
    else:
        cp2k_calc = cp2k_initializer.init_calc_potential()
    
    cp2k_calc.run()
    
    if abs(V) > eps:
        os.remove(global_const.cp2k_extpot_file)
