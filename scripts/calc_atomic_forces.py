# -*- coding: utf-8 -*-
#!/usr/bin/python

import sys, os
import numpy as np

from kpfm_sim_result_db import Result_db
from cp2k_init import CP2k_init
from cp2k_output_tools import get_output_from_file, get_forces_from_output
from axisym_pot_to_cube import axisym_pot_in_db_to_cube
from macro_gap_efield import calc_macro_gap_efield
import kpfm_sim_global_constants as global_const

from optparse import OptionParser

eps = 1.0e-13
r   = 6 # number of digits for rounding
#kpts = False
erm = "Usage: python calc_atomic_forces.py -p <project_path> -r <result_db_file> -g <global_res_db_file> (-k/--kpts if k-pouints in cp2k calc) (-l/--location <x> <y> <V>)"
debug = True

project_name = "afm_atomic_forces"

parser = OptionParser()
parser.add_option('-g', '--global_db_file', action ="store", type="string", nargs=1,
                    help='1 file expected <global_res_db_file>', default = "None")
parser.add_option('-r', '--result_db_file', action ="store", type="string", nargs=1,
                    help='1 file expected <result_db_file>', default = "None" )
parser.add_option('-p', '--project_path', action ="store", type="string", nargs=1,
                    help='1 directory expected <project_path>', default = "None" )
parser.add_option('-k', '--kpts', action='store_true', default = False, # not mandatory
                    help="if k-points are necesssary for the CP2K calc")
parser.add_option('-l', '--location', action='store', type="float", nargs=3, # not mandatory -> useful for large scans and/or parallel runs
                    help="Not mandatory: 3 floats expected <x> <y> <V> of point for which you want to calculate the atomic forces" , default = (float("nan"),float("nan")))
(options,args) = parser.parse_args()


if debug:
    print ("debug: options, args",options, args)

project_path       = options.project_path
result_db_file     = options.result_db_file
global_res_db_file = options.global_db_file
kpts               = options.kpts
point              = options.location
b_loc              = False # if point[0] == float("nan") else True

if (project_path == "None") or (result_db_file == "None") or (global_res_db_file == "None"):
    sys.exit( erm )
if point[0] != float("nan"):
    b_loc = True
    project_name += "x"+str(round(point[0],6))+"_y"+str(round(point[1],6))+"_V"+str(round(point[2],6))
    print("project name updated:",project_name)

    
wfn_file_name = project_name + global_const.cp2k_wfn_suffix(kpts=kpts)

if debug:
    print ("kpts",kpts)
    print ("wfn_file_name",wfn_file_name)

worker_dir = os.path.dirname(result_db_file)
worker_path = os.path.join(project_path, worker_dir)
result_db_path = os.path.join(project_path, result_db_file)
global_res_db_path = os.path.join(project_path, global_res_db_file)

results = Result_db(result_db_path)
global_results = Result_db(global_res_db_path)

with results:
    no_forces_scan_points = results.get_no_forces_scan_points()
    if b_loc:
        loc_points= results.get_all_s_scan_points(point[0], point[1], V=point[2])
        lpids = np.array(loc_points[:][0],dtype=np.int32) # loc point ids
        if debug:
            print ("debug: loc_points",loc_points)
            print ("debug: lpids", lpids)
            print ("debug: no_forces_scan_points ", no_forces_scan_points)
    previous_s = None
    for scan_point in no_forces_scan_points:
        if b_loc:
            if scan_point[0] in lpids:
                print("going to calculate:", scan_point)
            else:
                print("point",scan_point,"not on the right location:",point,"skipping")
                continue
        scan_point_id = scan_point[0]
        s = scan_point[1]
        V = scan_point[2]
        print("Calculating atomic forces for scan point {}".format(scan_point_id))
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
        #output = get_output_from_file(cp2k_calc.get_output_path())
        output = get_output_from_file(project_name+".out")
        forces = get_forces_from_output(output)
        results.write_atomic_forces(scan_point_id, forces)
        results.write_calc_forces_output(scan_point_id, output)
        os.remove(project_name+".out")
        try: # there is a posssibility, that CP2K does not print the wave file
            os.remove(wfn_file_name)
        except:
            print("wfn_file_name:",wfn_file_name,"does not exist, passing.")
