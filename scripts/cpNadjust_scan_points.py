# -*- coding: utf-8 -*-
#!/usr/bin/python

###################################################################
#                                                                 #
# script which allows you to change part of the model in database #
# especially useful when you want to change which atoms of        #
# already calculated system are supposed to be fixed .....        #
#                                                                 #
###################################################################

import sys
from kpfm_sim_result_db import Result_db
import numpy as np

if len(sys.argv) == 3:
    from_db_file = sys.argv[1]
    to_db_file = sys.argv[2]
else:
    sys.exit("Usage: python copy_scan_points.py from_db_file to_db_file")

from_db = Result_db(from_db_file)
to_db = Result_db(to_db_file)

gm = True
wf_path = None

with from_db:
    scan_points = from_db.get_all_scan_point_entries()
    with to_db:
        for scan_point in scan_points:
            from_id = scan_point[0]
            x = scan_point[1]
            y = scan_point[2]
            s = scan_point[3]
            V = scan_point[4]
            energy = scan_point[5]
            if to_db.get_scan_point_id(x, y, s, V) is None:
                to_id = to_db.write_scan_point(x, y, s, V, energy)
                print("Copying scan point {} from {} to scan point {} in {}".format(from_id,
                    from_db_file, to_id, to_db_file))
                    
                tmp, charges = from_db.extract_atoms_object(from_id, get_charges=True, get_model=gm);
                if gm:
                    atoms = tmp[0]
                    model_part = tmp[1]
                    is_fixed = tmp[2]
                    full_model, pos_in_part  = from_db.get_model_part()
                    print ("debug: model_part",model_part)
                    print ('debug: is_fixed', is_fixed)
                    print ('debug: full_model',full_model);
                    print ('debug: pos_in_part',pos_in_part);
                    #### -- CHANGE HERE BELLOW -- ####
                    ach = np.array([241,242,243,244,245]);# ach -=1 ;# atoms to change # already python logic
                    for ia in ach:
                        model_part[ia]=1 ;# 1-fixed tip; 2 middle tip; 3 tip apex ..
                        is_fixed[ia]=True
                    
                    ach = np.arange(248,777+1);# ach -=1 ;# atoms to change # already python logic
                    for ia in ach:
                        model_part[ia]=5 ;# 1-fixed tip; 2 middle tip; 3 tip apex, 4 sample centre, 5 sample bottom (in this example)
                        is_fixed[ia]=True
                    
                    print ("debug: model_part",model_part)
                    print ("debug: is_fixed",is_fixed)
                    #### --  changing part ends here -- ####
                    to_db.write_atoms(atoms, is_fixed, model_part, simplistic = True)
                    for i in range(len(full_model)):
                        to_db.write_model_part(full_model[i],pos_in_part[i])
                    gm = False;
                else:
                    atoms = tmp
                forces = from_db.get_atomic_forces(from_id)
                output = from_db.extract_output(from_id)
                #calc_forces_output = from_db.extract_output(from_id, "forces")
                #wf_path = from_db.get_wf_data_path(from_id)
                
                to_db.write_atomic_geo(to_id, atoms, charges)
                to_db.write_unit_cell(to_id, atoms)
                to_db.write_output_file(to_id, output)
                if forces is not None:
                    to_db.write_atomic_forces(to_id, forces)
                    to_db.write_calc_forces_output(to_id, calc_forces_output)
                if wf_path is not None:
                    to_db.write_wf_data_path(to_id, wf_path)

'''
def prepare_db_for_task(global_res_db_file, result_db_file, task_db_file):
    #
    #prepare_db_for_task(global_res_db_file, result_db_file, task_db_file)
    #adjust the result db file and the task db file with the all/last results from the global results file
    #only the important parts (scan points and geometry) copied
    #
    from_db = Result_db(global_res_db_file)
    to_db = Result_db(result_db_file)
    control_db = Result_db(task_db_file)
    gm = True
    with from_db:
        scan_points = from_db.get_all_scan_point_entries()
        #if scan_points == None:
        #    print("The global input file was not found; \n DEBUG: from_db",from_db,"\n DEBUG: global_res_db_file",global_res_db_file, "\n" )
        with to_db :
            with control_db :
                for scan_point in scan_points:
                    from_id = scan_point[0]
                    x = scan_point[1]
                    y = scan_point[2]
                    s = scan_point[3]
                    V = scan_point[4]
                    energy = scan_point[5]
                    if to_db.get_scan_point_id(x, y, s, V) is None:
                        to_id = to_db.write_scan_point(x, y, s, V, energy)
                        print("Copying scan point {} from {} to scan point {} in {}".format(from_id,
                            global_res_db_file, to_id, result_db_file))
                        tmp, charges = from_db.extract_atoms_object(from_id, get_charges=True, get_model=gm);
                        if gm:
                            atoms = tmp[0]	
                            model_part = tmp[1]
                            is_fixed = tmp[2]
                            full_model, pos_in_part  = from_db.get_model_part()
                            print ("debug: model_part",model_part)
                            print ('debug: is_fixed', is_fixed)
                            print ('debug: full_model',full_model);
                            print ('debug: pos_in_part',pos_in_part);
                            to_db.write_atoms(atoms, is_fixed, model_part, simplistic = True)
                            for i in range(len(full_model)):
                                to_db.write_model_part(full_model[i],pos_in_part[i])
                            gm = False;
                        else:
                            atoms = tmp
                        #print("debug: atoms",atoms)
                        #print("debug: atoms.positions", atoms.positions )
                        to_db.write_atomic_geo(to_id, atoms, charges)
                        to_db.write_unit_cell(to_id, atoms)
                    ####
                    if control_db.get_scan_point_id(x, y, s, V) is None:
                        control_id = control_db.write_scan_point(x, y, s, V, energy)
                        print("Copying scan point {} from {} to scan point {} in {}".format(from_id,
                            global_res_db_file, control_id, task_db_file))
    print()
    print("results and tasks db files updated")
'''
print("copied")
