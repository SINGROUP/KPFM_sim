# -*- coding: utf-8 -*-
#!/usr/bin/python

import sys
from kpfm_sim_result_db import Result_db
import os.path
from os import path
import shutil as shu

erm = "Usage: python copy_scan_points_and_wfn_files.py (optionally -k or --kpts if k-points are necessary) from_db_file <worker_path> to_db_file <results_wfn_path>"
kpts = False # do not change it, use the -k or --kpts flag before the other options
debug = True

if len(sys.argv) == 5:
    from_db_file = sys.argv[1]
    w1_pth = sys.argv[2]
    to_db_file = sys.argv[3]
    wo_pth = sys.argv[4]
elif len(sys.argv) == 6:
    if (sys.argv[1] == '-k') or (sys.argv[1] == '--kpts'):
        kpts = True
    else:
        sys.exit( erm )    
    from_db_file = sys.argv[2]
    w1_pth = sys.argv[3]
    to_db_file = sys.argv[4]
    wo_pth = sys.argv[5]
else:
    sys.exit( erm )

from_db = Result_db(from_db_file)
to_db = Result_db(to_db_file)

bwfc = True; ecl = []; # for error message, if problem with wfn copying
wf_ext="wfn" if not kpts else "kp"

if debug:
    print ("kpts",kpts)
    print ("wf_ext",wf_ext)

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
                    
                atoms, charges = from_db.extract_atoms_object(from_id, get_charges=True)
                forces = from_db.get_atomic_forces(from_id)
                output = from_db.extract_output(from_id)
                calc_forces_output = from_db.extract_output(from_id, "forces")
                wf_path = from_db.get_wf_data_path(from_id)
                
                to_db.write_atomic_geo(to_id, atoms, charges)
                to_db.write_unit_cell(to_id, atoms)
                to_db.write_output_file(to_id, output)
                if forces is not None:
                    to_db.write_atomic_forces(to_id, forces)
                    to_db.write_calc_forces_output(to_id, calc_forces_output)
                if wf_path is not None:
                    wf_path1 = wo_pth +"/scan_point-"+str(to_id)+"-RESTART."+wf_ext
                    f_path   = w1_pth +"/"+ wf_path
                    if debug:
                        print("wf_ext",wf_ext)
                        print("wf_path",wf_path)
                        print("wf_path1",wf_path1)
                        print("f_path",f_path)
                    try:
                        shu.copyfile(f_path,wf_path1)
                        print ("the wave-function file:",f_path,"was coppied to",wf_path1)
                    except:
                        bwfc = False
                        tmp  = "PROBLEM: CANNOT COPY - the wave-function file: "+f_path+" cannot be coppied to: "+wf_path1
                        print ( tmp )
                        ecl.append( tmp )
                    to_db.write_wf_data_path(to_id, wf_path1)

                                                
print()
if bwfc:
    print ("All the Wavefunctions files from:",w1_pth,"coppied succesfully." )
    print ("You should now run:")
    print ("rm -r "+w1_pth+"/wfn_data")
    print ("------------------------------")
    print ("to save the space")
else:
    print ("Some of the copying of wave function went wrong - You should have a look on these:")
    for ec in ecl:
        print (ec)
    print ("------------------------------")
    print ("the new wavefunction path is anyway written in the database")

print
print ("**** everything done ****")

