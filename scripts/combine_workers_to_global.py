# -*- coding: utf-8 -*-
#!/usr/bin/python

import sys
from kpfm_sim_result_db import Result_db
import numpy as np
import os.path
from os import path
import shutil as shu
#from optparse import OptionParser

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
#                                   #
# WAS importent if 2 workers do     #
# one line scan; now rather use:    #
# copy_scan_points_and_wfn_files.py #
# together with combine_results.sh  #
# in slurm scripts directory ....   #
#                                   #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

# where are the database files and worker_paths:

kpts = False # do not change it, use the -k or --kpts flag before the other options

erm="""
The script is supposed to combine results from 2 ... .db files into (global results) main db file. It also copies the wave-function files into the chosen directory and also changes the pathway in the (global results) main db file accordingly.
Usage:
python combine_workers_to_global.py (optionally -k or --kpts if k-points are necessary) <worker1_res_db_file> <worker1_path> <worker2_res_db_file> <worker2_path> <to: global_results_db_file> <global_results_wfn_path>
"""

if len(sys.argv) == 7:
    db_f_1 = sys.argv[1]
    w1_pth = sys.argv[2]
    db_f_2 = sys.argv[3]
    w2_pth = sys.argv[4]
    o_db_f = sys.argv[5]
    wo_pth = sys.argv[6]
elif len(sys.argv) == 8:
    if (sys.argv[1] == '-k') or (sys.argv[1] == '--kpts'):
        kpts = True
    else:
        sys.exit( erm )
    db_f_1 = sys.argv[2]
    w1_pth = sys.argv[3]
    db_f_2 = sys.argv[4]
    w2_pth = sys.argv[5]
    o_db_f = sys.argv[6]
    wo_pth = sys.argv[7]
else:
   sys.exit( erm )

# not 
#parser = OptionParser()
#parser.add_option( "-k","--kpts",       action="store_true", default=False, help="necessary if k-points are used inside the calculation" )
#(options, args) = parser.parse_args()
#opt_dict = vars(options)
#kpts = opt_dict.kpts # False

                           
# just to stay on the save side:
shu.copyfile(o_db_f,o_db_f+"-bak")

# originally for debugging #
#db_f_1 = 'afm_copy.db'
#db_f_2 = 'afm_copy2.db'
#o_db_f = 'glob_res/afm_copy.db'
#w1_pth = 'worker_1'
#w2_pth = 'worker_2'
# end orignal #


# debuging:

debug = False

# also for debugging
initial = False
i1f='afm_res.db'
i2f='afm_res2.db'
of ='glob_res/afm.db'
if initial:
    bo = path.exists(o_db_f)
    b1 = path.exists(db_f_1)
    b2 = path.exists(db_f_2)
    if bo:
        shu.copyfile(o_db_f,o_db_f+"-bak")
    shu.copyfile(of,o_db_f)
    if b1:
        shu.copyfile(db_f_1,db_f_1+"-bak")
    shu.copyfile(i1f,db_f_1)
    if b1:
        shu.copyfile(db_f_2,db_f_2+"-bak")
    shu.copyfile(i2f,db_f_2)
# ended debugging section

# openening the databses:

i_db_1 = Result_db(db_f_1)
i_db_2 = Result_db(db_f_2) 

to_db = Result_db(o_db_f)

# preparing indices variables to nicely order eveything :

i1s = []; i2s = []; 
x1s = []; x2s =	[];
y1s = []; y2s =	[];
s1s = []; s2s =	[];
V1s = []; V2s =	[];
bwfc = True; ecl = [];

# getting all the indices first, then going for the :

with i_db_1:
    scan_points_1 = i_db_1.get_all_scan_point_entries()
    for scan_point in scan_points_1:
        from_id = scan_point[0]
        x = scan_point[1]
        y = scan_point[2]
        s = scan_point[3]
        V = scan_point[4]
        i1s.append(from_id)
        x1s.append(x)
        y1s.append(y)
        s1s.append(s)
        V1s.append(V)
    # exit()
    with i_db_2:
        scan_points_2 = i_db_2.get_all_scan_point_entries()
        for scan_point in scan_points_2:
            from_id = scan_point[0]
            x = scan_point[1]
            y = scan_point[2]
            s = scan_point[3]
            V = scan_point[4]
            i2s.append(from_id)
            x2s.append(x)
            y2s.append(y)
            s2s.append(s)
            V2s.append(V)
        #
        if debug:
            print ('i1s:',i1s); print ('x1s:',x1s); print ('y1s:',y1s); print ('s1s:',s1s); print ('V1s:',V1s)
            print ('i2s:',i2s); print ('x2s:',x2s); print ('y2s:',y2s); print ('s2s:',s2s); print ('V2s:',V2s)
        #
        fsx = np.array(np.concatenate([x1s,x2s]))
        fsy = np.array(np.concatenate([y1s,y2s]))
        fss = np.array(np.concatenate([s1s,s2s]))
        fsV = np.absolute(np.concatenate([V1s,V2s])) # only absolute part, so 0.0V goes as first
        fsVs= np.sign(np.concatenate([V1s,V2s])) # signs of voltage
        #
        ind = np.lexsort((fsx,fsy,fsV,fsVs,-1*fss)) # we need to go from highest to lowest
        l1 = len(i1s)
        #
        if debug:
            print ('ind:',ind); print ([(fsx[i],fsy[i],fsV[i],fss[i]) for i in ind]);
            #print ("scan_points_o",scan_points_o  ); #print ("scan_points_1",scan_points_1  ); #print ("scan_points_2",scan_points_2  )
            print ('l1',l1)
            #exit()
       #
        for i in ind:
            b_1          = True          if i < l1 else False ;
            scan_point   = scan_points_1[i] if b_1 else scan_points_2[i-l1]
            from_db      = i_db_1           if b_1 else i_db_2
            from_db_file = db_f_1           if b_1 else db_f_2
            if debug:
                 print ('i', i, "i-l1",i-l1)
                 print ('b_1', b_1)
                 print ('scan_point', scan_point)
                 print ('from_db'   , from_db)
            from_id = scan_point[0]
            x = scan_point[1]
            y = scan_point[2]
            s = scan_point[3]
            V = scan_point[4]
            energy = scan_point[5]
            with to_db:
                if to_db.get_scan_point_id(x,y,s,V) is None:
                    to_id = to_db.write_scan_point(x, y, s, V, energy)
                    if debug:
                        print ("to_id",to_id)
                    print("Copying scan point {} from {} to scan point {} in {}".format(from_id,
                            from_db_file, to_id, o_db_f))
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
                    wf_ext="wfn" if not kpts else "kp"
                    if wf_path is not None:
                        # "scan_point-10_-RESTART.wfn"
                        wf_path1 = wo_pth +"/scan_point-"+str(to_id)+"-RESTART."+wf_ext
                        f_path   = w1_pth +"/"+ wf_path if b_1 else w2_pth +"/"+ wf_path
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
    print ("All the Wavefunctions files from:",w1_pth,"and",w2_pth,"coppied succesfully." )
    print ("You should now run:")
    print ("rm -r "+w1_pth+"/wfn_data")
    print ("rm -r "+w2_pth+"/wfn_data")
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
