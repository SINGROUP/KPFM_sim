#!/usr/bin/python3 -u

import os
import numpy as np
import sys
import glob
import shutil

from kpfm_sim_result_db import Result_db
from optparse import OptionParser

erm = "usage: python3 pot_path_to_result_db.py -o <output_db_file> and -i <input_db_file> or -e <path_to_folder_where_your_electric_fields_are>" 

parser = OptionParser()
parser.add_option('-o', '--output_db_file', action ="store", type="string", nargs=1, #not_mandatory
                    help='1 file expected <output_db_file>', default = "None" )
parser.add_option('-i', '--input_db_file', action ="store", type="string", nargs=1,
                    help='1 file expected <input_db_file>', default = "None")
parser.add_option('-e', '--ef_path', action='store', type="string", nargs=1, # not mandatory 
                    help="or 1 ef path expected: <worker path> <results_wfn_path>", default = "None")
(options,args) = parser.parse_args()

overshoot = 1000 # how many more indices can be above considered when going through the folder - 1000 seems to be pretty safe #

#print("D: options",options )

in_db_file   = options.input_db_file
db_file      = options.output_db_file
ef_pth       = options.ef_path

ef_pth     = ef_pth     if ef_pth     is not "None" else None
in_db_file = in_db_file if in_db_file is not "None" else None

ef_paths = [] # list where all the files will be stored

if db_file == "None":
    sys.exit(erm)
if (ef_pth is None) and (in_db_file is None):
    sys.exit(erm)
elif ef_pth is not None:
    all_f = os.listdir(ef_pth)
    mn = len(all_f) + overshoot 
    #print("D: all_f",all_f)
    for ifile in range(mn):
        f_name = 'field_'+str(ifile)+'_final_pot.cube' 
        #print("D: f_name",f_name)
        f_path = ef_pth+"/"+f_name
        #print("D: f_path",f_path)
        if f_name in all_f:
            ef_paths.append([ifile, f_path])
else:
    from_db = Result_db(in_db_file)
    with from_db:
        scan_points = from_db.get_all_scan_point_entries()
        for scan_point in scan_points:
             from_id = scan_point[0]
             pot_data_path = from_db.get_pot_data_path(from_id)
             if pot_data_path is not None:
                 ef_paths.append( [ from_id , pot_data_path ] )

#print("D: ef_paths",ef_paths)

ft_db = Result_db(db_file)

with ft_db:
    ft_db.create_tables()
    for i in range(len(ef_paths)) :
        idx, fp = ef_paths[i]
        #print ("D: idx,fp:",idx, fp)
        ft_db.write_pot_data_path(idx, fp)

print ("--done")
