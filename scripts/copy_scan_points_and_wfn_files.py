# -*- coding: utf-8 -*-
#!/usr/bin/python

import sys
from kpfm_sim_result_db import Result_db, copy_db_ft
import os.path
from os import path
import shutil as shu


from optparse import OptionParser

erm = "Usage: python copy_scan_points_and_wfn_files.py -i from_db_file -o to_db_file (-w <worker_path>  <results_wfn_path> ) (--kpts)"
debug = True
parser = OptionParser()
parser.add_option('-i', '--input_db_file', action ="store", type="string", nargs=1,
                    help='1 file expected <input_db_file>', default = "None")
parser.add_option('-o', '--output_db_file', action ="store", type="string", nargs=1,
                    help='1 file expected <output_db_file>', default = "None" )
parser.add_option('-k', '--kpts', action='store_true', default = False,
                    help="if k-points are necesssary for the CP2K calc")
parser.add_option('-w', '--wfn_paths', action='store', type="string", nargs=2, # not mandatory 
                    help="2 paths expected: <worker path> <results_wfn_path>", default = ("None","None"))
(options,args) = parser.parse_args()


from_db_file = options.input_db_file
to_db_file   = options.output_db_file
kpts         = options.kpts
w1_pth       = options.wfn_paths[0]
wo_pth       = options.wfn_paths[1]

w1_pth = w1_pth if w1_pth is not "None" else None
wo_pth = wo_pth if wo_pth is not "None" else None

if (from_db_file == "None") or (to_db_file == "None" ):
    print(erm + "\n going to exit ...")
    exit()  

#bwfc = True; ecl = []; # for error message, if problem with wfn copying
wf_ext="wfn" if not kpts else "kp"

if debug:
    print ("w1_pth, wo_pth :",w1_pth , wo_pth )
    print ("kpts",kpts)
    print ("wf_ext",wf_ext)

bwfc, ecl = copy_db_ft(from_db_file, to_db_file, w1_path=w1_pth, wo_path=wo_pth, wf_ext=wf_ext)

if bwfc:
    if w1_pth is not None:
        print ("All the Wavefunctions files from:",w1_pth,"coppied succesfully." )
        print ("You should now run:")
        print ("rm -r "+w1_pth+"/wfn_data")
        print ("------------------------------")
        print ("to save the space")
    else:
        print ("No Wavefunction copied, cheers!")
else:
    print ("Some of the copying of wave function went wrong - You should have a look on these:")
    for ec in ecl:
        print (ec)
    print ("------------------------------")
    print ("the new wavefunction path is anyway written in the database")

print
print ("**** everything done ****")

