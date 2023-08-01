#!/usr/bin/python3 -u

import os
import numpy as np
import sys
import glob
import shutil
from ase import atoms
import math
from optparse import OptionParser

from E_field.bias_to_cube_c import save_cube
from kpfm_sim_result_db import Result_db, prepare_db_for_small


# ***************** Description ******************************** #
#                                                                #
# This script will scan the directory with Electric field        #
# potentials from <e_field_folder> for potentials stored in npy  #
# or cube format. Those which are not in the cube formate will   #
# transfered to cube format using the cube_head here and         #
# geometries from the <db_file> .                                #
#                                                                #
# ****************** constants ********************************* #
#                                                                #
R_BOHR = 0.529772 # Bohr-radius in angstroms                     #
eV_to_au = 0.03674932540261308                                   #
#                                                                #
#                                                                #
# ***************** PARAMETERS ********************************* #
debug            = True;
write2db         = False;
# ------------7----- Main parameters ---------------------------- #

db_file        = "glob_res/kpfm.db"
e_field_folder = "E_field"

inx = 320 ; # amount of division points in /x/ direction #
iny = 625 ; # amount of division points in /y/ direction #
inz = 288 ; # amount of division points in /z/ direction #
ddd = np.array([[8.9118733969396402E-002,0.0,0.0],
                [0.0,9.0706854378510865E-002,0.0],
                [0.0,0.0,8.5754541751130356E-002]])
ddd *= R_BOHR # ddd is in BOHRs, while goemetry and this script anyway works in Angstroms @
cube_head="""  234    0.000000    0.000000    0.000000
  320    0.089119    0.000000    0.000000
  625    0.000000    0.090707    0.000000
  288    0.000000    0.000000    0.085755
"""
# To get these ideally run a CP2K calc. with:  #
#DFT.PRINT.V_HARTREE_CUBE.Section_parameters = "ON"
#DFT.PRINT.V_HARTREE_CUBE.Stride   = "1 1 1"
# and then copy the lines 3-6 from the ***-v_hartree.cube #

cubelist = glob.glob(e_field_folder+"/field_*_final_pot.cube")

npylist  = glob.glob(e_field_folder+"/field_*_final_pot.cubenpy.npy")

if debug:
    print("cubelist",cubelist)
    print("npylist",npylist)

cubenum=[]
for ic in cubelist:
    cubenum.append(int(ic.split("_")[2]))

cubenum = np.array(cubenum, dtype=np.int32)

def final_pot_name(idx):
    return e_field_folder+"/field_"+str(idx)+"_final_pot.cube" ; # g_file not needed any_more ; use of idx instead #


npynum=[]
for ic in npylist:
    i_tmp = int(ic.split("_")[2])
    print (i_tmp)
    if i_tmp in cubenum:
        print("cube file already exists")
    else:
        print(i_tmp, "to be recovered, following")
        ft_db = Result_db(db_file)
        with ft_db:
            geom = ft_db.extract_atoms_object(i_tmp, get_charges=False, get_model=False);
            if debug: 
                print ("debug:geom", geom)
        pos = geom.positions;
        n_at = len(pos)
        mol_xyz = np.zeros((n_at,4))
        mol_xyz[:,:3]=pos
        mol_xyz[:,3] =geom.get_atomic_numbers()
        g_or = np.array([0.,0.,0.])
        g_vec = ddd.copy()
        Varr = np.load(ic)
        save_cube(Varr,mol_xyz,g_or,g_vec,file_path=final_pot_name(i_tmp), cube_head=cube_head)
    npynum.append(i_tmp)

npynum = np.array(npynum, dtype=np.int32)

if debug:
    print(cubenum)
    print(npynum)

print ("Finished, good, bye!")
