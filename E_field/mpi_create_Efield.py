#!/usr/bin/python3 -u

import os
import numpy as np
import sys
import glob
import shutil
from ase import atoms
from ase.io import *
import math
from mpi4py import MPI

from bias_to_cube_c import create_biased_cube
from kpfm_sim_result_db import Result_db, prepare_db_for_small


# ***************** Description ******************************** #
#                                                                #
# This script will take a single geometry from a xxx.traj file   #
# and creates an (aproximate) potential between the metallic tip #
# and sample in the KPFM measurements.                           #
# Those here in-built functions can be part of outer scripts     #
#                                                                #
# ****************** constants ********************************* #
#                                                                #
R_BOHR = 0.529772 # Bohr-radius in angstroms                     #
eV_to_au = 0.03674932540261308                                   #
#                                                                #
#                                                                #
# ***************** PARAMETERS ********************************* #
create_potential = False;
write2db         = True;
# ----------------- Main parameters ---------------------------- #

glob_db_file= "glob_res/afm.db"
db_file="e_fields.db"

e_field_folder="E_field"

g_file = 'geometry_x0.0_y0.0.traj'; # from where the geometry is taken - in the future this is supposed to be from the data base
starting_idx = 2; # 0 already done !!!! # Max 21 per geometry_0.0_0.0

V_tip = 1.00 # at the moment this is tip-Voltage #

inx = 320 ; iny = 625 ; inz = 288 ;
ddd = np.array([[8.9118733969396402E-002,0.0,0.0],
                [0.0,9.0706854378510865E-002,0.0],
                [0.0,0.0,8.5754541751130356E-002]])
ddd *= R_BOHR # ddd is in BOHRs, while goemetry and this script anyway works in Angstroms @
cube_head="""  162    0.000000    0.000000    0.000000
  320    0.089119    0.000000    0.000000
  625    0.000000    0.090707    0.000000
  288    0.000000    0.000000    0.085755
"""

# **************** function definitions ********************* #

def final_pot_name(idx):
    return e_field_folder+"/"+g_file+"_"+str(idx)+"_final_pot.cube" ;

def give_atoms_object(ft_db):
    with ft_db:
        ft_db.create_tables()

# *************** the actual executives ********************* #

try: # make the folder where to store Efields, if it doesn't exists
    os.mkdir(e_field_folder)
    print("a folder made")
except:
    print("e field folder already exists")

max_id , max_pot_id = prepare_db_for_small(glob_db_file,db_file) # will create results 


# !!!!!! MPIS !!!!! 
if create_potential:
    comm = MPI.COMM_WORLD # communicator object containing all processes
    size = comm.Get_size()
    rank = comm.Get_rank()
    idx = max_pot_id+1+rank
    #
    if debug:
        print("I am rank %d in group of %d processes" % (rank, size))
        print ("rank",rank,"idx",idx)
    #
    ft_db = Result_db(db_file)
    with ft_db:
        geom = from_db.extract_atoms_object(from_id, get_charges=False, get_model=false);
        if debug: 
            print ("debug: geom")
    #geom=read(g_file,index=idx) #not necessary anymore!#
    create_biased_cube(geom,V_tip,final_pot_name=final_pot_name(idx),cube_head=cube_head);
    comm.barrier() # just to be absolutely sure
    MPI.Finalize() # will end up the mpi communication
# !!!!! End of mpi !!!!!

idx = 0 # TEST TEST TEST #


if write2db:
    with ft_db:
        #ft_db.create_tables()
        ft_db.write_pot_data_path(idx, final_pot_name(idx))

print ("done,done")
