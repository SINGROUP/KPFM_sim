#!/usr/bin/python3 -u

import os
import numpy as np
import sys
import glob
import shutil
from ase import atoms
import math
from optparse import OptionParser

from E_field.bias_to_cube_c import create_biased_cube
from kpfm_sim_result_db import Result_db, prepare_db_for_small


# ***************** Description ******************************** #
#                                                                #
# This script will prepare a separate db file <db_file> from     #
# your global repository <glob_db_file>                          #
# and creates an (aproximate) potential between the metallic tip #
# and sample in the KPFM measurements. It can create multiple    #
# potentials at the same if runned through MPI                   #
# UPDATE since Feb 18th, 2022 it uses procedure which size the   #
# atoms to 80% of the vdW sphere according to OPLS force-field.  # 
#                                                                #
#                    USAGE:                                      #
#                                                                #
# Copy to your running directory and adjust the main parameters: #
# <glob_db_file>, <db_file>, <e_field_folder>, <cc> ...          #
# <inx>, <iny>, <inz>, <ddd>, and ideally also the <cube head>   #
# there is a furher hint, how to get those parameters just       #
# bellow them.                                                   #
# In your running (slurm) script adjust the PYTHON_PATH to your  #
# KPFM_sim folder and then run:                                  #
# python3 mpi_create_Efield # for a serial run                   #
# or:                                                            #
# srun python3 mpi_create_Efield --mpi # for mpi4py parellized   #
# (use 'mpiexec -n XX' or 'mpirun -n XX' if your super-computer  #
#  do not know srun )                                            #
# You can also use a parser options -f "cube" xor -f "npy" xor   #
# -f "both", which tells you how you store the final results.    #
# for some runs we experienced problems while writing cube files #
#                                                                #
# Optionally:                                                    #
# Those here in-built functions can be part of outer scripts     #
#                                                                #
# ****************** constants ********************************* #
#                                                                #
R_BOHR = 0.529772 # Bohr-radius in angstroms                     #
eV_to_au = 0.03674932540261308                                   #
#                                                                #
#                                                                #
# ***************** PARAMETERS ********************************* #
debug            = False;
create_potential = True;
write2db         = True;
# ----------------- Main parameters ---------------------------- #

glob_db_file   = "glob_res/afm.db"
db_file        = "e_fields.db"
e_field_folder = "E_field"

cc = 11.0 # !!! IMPORTANT !!! the maximum "height" of the bottom (sample) atoms ; because of xzy - height refers to y in xyz file !!! IMPORTANT !!! #

inx = 320 ; # amount of division points in /x/ direction #
iny = 625 ; # amount of division points in /y/ direction #
inz = 288 ; # amount of division points in /z/ direction #
ddd = np.array([[8.9118733969396402E-002,0.0,0.0],
                [0.0,9.0706854378510865E-002,0.0],
                [0.0,0.0,8.5754541751130356E-002]])
ddd *= R_BOHR # ddd is in BOHRs, while goemetry and this script anyway works in Angstroms @
cube_head="""  162    0.000000    0.000000    0.000000
  320    0.089119    0.000000    0.000000
  625    0.000000    0.090707    0.000000
  288    0.000000    0.000000    0.085755
"""
# To get these ideally run a CP2K calc. with:  #
#DFT.PRINT.V_HARTREE_CUBE.Section_parameters = "ON"
#DFT.PRINT.V_HARTREE_CUBE.Stride   = "1 1 1"
# and then copy the lines 3-6 from the ***-v_hartree.cube #

# **************** parser and additional params ************* #

starting_idx = 0; # Leave 0 unless you know what you are doing and you want to calculate many   #

V_tip = 1.00 # leave to 1.00 ; at the moment this is tip-Voltage #

parser = OptionParser()
parser.add_option('--mpi', action='store_true', default = False, 
                    help="allows to run MPI")
parser.add_option('-f',"--data_format" , action="store" , type="string",
                      help="Specify the output format of the vector and scalar "
                      "field. Supported formats are: cube,npy or both", default="cube")
(options,args) = parser.parse_args()
mpi_b           = options.mpi; # boolen - run mpi? #
save            = options.data_format

ml = 100000 ; # maximal length of an array -- probably not needed anymore but not tested #
# note: the procedure could calculate the potential, just for some x & y point, but it is not adapted or tested, yet #

# **************** function definitions ********************* #

def final_pot_name(idx):
    return e_field_folder+"/field_"+str(idx)+"_final_pot.cube" ; # g_file not needed any_more ; use of idx instead #

#def give_atoms_object(ft_db):
#    with ft_db:
#        ft_db.create_tables()

def run_pre_calc(e_field_folder,glob_db_file,db_file):
    try: # make the folder where to store Efields, if it doesn't exists
        os.mkdir(e_field_folder)
        print("a folder made")
    except:
        print("e field folder already exists")
    
    no_pot_ids = prepare_db_for_small(glob_db_file,db_file) # will create results 
    return no_pot_ids

def one_create_potential(db_file,rank,idx, cc=11.0):
    '''
    db_file - from where you will take the geometry;
    rank - actually not really needed any more, left here for debugging
    index - index that have not been calculated yet and that is supposed to be calculated in this exact run
    cc - central height to distinguish between tip & sample
    '''
    #
    if debug:
        print("I am rank %d in group of %d processes; debug: index %d" % (rank, size, idx))
    ft_db = Result_db(db_file)
    with ft_db:
        geom = ft_db.extract_atoms_object(idx, get_charges=False, get_model=False);
        if debug: 
            print ("debug:geom", geom)
    create_biased_cube2(geom,V_tip,final_pot_name=final_pot_name(idx),cube_head=cube_head, cc=cc,idx=idx, save=save);
    return idx

def write_to_db(db_file,idx):
    ft_db = Result_db(db_file)
    with ft_db:
        #ft_db.create_tables()
        ft_db.write_pot_data_path(idx, final_pot_name(idx))

# *************** the actual executives ********************* #

if mpi_b:
    print ("running MPI")
    #
    from mpi4py import MPI
    #
    comm = MPI.COMM_WORLD # communicator object containing all processes
    size = comm.Get_size()
    rank = comm.Get_rank()
    l=0; idx0=-1 ##
    if rank == 0 :
        no_pot_ids = run_pre_calc(e_field_folder,glob_db_file,db_file)
        l = len(no_pot_ids)
        print ("Points for which the potential already was not calculated",no_pot_ids)
        print ("number of still not caclaculated potentials:", l)
        for ir in range(1,size): # sending to all possible ranks! -- ir #
            comm.send(int(no_pot_ids[starting_idx+ir]),dest=ir)
        idx0 = int(no_pot_ids[starting_idx])
    for ir in range(1,size): # sending to all possible ranks! -- ir #
        if rank == ir:
            idx0 = comm.recv(source=0)
    comm.barrier() # just to be absolutely sure
    l = comm.bcast(l, root=0)
    comm.barrier()
    if rank+starting_idx  <= l and idx0 > -1: # to check, that the everything is not already calculated
        if create_potential:
            idx = one_create_potential(db_file,rank,idx0, cc=cc) # only this is really parallelized
        if write2db and idx > 0:
            write_to_db(db_file,idx)
    print ("rank",rank,"idx",idx," -- DONE, DONE")

else:
    print ("running serial")
    size = 0
    rank = 0
    idx =  0 # !!!! just a check !!!! - if 0 and is not changed, then nothing is written #
    no_pot_ids = run_pre_calc(e_field_folder,glob_db_file,db_file)
    l = len(no_pot_ids)
    print ("Points for which the potential already was not calculated",no_pot_ids)
    print ("number of still not caclaculated potentials:", l)
    if rank+starting_idx <= l:
        idx0 = int(no_pot_ids[rank+starting_idx]) # here different idx is used, since idx is used to be shown that the one_create_potential is working fine
        if create_potential:
            idx = one_create_potential(db_file,rank,idx0, cc=cc)
        if write2db and idx > 0:
            write_to_db(db_file,idx)
    print ("rank",rank,"idx",idx," -- DONE, DONE")

print ("done,done")
