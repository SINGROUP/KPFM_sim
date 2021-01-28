# -*- coding: utf-8 -*-
#! /usr/bin/env python

import os, sys, shutil

import numpy as np
import sqlite3

from ase.constraints import FixAtoms
from ase.io import read

from pycp2k.cp2k import CP2K
from cp2k_common_params import set_common_params_diagonalize, set_geo_opt_params, set_smearing_params
from cp2k_init import CP2k_init
from cp2k_input_tools import set_poisson_solver
from cp2k_output_tools import get_output_from_file, get_energy_from_output,\
                            get_charges_from_output
from kpfm_sim_result_db import Result_db
import kpfm_sim_global_constants as global_const

# beware now the direction goes in /y/ #

# for only visual checking #

just_check = False

# We start with everything in /y/ and NOT periodic only in /y/   #

# Parameters 
project_name = "AFM_try"
z_lay_bot = 2.1 # hight of the layers
z_lay_top = z_lay_bot
no_bot_lay = 1
no_top_lay = 1
initial_macro_dist = 6.135 # the bottom most Cu atom vs. the heighest C atom of C6HxClBr ... #
# full_geo_n_cell_file = "geometry-in.in" # first atoms are tip, then sample goes 
full_geo_file = "input.xyz" # first atoms are tip, then sample goes 
full_cell_file = "input.lvs" # first atoms are tip, then sample goes 
pbc = [True,False,True]
n_sample_atoms = 2
sample_bottom_atoms_file = "bottom_fixed_indices.txt" # No -1  # //not needed, since already prepared in the "in" file #
tip_top_atoms_file = "top_fixed_indices.txt" # No -1  # //not needed, since already prepared in the "in" file #

kpts = True
nd = 6# number of digits for rounding

db_filename = "glob_res/opt.db"
cp2k_output_file = project_name + global_const.cp2k_out_suffix
result_xyz_file = project_name + global_const.cp2k_xyz_suffix
result_wfn_file = project_name + global_const.cp2k_wfn_suffix(kpts=kpts)
wfn_storage_rel_path = "wfn_data" 

# Scan point 
# x, y and s not used at the moment #
x = 0.0
y = 0.0
s = initial_macro_dist 
V = 0.0


#===============================================================================
# sample & getting to know, what is what #
#atoms = read(full_geo_n_cell_file) #, format='xyz', index=-1) # tip is 1st, then sample
atoms = read(full_geo_file, format='xyz', index=-1) # tip is 1st, then sample
cell  = np.genfromtxt(full_cell_file)
atoms.set_cell(cell)
n_at  = atoms.get_global_number_of_atoms()
n_tip_atoms = n_at - n_sample_atoms

tip_top_atom_inds = np.loadtxt(tip_top_atoms_file) if n_tip_atoms > 2 else [ 0 ] # -1 -> human to python logic
sample_bottom_atom_inds = np.loadtxt(sample_bottom_atoms_file) # human logic here
sample_bottom_atom_inds +=  -1 # human -> python logic 

tip_apex_atom_idx = n_sample_atoms -1  #originally set by /y/, here it is lowesr = last atom of the tip #
tip_atom_inds = [i for i in range(n_tip_atoms)] # atoms of the tip goes as 1st

sample_atom_inds = [i for i in range(n_tip_atoms,n_at)] # bellow ; the top atoms are /y/ above the bottom layer - bottom atom last #
sample_top_atom_inds = [atom.index for atom in atoms if (atom.position[1] > atoms.get_positions()[n_sample_atoms-1,1]+z_lay_bot*no_bot_lay+0.8 and 0 <= atom.index < n_sample_atoms)]

atoms.set_pbc(pbc)

if just_check:
    print("just checking indices and the whole model:")
    print("n_sample_atoms         ",n_sample_atoms)
    print("sample_atom_inds       ",sample_atom_inds)
    print("sample_bottom_atom_inds",sample_bottom_atom_inds)
    print("sample_top_atom_inds   ",sample_top_atom_inds)
    print("n_tip_atoms            ",n_tip_atoms)
    print("total amount of at     ",atoms.get_number_of_atoms())
    print("tip_atom_inds          ",tip_atom_inds )
    print("tip_top_atom_inds      ",tip_top_atom_inds)
    print("tip_apex_atom_idx      ",tip_apex_atom_idx)

'''
# I don't need any translations at the moment #


# Combine sample and tip
tip_translation = np.zeros(3)
tip_translation[0] = 0.5*sample_model.cell[0,0] - tip_model.positions[tip_apex_atom_idx, 0]
tip_translation[1] = vacuum + initial_macro_dist - tip_height
tip_translation[2] = 0.5*sample_model.cell[2,2] - tip_model.positions[tip_apex_atom_idx, 2]
tip_model.translate(tip_translation)
sample_model.extend(tip_model)
tip_top_atom_inds = tip_top_atom_inds + n_sample_atoms
tip_apex_atom_idx = tip_apex_atom_idx + n_sample_atoms
'''


# constrained atoms # 

fix_indices = [ *tip_top_atom_inds, *sample_bottom_atom_inds ]

fix_all = FixAtoms(fix_indices)
atoms.set_constraint(fix_all)

fixed_atoms = atoms.constraints[0].get_indices()

# Make sure your initial tip-sample model looks correct
if just_check:
    from ase.visualize import view
    print("just checking indices and the whole model:")
    print("n_sample_atoms         ",n_sample_atoms)
    print("sample_atom_inds       ",sample_atom_inds)
    print("sample_bottom_atom_inds",sample_bottom_atom_inds)
    print("sample_top_atom_inds   ",sample_top_atom_inds)
    print("n_tip_atoms            ",n_tip_atoms)
    print("total amount of at     ",atoms.get_number_of_atoms())
    print("tip_atom_inds          ",tip_atom_inds )
    print("tip_top_atom_inds      ",tip_top_atom_inds)
    print("tip_apex_atom_idx      ",tip_apex_atom_idx)
    print("fix_indices            ",fix_indices)
    print("fixed_atoms            ",fixed_atoms)
    print("cell                   ",atoms.get_cell())
    print("pbc                    ",atoms.get_pbc())
    print("***********************")
    view(atoms)
    sys.exit()

print ("running the real system")

# Initialize the PYCP2k calculator object
cp2k_initializer = CP2k_init(project_name, atoms)
cp2k_calc = cp2k_initializer.init_desc_tip()

# Run CP2k to optimize geometry
cp2k_calc.run()


#===============================================================================
# IMPORTANT!
# This last part of this script is essential for setting up the result database
# file. It would probably be a good idea to move this to a separate module and
# call it from there. Or you can just copy it to your initialization scripts.

relaxed_model = read(result_xyz_file, format='xyz')
relaxed_model.cell = atoms.get_cell()
relaxed_model.pbc = atoms.get_pbc()

is_fixed = [(index in fixed_atoms) for index in range(relaxed_model.get_number_of_atoms())]

# Save information about the "role" of each atom in the model
atom_belongs_to = []
for atom_i in range(relaxed_model.get_number_of_atoms()):
    if atom_i in sample_atom_inds:
        if atom_i in sample_bottom_atom_inds:
            atom_belongs_to.append(("sample", "bottom"))
        elif atom_i in sample_top_atom_inds:
            atom_belongs_to.append(("sample", "top"))
        else:
            atom_belongs_to.append(("sample", "center"))
    elif atom_i in tip_atom_inds:
        if atom_i in tip_top_atom_inds:
            atom_belongs_to.append(("tip", "top"))
        elif atom_i == tip_apex_atom_idx:
            atom_belongs_to.append(("tip", "apex"))
        else:
            atom_belongs_to.append(("tip", "center"))
    else:
        atom_belongs_to.append(("unknown", "unknown"))

cp2k_output = get_output_from_file(cp2k_output_file)
energy = get_energy_from_output(cp2k_output)
charges = get_charges_from_output(cp2k_output)

for atom_i, atom_role in enumerate(atom_belongs_to):
    print ("atom {}: belongs to = {}; additional info = {}".format(atom_i, atom_role[0], atom_role[1]))
#sys.exit()

#===============================================================================
# Collect results to a database

results = Result_db(db_filename)
with results:
    print("debug results:", results)
    scan_point_id = results.write_scan_point( round(x,nd), round(y,nd), round(s,nd), round(V,nd), energy)
    print ("debug scan_point_id", scan_point_id)
    wfn_storage_file_name = global_const.wfn_storage_prefix + repr(scan_point_id) + \
                                global_const.cp2k_wfn_suffix(kpts=kpts)

    wfn_storage_path = os.path.join(wfn_storage_rel_path, wfn_storage_file_name)
    try:
        os.mkdir(wfn_storage_rel_path)
    except OSError:
        pass
    shutil.copy(result_wfn_file, wfn_storage_path)
    results.write_wf_data_path(scan_point_id, wfn_storage_path)
    results.write_output_file(scan_point_id, cp2k_output)
    results.write_atoms(relaxed_model, is_fixed, atom_belongs_to)
    results.write_atomic_geo(scan_point_id, relaxed_model, charges)
    results.write_unit_cell(scan_point_id, relaxed_model)
