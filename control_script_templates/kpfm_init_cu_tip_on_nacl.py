# -*- coding: utf-8 -*-
#! /usr/bin/env python

import os, sys, shutil

import numpy as np
import sqlite3

#from ase.visualize import view
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

# Parameters
project_name = "nacl_slab_cu_tip"
a_copper = 3.65
a_nacl = 5.73
vacuum = 6.0
initial_macro_dist = 25.0
sample_geo_file = "nacl_slab_10x4x10-pos-1.xyz"
sample_bottom_atoms_file = "bottom_fixed_indices.txt"
tip_geo_file = "cu_tip_100-pos-1.xyz"
tip_top_atoms_file = "top_fixed_indices.txt"
db_filename = "kpfm_cu_tip_on_nacl.db"
cp2k_output_file = project_name + global_const.cp2k_out_suffix
result_xyz_file = project_name + global_const.cp2k_xyz_suffix
result_wfn_file = project_name + global_const.cp2k_wfn_suffix
wfn_storage_rel_path = "wfn_data" 

# Scan point
x = 0.0
y = 0.0
s = initial_macro_dist
V = 0.0


#===============================================================================
# Create the structure here as an ASE Atoms object. It will be used later on to
# automatically create entries in the CP2K input.

# Load prerelaxed atomic models from files
sample_model = read(sample_geo_file, format='xyz')
sample_bottom_atom_inds = np.loadtxt(sample_bottom_atoms_file)
sample_bottom_atom_inds = sample_bottom_atom_inds - 1
tip_model = read(tip_geo_file, format='xyz')
tip_top_atom_inds = np.loadtxt(tip_top_atoms_file)
tip_top_atom_inds = tip_top_atom_inds - 1

# Extract relevant information from models
sample_model.center(vacuum=0.0)
sample_height = sample_model.cell[1, 1]
n_sample_atoms = sample_model.get_number_of_atoms()
sample_atom_inds = [i for i in range(n_sample_atoms)]
sample_top_atom_inds = [atom.index for atom in sample_model if atom.position[1] > sample_height-0.25*a_nacl]

tip_model.center(vacuum=0.0)
tip_height = tip_model.cell[1, 1]
n_tip_atoms = tip_model.get_number_of_atoms()
tip_atom_inds = [i + n_sample_atoms for i in range(n_tip_atoms)]
tip_apex_atom_idx = 0
for atom in tip_model:
    if atom.position[1] == 0.0:
        tip_apex_atom_idx = atom.index
        break

# Set correct unit cell and move sample to correct position
sample_model.cell[0, 0] = sample_model.cell[0, 0] + 0.5 * a_nacl
sample_model.cell[1, 1] = initial_macro_dist + 2.0 * vacuum
sample_model.cell[2, 2] = sample_model.cell[2, 2] + 0.5 * a_nacl
sample_model.set_pbc((True, False, True))
sample_model.translate((0, vacuum, 0))

# Combine sample and tip
tip_translation = np.zeros(3)
tip_translation[0] = 0.5*sample_model.cell[0,0] - tip_model.positions[tip_apex_atom_idx, 0]
tip_translation[1] = vacuum + initial_macro_dist - tip_height
tip_translation[2] = 0.5*sample_model.cell[2,2] - tip_model.positions[tip_apex_atom_idx, 2]
tip_model.translate(tip_translation)
sample_model.extend(tip_model)
tip_top_atom_inds = tip_top_atom_inds + n_sample_atoms
tip_apex_atom_idx = tip_apex_atom_idx + n_sample_atoms

# Fix atoms in the top of the tip and in the bottom of the sample
fixed_atoms = np.append(sample_bottom_atom_inds, tip_top_atom_inds)
fixed_atoms_int = [int(fixed_atom) for fixed_atom in fixed_atoms]
fix_bulk = FixAtoms(fixed_atoms_int)
sample_model.set_constraint(fix_bulk)

atoms = sample_model

# Make sure your initial tip-sample model looks correct
#view(atoms)
#sys.exit()


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

is_fixed = [(index in fixed_atoms_int) for index in range(relaxed_model.get_number_of_atoms())]

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

#for atom_i, atom_role in enumerate(atom_belongs_to):
    #print "atom {}: belongs to = {}; additional info = {}".format(atom_i, atom_role[0], atom_role[1])
#sys.exit()

#===============================================================================
# Collect results to a database

results = Result_db(db_filename)
with results:
    scan_point_id = results.write_scan_point(x, y, s, V, energy)
    wfn_storage_file_name = global_const.wfn_storage_prefix + repr(scan_point_id) + \
                                global_const.cp2k_wfn_suffix

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
