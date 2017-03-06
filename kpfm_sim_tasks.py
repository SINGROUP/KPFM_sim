# -*- coding: utf-8 -*-

import os, shutil
from abc import ABCMeta, abstractmethod

from ase.io import read

from kpfm_sim_result_db import Result_db
from cp2k_output_tools import get_output_from_file, get_energy_from_output,\
                                get_charges_from_output
from axisym_pot_to_cube import axisym_pot_in_db_to_cube
from macro_gap_efield import calc_macro_gap_efield
import kpfm_sim_global_constants as global_const

eps = 1.0e-13


class Abstract_task(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, x, y, s, V, result_db_file, global_res_db_file, state, slurm_id):
        self.x = x
        self.y = y
        self.s = s
        self.V = V
        self.result_db_file = result_db_file
        self.global_res_db_file = global_res_db_file
        self.state = state
        self.slurm_id = slurm_id
        self.calc_initialized = False
        self.E_per_V = None


    @abstractmethod
    def get_db_record(self):
        pass


    @abstractmethod
    def init_calculation(self, task_name, project_path, worker_path):
        self.task_name = task_name
        self.project_path = project_path
        self.worker_path = worker_path
        # Test for absolute path for backward compatibility
        if os.path.isabs(self.result_db_file):
            result_db_path = self.result_db_file
            global_res_db_path = self.global_res_db_file
        else:
            result_db_path = os.path.join(self.project_path, self.result_db_file)
            global_res_db_path = os.path.join(self.project_path, self.global_res_db_file)
        self.results = Result_db(result_db_path)
        self.global_results = Result_db(global_res_db_path)


    @abstractmethod
    def next_step(self):
        pass


    def get_atoms_object(self):
        if self.calc_initialized:
            return self.atoms
        else:
            raise Exception("Could not get atoms object from task because the calculation was not initialized.")


    def get_restart_data(self):
        if self.slurm_id is not None:
            restart_data_path = os.path.join(self.worker_path, repr(self.slurm_id))
        else:
            return False
        try:
            shutil.copy(os.path.join(restart_data_path, self.task_name +
                        global_const.cp2k_restart_suffix), '.')
            shutil.copy(os.path.join(restart_data_path, self.task_name +
                        global_const.cp2k_out_suffix), '.')
            shutil.copy(os.path.join(restart_data_path, self.task_name +
                        global_const.cp2k_wfn_suffix), '.')
            print "\nRestart files found\n"
        except IOError:
            return False
        return True


    def get_db_update(self):
        db_update = [self.state, self.x, self.y, self.s, self.V]
        return db_update


    def update_atoms_object(self, xyz_file):
        new_atoms = read(xyz_file, format='xyz')
        self.atoms.positions = new_atoms.get_positions()
        os.remove(xyz_file)


    def write_step_results_to_db(self, cp2k_output_path):
        wfn_file_name = self.task_name + global_const.cp2k_wfn_suffix
        cp2k_output = get_output_from_file(cp2k_output_path)
        energy = get_energy_from_output(cp2k_output)
        charges = get_charges_from_output(cp2k_output)
        
        with self.results:
            scan_point_id = self.results.write_scan_point(self.x, self.y, self.s, self.V, energy)
            if scan_point_id is None:
                raise Exception("Tried to write scan point that already exists in the results database.")
            self.results.write_atomic_geo(scan_point_id, self.atoms, charges)
            self.results.write_unit_cell(scan_point_id, self.atoms)
            self.results.write_output_file(scan_point_id, cp2k_output)
            wfn_rel_storage_path = self.__store_wf_data(scan_point_id, wfn_file_name)
            self.results.write_wf_data_path(scan_point_id, wfn_rel_storage_path)
        
        os.remove(cp2k_output_path)
        os.remove(self.task_name + global_const.cp2k_restart_suffix)


    def is_scan_point_in_db(self):
        with self.results:    
            scan_point_id = self.results.get_scan_point_id(self.x, self.y, self.s, self.V)
            if scan_point_id is None:
                return False
            else:
                return True


    def __store_wf_data(self, scan_point_id, wfn_file_name):
        wfn_storage_base_path = os.path.join(self.project_path, global_const.wfn_storage_folder)
        try:
            os.mkdir(wfn_storage_base_path)
        except OSError:
            pass
        wfn_storage_file_name = global_const.wfn_storage_prefix + repr(scan_point_id) + \
                                "_" + os.path.basename(self.worker_path) + global_const.cp2k_wfn_suffix
        wfn_storage_path = os.path.join(wfn_storage_base_path, wfn_storage_file_name)
        shutil.copy(wfn_file_name, wfn_storage_path)
        wfn_rel_storage_path = os.path.join(global_const.wfn_storage_folder, wfn_storage_file_name)
        return wfn_rel_storage_path


class Descend_tip_task(Abstract_task):
    def __init__(self, x, y, s, V, s_start, s_end, s_step, result_db_file,
                global_res_db_file, state = global_const.state_planned, slurm_id = None):
        Abstract_task.__init__(self, x, y, s, V, result_db_file, global_res_db_file, state, slurm_id)
        self.task_type = global_const.task_descend_tip
        self.s_start = s_start
        self.s_end = s_end
        self.s_step = s_step


    # Move tip down
    def descend_tip(self, s_step = None):
        if s_step is None:
            s_step = self.s_step
        for atom in self.atoms:
            if atom.index in self.tip_atom_inds:
                atom.position[1] = atom.position[1] - s_step


    # Translate tip in x-y (x-z) plane
    def translate_tip(self, translation_vec):
        for atom in self.atoms:
            if atom.index in self.sample_atom_inds:
                atom.position[0] = atom.position[0] + translation_vec[0]
                atom.position[2] = atom.position[2] + translation_vec[1]


    # Get the initial atomic model from the results database and move the tip
    # to its current location
    def init_calculation(self, task_name, project_path, worker_path):
        Abstract_task.init_calculation(self, task_name, project_path, worker_path)
        
        wfn_file_name = self.task_name + global_const.cp2k_wfn_suffix
        
        # Determine here which scan point to get
        with self.global_results:
            global_scan_points = self.global_results.get_larger_s_scan_points(self.x, self.y, self.s, self.V)
        with self.results:
            scan_points = self.results.get_larger_s_scan_points(self.x, self.y, self.s, self.V)
        
        if scan_points or global_scan_points:
            if scan_points:
                init_source_point = scan_points[0]
                init_source = self.results
                print "Initial larger s scan point from local results db."
            else:
                init_source_point = global_scan_points[0]
                init_source = self.global_results
                print "Initial larger s scan point from global results db."

            init_scan_point_id = init_source_point[0]
            init_s_step = init_source_point[1] - self.s
            with init_source:
                self.atoms = init_source.extract_atoms_object(init_scan_point_id)
                if self.atoms is None:
                    raise Exception("Could not obtain initial atoms object from the results database.")
                self.tip_atom_inds = init_source.get_model_part_atom_ids("tip")
                self.sample_atom_inds = init_source.get_model_part_atom_ids("sample")
                init_source.extract_wf_data(init_scan_point_id, wfn_file_name, self.project_path)
                self.descend_tip(init_s_step)
            
        elif abs(self.V) > eps:
            with self.results:
                init_scan_point_id = self.results.get_scan_point_id(self.x, self.y, self.s, 0.0)
            if init_scan_point_id is None:
                with self.global_results:
                    init_scan_point_id = self.global_results.get_scan_point_id(self.x, self.y, self.s, 0.0)
                if init_scan_point_id is None:
                    raise Exception("Tried to init bias voltage tuning at tip position that does not exist.")
                init_source = self.global_results
            else:
                init_source = self.results
            
            with init_source:
                self.atoms = init_source.extract_atoms_object(init_scan_point_id)
                if self.atoms is None:
                    raise Exception("Could not obtain initial atoms object from the results database.")
                self.tip_atom_inds = init_source.get_model_part_atom_ids("tip")
                self.sample_atom_inds = init_source.get_model_part_atom_ids("sample")
                init_source.extract_wf_data(init_scan_point_id, wfn_file_name, self.project_path)
        
        elif abs(self.V) < eps:
            with self.results:
                scan_points = self.results.get_all_s_scan_points(0.0, 0.0)
            if not scan_points:
                with self.global_results:
                    scan_points = self.global_results.get_all_s_scan_points(0.0, 0.0)
                if not scan_points:
                    raise Exception("Could not find suitable initial scan point from the results databases.")
                init_source = self.global_results
                print "Initial translated scan point from global results db."
            else:
                init_source = self.results
                print "Initial translated scan point from local results db."
            init_scan_point = scan_points[-1]
            init_scan_point_id = init_scan_point[0]
            init_s = init_scan_point[1]
            if init_s < self.s:
                raise Exception("Could not find suitable initial scan point from the results databases.")
            tip_translation = (self.x, self.y)    
            with init_source:
                self.atoms = init_source.extract_atoms_object(init_scan_point_id)
                if self.atoms is None:
                    raise Exception("Could not obtain initial atoms object from the results database.")
                self.tip_atom_inds = init_source.get_model_part_atom_ids("tip")
                self.sample_atom_inds = init_source.get_model_part_atom_ids("sample")
                self.translate_tip(tip_translation)
                if init_s-self.s > eps:
                    self.descend_tip(init_s-self.s)
        
        else:
            raise Exception("Could not find suitable initial scan point from the results databases.")
        
        if abs(self.V) > eps:
            #TODO: Add here alternative for creating piecewise linear/constant
            # electrostatic potential to be used with metallic tip+substrate systems
            # (call create_linear_pot_cube function in piecewise_linear_potential module).
            
            # Calculate average electric field in the gap of the macroscopic model from
            # the electrostatic potential
            if global_const.use_uniform_efield:
                self.E_per_V = calc_macro_gap_efield(self.s, self.global_results)
            # Extracts external electrostatic potential from the results database
            # interpolating with respect to values of s and writes it to file "pot.cube"
            else:
                with self.global_results:
                    axisym_pot_in_db_to_cube(global_const.cp2k_extpot_file, self.s,
                                                self.atoms, self.global_results)
        
        self.calc_initialized = True


    def next_step(self):
        if self.calc_initialized:
            if self.s - self.s_step >= self.s_end-eps:
                self.descend_tip()
                self.s = self.s - self.s_step
                if abs(self.V) > eps:
                    # Calculate average electric field in the gap of the macroscopic model from
                    # the electrostatic potential
                    if global_const.use_uniform_efield:
                        self.E_per_V = calc_macro_gap_efield(self.s, self.global_results)
                    # Extracts external electrostatic potential from the results database
                    # interpolating with respect to values of s and writes it to file "pot.cube"
                    else:
                        with self.global_results:
                            axisym_pot_in_db_to_cube(global_const.cp2k_extpot_file, self.s,
                                                        self.atoms, self.global_results)
                return True
            else:
                return False
        else:
            raise Exception("Task tried to take next step, but the calculation was not initialized.")


    def get_db_record(self):
        db_record = [self.slurm_id, self.task_type, self.state, self.x, self.y, self.s,
                    self.V, self.s_start, self.s_end, self.s_step, self.result_db_file,
                    self.global_res_db_file]
        return db_record


class Tune_bias_task(Abstract_task):
    def __init__(self, x, y, s, V, V_start, V_end, V_step, result_db_file,
                global_res_db_file, state = global_const.state_planned, slurm_id = None):
        Abstract_task.__init__(self, x, y, s, V, result_db_file, global_res_db_file, state, slurm_id)
        self.task_type = global_const.task_tune_bias
        self.V_start = V_start
        self.V_end = V_end
        self.V_step = V_step


    def init_calculation(self, task_name, project_path, worker_path):
        Abstract_task.init_calculation(self, task_name, project_path, worker_path)
        
        wfn_file_name = self.task_name + global_const.cp2k_wfn_suffix
        
        if abs(self.V-self.V_start) < eps:
            with self.results:
                init_scan_point_id = self.results.get_scan_point_id(self.x, self.y, self.s, 0.0)
            if init_scan_point_id is None:
                with self.global_results:
                    init_scan_point_id = self.global_results.get_scan_point_id(self.x, self.y, self.s, 0.0)
                if init_scan_point_id is None:
                    raise Exception("Tried to init bias voltage tuning at tip position that does not exist.")
                init_source = self.global_results
            else:
                init_source = self.results
        else:
            with self.results:
                init_scan_point_id = self.results.get_scan_point_id(self.x, self.y, self.s, self.V-self.V_step)
            if init_scan_point_id is None:
                raise Exception("Tried to init bias voltage tuning at tip position that does not exist.")
            init_source = self.results
                
        with init_source:
            self.atoms = init_source.extract_atoms_object(init_scan_point_id)
            init_source.extract_wf_data(init_scan_point_id, wfn_file_name, self.project_path)
        
        # Calculate average electric field in the gap of the macroscopic model from
        # the electrostatic potential
        if global_const.use_uniform_efield:
            self.E_per_V = calc_macro_gap_efield(self.s, self.global_results)
        # Extracts external electrostatic potential from the results database
        # interpolating with respect to values of s and writes it to file "pot.cube"
        else:
            with self.global_results:
                axisym_pot_in_db_to_cube(global_const.cp2k_extpot_file, self.s,
                                            self.atoms, self.global_results)
        
        self.calc_initialized = True


    def next_step(self):
        if self.calc_initialized:
            if self.V_step > 0:
                if self.V + self.V_step <= self.V_end+eps:
                    self.V = self.V + self.V_step
                    if abs(self.V) > eps:
                        return True
                    else:
                        return self.next_step()
            elif self.V_step < 0:
                if self.V + self.V_step >= self.V_end-eps:
                    self.V = self.V + self.V_step
                    if abs(self.V) > eps:
                        return True
                    else:
                        return self.next_step()
            return False
        else:
            raise Exception("Task tried to take next step, but the calculation was not initialized.")


    def get_db_record(self):
        db_record = [self.slurm_id, self.task_type, self.state, self.x, self.y, self.s,
                    self.V, self.V_start, self.V_end, self.V_step, self.result_db_file,
                    self.global_res_db_file]
        return db_record
