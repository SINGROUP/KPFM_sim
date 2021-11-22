# -*- coding: utf-8 -*-

import os, shutil
from abc import ABCMeta, abstractmethod

from ase.io import read

from kpfm_sim_result_db import Result_db
from cp2k_output_tools import get_output_from_file, get_energy_from_output,\
                                get_charges_from_output, get_forces_from_output
import kpfm_sim_global_constants as global_const

import numpy as np

style_metallic = True # new style of dealing with electrostatics, ommitting the original procedures and cython#

if not style_metallic:
    from axisym_pot_to_cube import axisym_pot_in_db_to_cube
    from macro_gap_efield import calc_macro_gap_efield


eps = 1.0e-13

nd = 6# number of digits for rounding

debug = False

def prepare_db_for_task(global_res_db_file, result_db_file, task_db_file):
    '''
    --- now moved to -> kpfm_sim_result_db.py <- stayed here for the dependency reasons ---
    prepare_db_for_task(global_res_db_file, result_db_file, task_db_file)
    adjust the result db file and the task db file with the all/last results from the global results file
    only the important parts (scan points and geometry) copied
    '''
    from_db = Result_db(global_res_db_file)
    to_db = Result_db(result_db_file)
    control_db = Result_db(task_db_file)
    gm = True
    with from_db:
        scan_points = from_db.get_all_scan_point_entries()
        #if scan_points == None:
        #    print("The global input file was not found; \n DEBUG: from_db",from_db,"\n DEBUG: global_res_db_file",global_res_db_file, "\n" )
        with to_db :
            with control_db :
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
                            global_res_db_file, to_id, result_db_file))
                        tmp, charges = from_db.extract_atoms_object(from_id, get_charges=True, get_model=gm);
                        pot_data_path = from_db.get_pot_data_path(from_id)
                        if gm:
                            atoms = tmp[0]	
                            model_part = tmp[1]
                            is_fixed = tmp[2]
                            full_model, pos_in_part  = from_db.get_model_part()
                            if debug:
                                print ("debug: model_part",model_part)
                                print ('debug: is_fixed', is_fixed)
                                print ('debug: full_model',full_model);
                                print ('debug: pos_in_part',pos_in_part);
                                print ('debug: pot_data_path',pot_data_path)
                            to_db.write_atoms(atoms, is_fixed, model_part, simplistic = True)
                            for i in range(len(full_model)):
                                to_db.write_model_part(full_model[i],pos_in_part[i])
                            gm = False;
                        else:
                            atoms = tmp
                        #print("debug: atoms",atoms)
                        #print("debug: atoms.positions", atoms.positions )
                        to_db.write_atomic_geo(to_id, atoms, charges)
                        to_db.write_unit_cell(to_id, atoms)
                        if pot_data_path is not None:
                            to_db.write_pot_data_path(to_id,pot_data_path)
                    ####
                    if control_db.get_scan_point_id(x, y, s, V) is None:
                        control_id = control_db.write_scan_point(x, y, s, V, energy)
                        print("Copying scan point {} from {} to scan point {} in {}".format(from_id,
                            global_res_db_file, control_id, task_db_file))
    print()
    print("results and tasks db files updated")


def copy_old_files_in_wrkdir(task_name):
    '''
    copy files in the directory when task is restarted - useful for a 
    '''
    xyz_file_name = task_name + global_const.cp2k_xyz_suffix
    output_file = task_name + global_const.cp2k_out_suffix
    input_file = task_name + global_const.cp2k_input_suffix
    # xyz file
    try:
        shutil.copyfile(xyz_file_name+"-2",xyz_file_name+"-3")
    except:
        pass # do nothing
    try:
        shutil.copyfile(xyz_file_name+"-1",xyz_file_name+"-2")
    except:
        pass # do nothing
    try:
        shutil.copyfile(xyz_file_name,xyz_file_name+"-1")
        if debug:
            print("debug: an old xyz file copied")
    except:
        pass # do nothing
    # output file
    try:
        shutil.copyfile(output_file+"-2",output_file+"-3")
    except:
        pass # do nothing
    try:
        shutil.copyfile(output_file+"-1",output_file+"-2") 
    except:
        pass # do nothing
    try:
        shutil.copyfile(output_file,output_file+"-1")
       	if debug:
       	    print("debug: an old output file copied")
    except:
        pass # do nothing
    # input file
    try:
        shutil.copyfile(input_file+"-2",input_file+"-3")
    except:
        pass # do nothing
    try:
        shutil.copyfile(input_file+"-1",input_file+"-2")
    except:
        pass # do nothing
    try:
        shutil.copyfile(input_file,input_file+"-1")
        if debug:
            print("debug: an old input file copied")
    except:
        pass # do nothing
    #-- end --#
    
## ********* THE ACTUAL TASKS CLASSES *********** ##
class Abstract_task(object, metaclass=ABCMeta):
    @abstractmethod
    def __init__(self, x, y, s, V, result_db_file, global_res_db_file, state, slurm_id, kpts=False, wfn=True):
        self.x = round(x,nd)
        self.y = round(y,nd)
        self.s = round(s,nd)
        self.V = round(V,nd)
        self.result_db_file = result_db_file
        self.global_res_db_file = global_res_db_file
        self.state = state
        self.slurm_id = slurm_id
        self.calc_initialized = False
        self.E_per_V = None
        self.kpts = kpts
        self.wfn = wfn
        self.tip_save_ind = None ## they are supposed to be named *safe* ##
        self.sample_save_ind = None


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

    def get_save_inds(self):
        '''
        will return two indices - tip top (fixed) atom and sample bottom (fixed) atom - so the geometries xyz can be compared
        ''' 
        if self.calc_initialized:
            return [self.tip_save_ind, self.sample_save_ind]
        else:
            raise Exception("Could not get atoms object from task because the calculation was not initialized.")

    def check_restart(self,pos,save_inds,xyz_file_name,zer=10**(-1*nd)):
        '''
        atomic positions = check_restart(save_inds,pos_z,xyz_file_name,zer):
        routinte to check if we should use the final geometry saved in the xyz_file_name using indices of atoms in save_inds is the same as positions in pos_z
        zer - maximal difference ~ 10^-6
        '''
        si = save_inds
        tp = pos[si[0]] # tip_positions
        sp = pos[si[1]] # sample_position
        if debug:
            print ("DEBUG: si, tp, sp", si, tp, sp );
        #
        try:
            neco = read(xyz_file_name)
            npos = neco.positions
            txyz = npos[si[0]]
            sxyz = npos[si[1]]
            if debug:
                print ("DEBUG: txyz, sxyz", txyz, sxyz );
                print ("DEBUG: dif tip", np.linalg.norm(tp - txyz) )
                print ("DEBUG: dif sample", np.linalg.norm(sp - sxyz) )
                print ("DEBUG: zero", zer)
            #
            if (np.linalg.norm(tp - txyz) < zer) and (np.linalg.norm(sp - sxyz) < zer):
                print("better geometry found, starting from the xyz file ...")
                if True:
                    nlnp = np.linalg.norm(pos-npos,axis=1)
                    print ("DEBUG: np.linalg.norm(pos-npos)",nlnp[(nlnp > zer)])
                return npos #True
            else:
                print("differences in geometries, going with the task-made geometry")
                return pos # False
        except:
            print("xyz file have not been found, going with the task-made geometry")
            return pos #False


    def get_atoms_object(self,save_inds=None,xyz_file_name=None):
        '''
        routine for obtaining the atoms object. If there exist xyz file with newer geometry, it is checked if it is from the same step and it updates the geometry from the xyz file
        '''
        if self.calc_initialized:
            ao = self.atoms
            if (save_inds is not None) and (xyz_file_name is not None):
                pos = self.check_restart(ao.positions,save_inds,xyz_file_name)
                ao.positions = pos
            #
            return ao
        else:
            raise Exception("Could not get atoms object from task because the calculation was not initialized.")


    # ---  Maybe to be completely omitted  ---
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
                        global_const.cp2k_wfn_suffix(kpts=self.kpts)), '.')
            print("\nRestart files found\n")
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


    def write_step_results_to_db(self, cp2k_output_path, kpts=False,bForces=False,wfnStore=True):
        if wfnStore:
            print ("debug: kpts =",kpts,"; self.kpts",self.kpts)
            wfn_file_name = self.task_name + global_const.cp2k_wfn_suffix(kpts=kpts)
            print ("debug: wfn_file_name" , wfn_file_name )
        cp2k_output = get_output_from_file(cp2k_output_path)
        energy = get_energy_from_output(cp2k_output)
        charges = get_charges_from_output(cp2k_output)
        forces = get_forces_from_output(cp2k_output) if bForces else None
        
        with self.results:
            scan_point_id = self.results.write_scan_point(self.x, self.y, self.s, self.V, energy)
            if scan_point_id is None:
                raise Exception("Tried to write scan point that already exists in the results database.")
            self.results.write_atomic_geo(scan_point_id, self.atoms, charges)
            self.results.write_unit_cell(scan_point_id, self.atoms)
            self.results.write_output_file(scan_point_id, cp2k_output)
            if forces is not None:
                self.results.write_atomic_forces(scan_point_id, forces)
            if wfnStore:
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
                                "_" + os.path.basename(self.worker_path) + global_const.cp2k_wfn_suffix(kpts=self.kpts)
        if debug:
            print("debug: wfn_file_name",wfn_file_name)
            print("debug: wfn_storage_file_name",wfn_storage_file_name)
            print("debug: self.kpts",self.kpts)
        wfn_storage_path = os.path.join(wfn_storage_base_path, wfn_storage_file_name)
        shutil.copy(wfn_file_name, wfn_storage_path)
        wfn_rel_storage_path = os.path.join(global_const.wfn_storage_folder, wfn_storage_file_name)
        return wfn_rel_storage_path


class Descend_tip_task(Abstract_task):
    def __init__(self, x, y, s, V, s_start, s_end, s_step, result_db_file,
                global_res_db_file, state = global_const.state_planned, slurm_id = None, kpts=False, wfn=True):
        if debug:
            print("debug: x,",x,"y",y,"s",s,"V",V,"s_start",s_start,"s_end",s_end,"s_step",s_step)
        Abstract_task.__init__(self, x, y, s, V, result_db_file, global_res_db_file, state, slurm_id, kpts=kpts, wfn=wfn)
        self.task_type = global_const.task_descend_tip
        self.s_start = round(s_start,nd)
        self.s_end = round(s_end,nd)
        self.s_step = round(s_step,nd)
        self.kpts = kpts


    # Move tip down
    def descend_tip(self, s_step = None):
        if s_step is None:
            s_step = self.s_step
        for atom in self.atoms:
            if atom.index in self.tip_atom_inds:
                atom.position[1] = atom.position[1] - s_step


    # Translate tip in x-y (x-z) plane # - this was the original procedure:
    def translate_sample(self, translation_vec):
        for atom in self.atoms:
            if atom.index in self.sample_atom_inds:
                atom.position[0] = atom.position[0] + translation_vec[0]
                atom.position[2] = atom.position[2] + translation_vec[1]

    # Translate tip in x-y (x-z) plane
    def translate_tip(self, translation_vec):
        for atom in self.atoms:
            if atom.index in self.tip_atom_inds:
                atom.position[0] = atom.position[0] + translation_vec[0]
                atom.position[2] = atom.position[2] + translation_vec[1]


    # Translate tip in z (y) axis according to the s_start
    def start_tip(self):
        if debug:
            print("debug in start_tip; self.s",self.s,"self.s_start",self.s_start,"s_start-s =",self.s_start-self.s)
        if self.s != self.s_start:
            print("s_start differs from s, moving the tip according to s_start")
            for atom in self.atoms:
                if atom.index in self.tip_atom_inds:
                    atom.position[1] = atom.position[1] + self.s_start-self.s
            self.s=self.s_start
        if debug:
            print("debug tip atoms adjusted")


    # Get the initial atomic model from the results database and move the tip
    # to its current location
    def init_calculation(self, task_name, project_path, worker_path):
        Abstract_task.init_calculation(self, task_name, project_path, worker_path)
        if debug:
            print ("DEBUG: inside init calculation")
        
        wfn_file_name = self.task_name + global_const.cp2k_wfn_suffix(kpts=self.kpts)
        if debug:
            print ("wfn_file_name",wfn_file_name)
        
        # Determine here which scan point to get
        with self.global_results:
            global_scan_points = self.global_results.get_larger_s_scan_points(self.x, self.y, self.s, self.V)
        with self.results:
            scan_points = self.results.get_larger_s_scan_points(self.x, self.y, self.s, self.V)
        
        if scan_points or global_scan_points:
            if scan_points:
                init_source_point = scan_points[0]
                init_source = self.results
                print("Initial larger s scan point from local results db.")
            else:
                init_source_point = global_scan_points[0]
                init_source = self.global_results
                print("Initial larger s scan point from global results db.")

            init_scan_point_id = init_source_point[0]
            if debug:
                print ("INIT-CALCULATION: init_scan_point_id", init_scan_point_id)
            init_s_step = round(init_source_point[1] - self.s,nd)
            if debug:
                print ("INIT-CALCULATION: init_s_step", init_s_step)
            with init_source:
                self.atoms = init_source.extract_atoms_object(init_scan_point_id)
                if debug:
                    print ("INIT-CALCULATION: init_source",init_source)
                    print ("INIT-CALCULATION: self.atoms", self.atoms)
                if self.atoms is None:
                    raise Exception("Could not obtain initial atoms object from the results database.")
                self.tip_atom_inds = init_source.get_model_part_atom_ids("tip")
                self.sample_atom_inds = init_source.get_model_part_atom_ids("sample")
                try:
                    print ("DEBUG: trying get save index (only 1) for checking restart - tip top")
                    self.tip_save_ind = init_source.get_model_part_atom_ids("tip",position_in_part="top")[0]
                except:
                    print ("DEBUG: cannot do: self.tip_save_inds = init_source.get_model_part_atom_ids('tip',position_in_part='top'")
                try:
                    print ("DEBUG: trying get save index (only 1) for checking restart - sample bottom")
                    self.sample_save_ind = init_source.get_model_part_atom_ids("sample",position_in_part="bottom")[0]
                except:
                    print ("DEBUG: cannot do: self.tip_save_inds = init_source.get_model_part_atom_ids('sample',position_in_part='bottom'")

                init_source.extract_wf_data(init_scan_point_id, wfn_file_name, self.project_path)
                #self.start_tip() moved to run_task.py
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
            # now we want to use already calculated geometries and do not do any optimization #
            # The original scripts are left here 
            with self.results:
                if style_metallic:
                    scan_points = self.results.get_all_s_scan_points(self.x,self.y)
                else:
                    can_points = self.results.get_all_s_scan_points(0.0, 0.0)
            if not scan_points:
                with self.global_results:
                    if style_metallic:
                        scan_points = self.global_results.get_all_s_scan_points(self.x,self.y)
                    else:
                        can_points = self.global_results.get_all_s_scan_points(0.0, 0.0)
                if not scan_points:
                    raise Exception("Could not find suitable initial scan point from the results databases.")
                init_source = self.global_results
                print("Initial translated scan point from global results db.")
            else:
                init_source = self.results
                print("Initial translated scan point from local results db.")
            if debug:
                print ("scan_points", scan_points)
            init_scan_point = scan_points[-1]
            if debug:
                print ("init_scan_point", init_scan_point)
            init_scan_point_id = init_scan_point[0]
            if debug:
                print ("init_scan_point_id", init_scan_point_id)
            init_s = init_scan_point[1]
            if debug:
                print ("init_s", init_s)
            if init_s < self.s:
                raise Exception("Could not find suitable initial scan point from the results databases.")
            tip_translation = (self.x, self.y)    
            with init_source:
                self.atoms = init_source.extract_atoms_object(init_scan_point_id)
                if self.atoms is None:
                    print ("==== something went wrong, so we are trying to get data from database with None id ==== ")
                    self.atoms = init_source.extract_atoms_object(None)
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
            if global_const.use_uniform_efield and not style_metallic:
                self.E_per_V = calc_macro_gap_efield(self.s, self.global_results)
            # Extracts external electrostatic potential from the results database
            # interpolating with respect to values of s and writes it to file "pot.cube"
            else:
                with self.global_results:
                    if not style_metallic:
                        axisym_pot_in_db_to_cube(global_const.cp2k_extpot_file, self.s,
                                                self.atoms, self.global_results)
                    else: # style_metallic: #
                        print ("going to recall the potential from precalculated E_field procedures" )
                        # TO BE WRITTEN -- HERE !!! #
        
        self.calc_initialized = True


    def next_step(self):
        if self.calc_initialized:
            if self.s - self.s_step >= self.s_end-eps:
                self.descend_tip()
                self.s =round( self.s - self.s_step,nd)
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
                global_res_db_file, state = global_const.state_planned, slurm_id = None, kpts=False):
        Abstract_task.__init__(self, x, y, s, V, result_db_file, global_res_db_file, state, slurm_id, kpts=kpts)
        self.task_type = global_const.task_tune_bias
        self.V_start = round(V_start,nd)
        self.V_end = round(V_end,nd)
        self.V_step = round(V_step,nd)
        self.kpts = kpts


    def init_calculation(self, task_name, project_path, worker_path):
        Abstract_task.init_calculation(self, task_name, project_path, worker_path)
        
        wfn_file_name = self.task_name + global_const.cp2k_wfn_suffix(kpts=self.kpts)
        
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
                    self.V = round(self.V + self.V_step,nd)
                    if abs(self.V) > eps:
                        return True
                    else:
                        return self.next_step()
            elif self.V_step < 0:
                if self.V + self.V_step >= self.V_end-eps:
                    self.V = round(self.V + self.V_step,nd)
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
