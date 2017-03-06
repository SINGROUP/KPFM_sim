# coding: utf-8

import sys, shutil, os

from kpfm_sim_task_db import Task_control_db
import kpfm_sim_global_constants as global_const

memory_limit_str = "memory limit"
time_limit_str = "TIME LIMIT"
memory_limit_cause = "mem_limit"
time_limit_cause = "time_limit"


def cp2k_error_handling(task_id, task_name, slurm_id, project_path, worker_path, task_db):
    print "***\ncp2k_error_handling was called\n***"
    common_error_handling(task_id, task_name, slurm_id, project_path, worker_path, task_db)
    
    print "***\nCP2k calculation failed.\n***"
    sys.exit(0)


def task_termination_handling(slurm_id, project_path, task_db_file):
    print "***\ntask_termination_handling was called\n***"
    try:
        task_id, task_name = read_task_info()
    except IOError:
        raise Exception("Task info file containing task_id and task_name was not found.\n"
                "Either task was not found from the database in the first place or\n"
                "there is something wrong with writing or reading the info file.")
    
    worker_dir = os.path.dirname(task_db_file)
    worker_path = os.path.join(project_path, worker_dir)
    task_db_file = os.path.basename(task_db_file)
    task_db_path = os.path.join(worker_path, task_db_file)
    task_db = Task_control_db(task_db_path)
    
    common_error_handling(task_id, task_name, slurm_id, project_path, worker_path, task_db)
        
    print "***\nTask was terminated.\n***"
    sys.exit(0)


def common_error_handling(task_id, task_name, slurm_id, project_path, worker_path, task_db):
    slurm_output_path = os.path.join(project_path, "slurm-{}.out".format(slurm_id))
    termination_cause = find_cause_of_termination(slurm_output_path)
    try:
        with task_db:
            if termination_cause == time_limit_cause:
                task_db.update_task_state(task_id, "waiting")
            else:
                task_db.update_task_state(task_id, "stopped ({})".format(termination_cause))
    except:
        print "Error in updating the task state in error handling"
    
    lsdir = os.listdir('.')
    for elem in lsdir:
        print elem
    
    try:
        restart_data_path = os.path.join(worker_path, repr(slurm_id))
        os.mkdir(restart_data_path)
    except OSError:
        print "***\nCould not make directory for restart data to {}.\n***".format(restart_data_path)
        sys.exit(1)
    
    try:
        shutil.copy(task_name + global_const.cp2k_input_suffix,
                    restart_data_path)
        shutil.copy(task_name + global_const.cp2k_out_suffix,
                    restart_data_path)
        print "\nInput and output files saved to {}\n".format(restart_data_path)
    except IOError:
        print "***\nError when copying output file.\n***"
        sys.exit(1)
    
    try:
        shutil.copy(task_name + global_const.cp2k_restart_suffix,
                    restart_data_path)
        shutil.copy(task_name + global_const.cp2k_wfn_suffix,
                    restart_data_path)
        print "\nRestart files saved to {}\n".format(restart_data_path)
    except IOError:
        print "***\nError when copying restart files.\n***"
        sys.exit(1)


def find_cause_of_termination(slurm_file_path):
    termination_cause = "unknown"
    with open(slurm_file_path, 'r') as slurm_file:
        for line in slurm_file:
            if memory_limit_str in line:
                termination_cause = memory_limit_cause
                break
            if time_limit_str in line:
                termination_cause = time_limit_cause
                break
    return termination_cause


def write_task_info(task_id, task_name):
    with open(global_const.task_info_file, 'w') as info_file:
        info_file.write("{} {}".format(task_id, task_name))


def read_task_info():
    with open(global_const.task_info_file, 'r') as info_file:
        line = info_file.readline()
        split_line = line.split()
        task_id = int(split_line[0])
        task_name = split_line[1]
    return task_id, task_name
