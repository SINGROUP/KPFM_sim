# coding: utf-8
#!/usr/bin/python

import sys, os
from subprocess import check_output, CalledProcessError

from kpfm_sim_tasks import Descend_tip_task, Tune_bias_task
from kpfm_sim_error_handling import cp2k_error_handling, write_task_info
from kpfm_sim_task_db import Task_control_db
import kpfm_sim_global_constants as global_const

from cp2k_init import CP2k_init, init_cp2k_restart


task_name = "kpfm_task"
task_types = ["descend_tip", "tune_bias"]
xyz_file_name = task_name + global_const.cp2k_xyz_suffix
restart_file = task_name + global_const.cp2k_restart_suffix
output_file = task_name + global_const.cp2k_out_suffix
default_state_constraint = global_const.state_planned

if len(sys.argv) == 4:
    task_db_file = sys.argv[1]
    project_path = sys.argv[2]
    slurm_id = int(sys.argv[3])
    task_type_constraint = None
    task_state_constraint = default_state_constraint
elif len(sys.argv) == 5:
    task_db_file = sys.argv[1]
    project_path = sys.argv[2]
    slurm_id = int(sys.argv[3])
    task_type_constraint = sys.argv[4]
    task_state_constraint = default_state_constraint
elif len(sys.argv) == 6:
    task_db_file = sys.argv[1]
    project_path = sys.argv[2]
    slurm_id = int(sys.argv[3])
    task_type_constraint = sys.argv[4]
    task_state_constraint = sys.argv[5]
else:
   sys.exit("Usage: python run_task.py <task_db_file> <project_path> <slurm_id> [type_constraint] [status_constraint]\n" \
            "Available task types: {}, {}".format(task_types[0], task_types[1]))


worker_dir = os.path.dirname(task_db_file)
worker_path = os.path.join(project_path, worker_dir)
task_db_file = os.path.basename(task_db_file)
task_db_path = os.path.join(worker_path, task_db_file)

task_db = Task_control_db(task_db_path)

# Make sure that the fetched task is reserved so that it is not fetched by
# another slurm job at the same time
with task_db:
    task_id, task = task_db.reserve_active_task(task_type_constraint, task_state_constraint)
    if task_id is None:
        raise Exception("No active tasks found using the constraints: "
                        "type = {}, state = {}.".format(task_type_constraint, task_state_constraint))
write_task_info(task_id, task_name)
task.init_calculation(task_name, project_path, worker_path)

cp2k_restart_exists = False
print "Type of the task = {}, state of the task = {}".format(task.task_type, task.state)
if task.state == global_const.state_waiting:
    cp2k_restart_exists = task.get_restart_data()

is_steps_left = True

# If CP2k calculation was terminated before, do a restart
if cp2k_restart_exists:
    cp2k_restart_calc = init_cp2k_restart(task_name, restart_file)
    task.slurm_id = slurm_id
    task.state = global_const.state_running
    with task_db:
        task_db.update_task_slurm_id(task_id, slurm_id)
        task_db.update_task_state(task_id, global_const.state_running)
    try:
        cp2k_restart_calc.run()
    except CalledProcessError:
        cp2k_error_handling(task_id, task_name, slurm_id, project_path, worker_path, task_db)
    try:
        task.update_atoms_object(xyz_file_name)
    except IOError:
        cp2k_error_handling(task_id, task_name, slurm_id, project_path, worker_path, task_db)
    task.write_step_results_to_db(cp2k_restart_calc.get_output_path())
    is_steps_left = task.next_step()
    with task_db:
        task_db.update_task(task_id, task)

# Calculate new steps of the loaded task until all done or time ends.
while is_steps_left:
    # Check if the scan point that is going to be calculated already
    # exists in the results database. In that case, skip it.
    scan_point_exists = task.is_scan_point_in_db()
    if scan_point_exists:
        print "Scan point x = {}, y = {}, s = {}, V = {} was found in the results database. " \
                "Skipping...".format(task.x, task.y, task.s, task.V)
        is_steps_left = task.next_step()
        continue
    
    # Check which task type was loaded from the database and choose
    # correct CP2k initialization based on that
    cp2k_initializer = CP2k_init(task_name, task.get_atoms_object())
    if isinstance(task, Descend_tip_task):
        cp2k_calc = cp2k_initializer.init_desc_tip(task.V, E_per_V=task.E_per_V)
    elif isinstance(task, Tune_bias_task):
        cp2k_calc = cp2k_initializer.init_tune_bias(task.V, E_per_V=task.E_per_V)
    else:
        raise Exception("No CP2k initialization implemented for the given task type.")
    task.slurm_id = slurm_id
    task.state = global_const.state_running
    with task_db:
        task_db.update_task_slurm_id(task_id, slurm_id)
        task_db.update_task_state(task_id, global_const.state_running)
    try:
        print "Running simulation at scan point x = {}, y = {}, s = {}, V = {}".format(task.x,
                task.y, task.s, task.V)
        #cp2k_calc.write_input_file()
        cp2k_calc.run()
    except CalledProcessError:
        cp2k_error_handling(task_id, task_name, slurm_id, project_path, worker_path, task_db)
    try:
        task.update_atoms_object(xyz_file_name)
    except IOError:
        cp2k_error_handling(task_id, task_name, slurm_id, project_path, worker_path, task_db)
    task.write_step_results_to_db(cp2k_calc.get_output_path())
    is_steps_left = task.next_step()
    with task_db:
        task_db.update_task(task_id, task)

with task_db:
    task_db.move_task_to_completed(task_id)
