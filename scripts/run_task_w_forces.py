# coding: utf-8
#!/usr/bin/python

import sys, os
from subprocess import check_output, CalledProcessError

from kpfm_sim_tasks import Descend_tip_task, Tune_bias_task, copy_old_files_in_wrkdir
from kpfm_sim_error_handling import cp2k_error_handling, write_task_info
from kpfm_sim_task_db import Task_control_db
import kpfm_sim_global_constants as global_const

from optparse import OptionParser

from cp2k_init import CP2k_init, init_cp2k_restart

def get_output_path(wrk_dir,pr_nm):
        return wrk_dir + "/" + pr_nm + ".out"

# write out things for debugging : #

debug = False

# Main settings: #

task_name = "afm_task"
task_types = ["descend_tip"] #, "tune_bias"]

# end basic settings#

xyz_file_name = task_name + global_const.cp2k_xyz_suffix
restart_file = task_name + global_const.cp2k_restart_suffix
output_file = task_name + global_const.cp2k_out_suffix
default_state_constraint = global_const.state_planned
#kpts = False # False 


# setting paths, constrains and k-points:
#sys.exit("Usage: python run_task.py (optionally -k or --kpts if k-points are necessary) <task_db_file> <project_path> <slurm_id> [type_constraint] [status_constraint]\n" \
#            "Available task types: {}, {}".format(task_types[0], task_types[1]))
erm="2 *.db files are expected \n usage: python run_task.py (optionally -k or --kpts if k-points are necessary) -f <task_db_file> <project_path> -s <slurm_id> -t [type_constraint] -c [status_constraint]\n"

parser = OptionParser()
parser.add_option('-f', '--files', action ="store", type="string", nargs=2,
                    help='2 files expected <taks_db_file> <project_path>')
parser.add_option('-s', '--slurm_id', action ="store", type="string", nargs=1,
                    help='1 slurm_id expexted <slurm_id>')
parser.add_option('-t', '--type_constraint', action ="store", type="string",
                    help='optional type constraint [type constraint]')
parser.add_option('-c', '--status_constraint', action ="store", type="string",
                    help='optional status constraint [status constraint]')
parser.add_option('-k', '--kpts', action='store_true', default = False,
                    help="if k-points are necesssary for the CP2K calc")
parser.add_option('-w', '--no_wfn', action='store_false', default = True, 
                    help="do not store the wfn or kp file for late recalculations")
parser.add_option('-n', '--no_forces', action='store_false', default = True, 
                    help="do not store the forces")
(options,args) = parser.parse_args()

if debug:
    print("options:",options)
    print("args:",args)
    print("options.files",options.files)
    print("options.slurm_id",options.slurm_id)
    print("options.kpts",options.kpts)
    print("options.no_wfn",options.no_wfn)
    print("options.no_forces",options.no_forces)

if (options.files == None) or (options.slurm_id == None):
    sys.exit(erm)

kpts = options.kpts
bForces = options.no_forces
wfnStore = options.no_wfn

task_db_file = options.files[0]
project_path = options.files[1]
slurm_id = int(options.slurm_id)
task_type_constraint = options.type_constraint
task_state_constraint = options.status_constraint

worker_dir = os.path.dirname(task_db_file)
worker_path = os.path.join(project_path, worker_dir)
task_db_file = os.path.basename(task_db_file)
task_db_path = os.path.join(worker_path, task_db_file)

task_db = Task_control_db(task_db_path)

if debug: 
    print("kpts",kpts)
    print("worker_dir  ",worker_dir )
    print("worker_path ",worker_path)
    print("task_db_file",task_db_file)
    print("task_db_path",task_db_path)
    print("task_db     ",task_db_path)
    #sys.exit()

# Make sure that the fetched task is reserved so that it is not fetched by
# another slurm job at the same time
with task_db:
    task_id, task = task_db.reserve_active_task(task_type_constraint, task_state_constraint)
    if task_id is None:
        raise Exception("No active tasks found using the constraints: "
                        "type = {}, state = {}.".format(task_type_constraint, task_state_constraint))
write_task_info(task_id, task_name)
task.init_calculation(task_name, project_path, worker_path)

cp2k_restart_exists = False ## We try to put there automatic restart ##
print("Type of the task = {}, state of the task = {}".format(task.task_type, task.state))
if task.state == global_const.state_waiting:
    cp2k_restart_exists = task.get_restart_data()

is_steps_left = True
was_not_shift = True

## -- Going to be automatic -- ##
## If CP2k calculation was terminated before, do a restart
# if cp2k_restart_exists:
#    cp2k_restart_calc = init_cp2k_restart(task_name, restart_file)
#    task.slurm_id = slurm_id
#    task.state = global_const.state_running
#    with task_db:
#        task_db.update_task_slurm_id(task_id, slurm_id)
#        task_db.update_task_state(task_id, global_const.state_running)
#    try:
#        cp2k_restart_calc.run()
#    except CalledProcessError:
#        cp2k_error_handling(task_id, task_name, slurm_id, project_path, worker_path, task_db)
#    try:
#        task.update_atoms_object(xyz_file_name)
#    except IOError:
#        cp2k_error_handling(task_id, task_name, slurm_id, project_path, worker_path, task_db)
#    task.write_step_results_to_db(cp2k_restart_calc.get_output_path(), kpts=kpts, bForces = bForces, wfnStore = wfnStore )
#    is_steps_left = task.next_step()
#    with task_db:
#        task_db.update_task(task_id, task)
## -- Going to be automatic --  ##

# Calculate new steps of the loaded task until all done or time ends.
while is_steps_left:
    # Check if the scan point that is going to be calculated already
    # exists in the results database. In that case, skip it.
    if debug:
        print("debug: task.s", task.s, "task.s_start",task.s_start)
    if (task.s != task.s_start) and was_not_shift :
        print("shifting to starting point, debug: was_not_shift", was_not_shift )
        was_not_shift=False;
        task.start_tip()
        continue
    scan_point_exists = task.is_scan_point_in_db()
    if scan_point_exists:
        print("Scan point x = {}, y = {}, s = {}, V = {} was found in the results database. " \
                "Skipping...".format(task.x, task.y, task.s, task.V))
        is_steps_left = task.next_step()
        continue
    
    # copy files in the directory (if exists to older ones - for possible debugging or contol)
    copy_old_files_in_wrkdir(task_name)
    # Check which task type was loaded from the database and choose
    # correct CP2k initialization based on that
    cp2k_initializer = CP2k_init(task_name, task.get_atoms_object(save_inds=task.get_save_inds(),xyz_file_name=xyz_file_name))
    #print ("DEBUG: task.atoms",task.atoms ); exit()
   
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
        print("Running simulation at scan point x = {}, y = {}, s = {}, V = {}".format(task.x,
                task.y, task.s, task.V))
        #cp2k_calc.write_input_file()
        cp2k_calc.run()
    except CalledProcessError:
        cp2k_error_handling(task_id, task_name, slurm_id, project_path, worker_path, task_db)
    try:
        task.update_atoms_object(xyz_file_name)
    except IOError:
        cp2k_error_handling(task_id, task_name, slurm_id, project_path, worker_path, task_db)
    task.write_step_results_to_db( task_name+".out", kpts=kpts, bForces = bForces, wfnStore = wfnStore )
    #task.write_step_results_to_db(get_output_path(worker_dir,task_name))
    is_steps_left = task.next_step()
    with task_db:
        task_db.update_task(task_id, task)

with task_db:
    task_db.move_task_to_completed(task_id)
