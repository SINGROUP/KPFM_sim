# -*- coding: utf-8 -*-
#!/usr/bin/python

import os, sys

#from KPFM_sim.kpfm_sim_tasks import Descend_tip_task
#from KPFM_sim.kpfm_sim_task_db import Task_control_db
#from KPFM_sim.kpfm_sim_result_db import Result_db
from kpfm_sim_tasks import Descend_tip_task
from kpfm_sim_task_db import Task_control_db
from kpfm_sim_result_db import Result_db

# set this part :

x = 0.0
y = 0.0
s = 6.135
V = 0.0
s_start = s
s_end = 2.00
s_step = 0.2

# don't touch the ater part, unless you know what you are doing

erm = "Usage: python plan_descend_tip_task.py (optionally -k or --kpts if k-points are necessary) <task_db_file> <result_db_file> <global_res_db_file>"
#nd = 6# number of digits for rounding -- not needed rounding of x, y, s, V, s_... is in : kpfm_sim_tasks.py
kpts = False
#print("debug: sys.argv:", sys.argv, "len(sys.argv)", len(sys.argv))

if len(sys.argv) == 4:
    task_db_file = sys.argv[1]
    result_db_file = sys.argv[2]
    global_res_db_file = sys.argv[3]
elif len(sys.argv) == 5:
    if (sys.argv[1] == "-k") or (sys.argv[1] == "--kpts"):
        kpts = True
    else:
        sys.exit(erm)
    task_db_file = sys.argv[2]
    result_db_file = sys.argv[3]
    global_res_db_file = sys.argv[4]
else:
    sys.exit(erm)

if not os.path.isfile(result_db_file):
    local_results = Result_db(result_db_file)
    with local_results:
        local_results.copy_atoms_data(global_res_db_file)

"""
class Descend_tip_task(Abstract_task):
    def __init__(self, x, y, s, V, s_start, s_end, s_step, result_db_file,
                    global_res_db_file, state = global_const.state_planned, slurm_id = None):
"""

new_task = Descend_tip_task(x, y, s, V, s_start, s_end, s_step, result_db_file, global_res_db_file, kpts=kpts)

task_db = Task_control_db(task_db_file)
with task_db:
    task_db.write_task(new_task)
