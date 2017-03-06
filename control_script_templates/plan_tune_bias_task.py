# -*- coding: utf-8 -*-
#!/usr/bin/python

import os, sys

from kpfm_sim_tasks import Tune_bias_task
from kpfm_sim_task_db import Task_control_db
from kpfm_sim_result_db import Result_db

a_nacl = 5.73
x = 0.0#0.5*a_nacl
y = 0.0
s = 18.0
V_start = 0.5
V_end = 2.0
V_step = 0.5
V = V_start

if len(sys.argv) == 4:
    task_db_file = sys.argv[1]
    result_db_file = sys.argv[2]
    global_res_db_file = sys.argv[3]
else:
    sys.exit("Usage: python plan_tune_bias_task.py <task_db_file> <result_db_file> <global_res_db_file>")

if not os.path.isfile(result_db_file):
    local_results = Result_db(result_db_file)
    with local_results:
        local_results.copy_atoms_data(global_res_db_file)

new_task = Tune_bias_task(x, y, s, V, V_start, V_end, V_step, result_db_file, global_res_db_file)

task_db = Task_control_db(task_db_file)
with task_db:
    task_db.write_task(new_task)
