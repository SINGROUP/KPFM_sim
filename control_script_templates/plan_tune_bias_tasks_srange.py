# -*- coding: utf-8 -*-
#!/usr/bin/python

import os, sys

from kpfm_sim_tasks import Tune_bias_task
from kpfm_sim_task_db import Task_control_db
from kpfm_sim_result_db import Result_db

eps = 1.0e-13

a_nacl = 5.73
x = 0.5*a_nacl
y = 0.0

s_min = 19.0
s_max = 25.0

V_start = 0.5
V_end = 2.0
V_step = 0.5
V = V_start

if len(sys.argv) == 5:
    task_db_file = sys.argv[1]
    result_db_file = sys.argv[2]
    global_res_db_file = sys.argv[3]
    init_db = sys.argv[4]
else:
    sys.exit("Usage: python {} task_db_file result_db_file global_res_db_file init_db(local/global)".format(sys.argv[0]))

if init_db == "local":
    results = Result_db(result_db_file)
elif init_db == "global":
    results = Result_db(global_res_db_file)
    if not os.path.isfile(result_db_file):
        local_results = Result_db(result_db_file)
        with local_results:
            local_results.copy_atoms_data(global_res_db_file)
else:
    print("init_db input can be either \'local\' or \'global\'!")
    sys.exit(2)

with results:
    srange_scan_points = results.get_s_range_scan_points(x, y, s_min, s_max)
    bias_tuned_scan_points = results.get_s_range_scan_points(x, y, s_min, s_max, V_start)

bias_tuned_s_list = [scan_point[1] for scan_point in bias_tuned_scan_points]
bias_tuned_s_list.reverse()
if bias_tuned_s_list:
    s_bias_tuned = bias_tuned_s_list.pop()
else:
    s_bias_tuned = None

s_list = []
for scan_point in srange_scan_points:
    s = scan_point[1]
    if s_bias_tuned is None or s < s_bias_tuned-eps:
        s_list.append(s)
    else:
        if bias_tuned_s_list:
            s_bias_tuned = bias_tuned_s_list.pop()
        else:
            s_bias_tuned = None

task_db = Task_control_db(task_db_file)
with task_db:
    for s in s_list:
        new_task = Tune_bias_task(x, y, s, V, V_start, V_end, V_step, result_db_file, global_res_db_file)
        task_db.write_task(new_task)
