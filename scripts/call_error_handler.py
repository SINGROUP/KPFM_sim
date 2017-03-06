# coding: utf-8
#!/usr/bin/python

import sys
from kpfm_sim_error_handling import task_termination_handling

if len(sys.argv) == 4:
    slurm_id = int(sys.argv[1])
    project_path = sys.argv[2]
    task_db_file = sys.argv[3]
else:
    sys.exit("Not enough arguments given to error handler caller.")

task_termination_handling(slurm_id, project_path, task_db_file)
