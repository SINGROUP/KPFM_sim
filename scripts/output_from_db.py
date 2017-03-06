# coding: utf-8
#!/usr/bin/python

import os, sys
import sqlite3

from kpfm_sim_task_db import Task_control_db
from kpfm_sim_result_db import Result_db

if len(sys.argv) != 4:
   sys.exit("Usage: python output_from_db.py row_id output_type db_filename\n"
            "Available output types: geo_opt, forces")
row_id = int(sys.argv[1])
output_type = sys.argv[2]
db_filename = sys.argv[3]

result_db = Result_db(db_filename)

if output_type == "geo_opt":
    with result_db:
        output = result_db.extract_output(row_id)
    with open("geo_opt.out", 'w') as out_file:
        out_file.write(output)
elif output_type == "forces":
    with result_db:
        output = result_db.extract_output(row_id, "forces")
    with open("calc_forces.out", 'w') as out_file:
        out_file.write(output)
else:
    raise Exception("The given output type does not exist.")
