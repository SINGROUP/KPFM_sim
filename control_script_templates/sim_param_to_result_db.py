# coding: utf-8
#!/usr/bin/python

import sys
import sqlite3

if len(sys.argv) == 2:
    result_db_file = sys.argv[1]
else:
    sys.exit("Usage: python sim_param_to_result_db.py result_db_file")

parameters = {}

parameters["cutoff"] = 1000
parameters["vacuum"] = 6.0

with sqlite3.connect(result_db_file) as result_db_con:
    cur = result_db_con.cursor()
    
    cur.execute("CREATE TABLE IF NOT EXISTS sim_parameters(key TEXT, value TEXT)")
    result_db_con.commit()
    
    for key, value in parameters.items():
        cur.execute("INSERT INTO sim_parameters VALUES(?, ?)", (key, value))
    result_db_con.commit()
