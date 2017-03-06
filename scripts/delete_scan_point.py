# coding: utf-8
#!/usr/bin/python

import os, sys
import sqlite3

from kpfm_sim_result_db import Result_db

if len(sys.argv) != 3:
   sys.exit("Usage: python delete_scan_point.py scan_point_id db_filename")
scan_point_id = int(sys.argv[1])
db_filename = sys.argv[2]

result_db = Result_db(db_filename)
with result_db:
    result_db.delete_scan_point(scan_point_id)
