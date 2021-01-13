# coding: utf-8
#!/usr/bin/python

import os, sys
import sqlite3

from kpfm_sim_result_db import Result_db

eps = 1.0e-13

if (len(sys.argv) == 5) or (len(sys.argv) == 7):
    x_tip = float(sys.argv[1])
    y_tip = float(sys.argv[2])
    s_tip = float(sys.argv[3])
    db_filename = sys.argv[4]
    is_V_range = False
    if len(sys.argv) == 7:
        V_min = float(sys.argv[5])-eps
        V_max = float(sys.argv[6])+eps
        is_V_range = True
else:
    sys.exit("Usage: python delete_scan_point.py <x> <y> <s> <db_filename> [V_min] [V_max]")

result_db = Result_db(db_filename)
with result_db:
    scan_points = result_db.get_all_V_scan_points(x_tip, y_tip, s_tip)
    if not is_V_range:
        for scan_point in scan_points:
            print("Deleting scan point {} (V={})".format(scan_point[0], scan_point[1]))
            result_db.delete_scan_point(scan_point[0])
    else:
        for scan_point in scan_points:
            if (scan_point[1] >= V_min) and (scan_point[1] <= V_max):
                print("Deleting scan point {} (V={})".format(scan_point[0], scan_point[1]))
                result_db.delete_scan_point(scan_point[0])
