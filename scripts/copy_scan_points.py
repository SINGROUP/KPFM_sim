# -*- coding: utf-8 -*-
#!/usr/bin/python

import sys
from kpfm_sim_result_db import Result_db

if len(sys.argv) == 3:
    from_db_file = sys.argv[1]
    to_db_file = sys.argv[2]
else:
    sys.exit("Usage: python copy_scan_points.py from_db_file to_db_file")

from_db = Result_db(from_db_file)
to_db = Result_db(to_db_file)

with from_db:
    scan_points = from_db.get_all_scan_point_entries()
    with to_db:
        for scan_point in scan_points:
            from_id = scan_point[0]
            x = scan_point[1]
            y = scan_point[2]
            s = scan_point[3]
            V = scan_point[4]
            energy = scan_point[5]
            if to_db.get_scan_point_id(x, y, s, V) is None:
                to_id = to_db.write_scan_point(x, y, s, V, energy)
                print("Copying scan point {} from {} to scan point {} in {}".format(from_id,
                    from_db_file, to_id, to_db_file))
                    
                atoms, charges = from_db.extract_atoms_object(from_id, get_charges=True)
                forces = from_db.get_atomic_forces(from_id)
                output = from_db.extract_output(from_id)
                calc_forces_output = from_db.extract_output(from_id, "forces")
                wf_path = from_db.get_wf_data_path(from_id)
                
                to_db.write_atomic_geo(to_id, atoms, charges)
                to_db.write_unit_cell(to_id, atoms)
                to_db.write_output_file(to_id, output)
                if forces is not None:
                    to_db.write_atomic_forces(to_id, forces)
                    to_db.write_calc_forces_output(to_id, calc_forces_output)
                if wf_path is not None:
                    to_db.write_wf_data_path(to_id, wf_path)
