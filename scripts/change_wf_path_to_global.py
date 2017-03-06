# -*- coding: utf-8 -*-
#!/usr/bin/python

import os, sys
from kpfm_sim_result_db import Result_db
import kpfm_sim_global_constants as global_const

if len(sys.argv) == 3:
    result_db_file = sys.argv[1]
    worker_num = sys.argv[2]
else:
    sys.exit("Usage: python {} result_db_file worker_num".format(sys.argv[0]))

result_db = Result_db(result_db_file)
result_db_dir = os.path.dirname(result_db_file)
try:
    os.mkdir(global_const.wfn_storage_folder)
except OSError:
    pass

with result_db:
    scan_points = result_db.get_all_scan_point_entries()
    for scan_point in scan_points:
        scan_point_id = scan_point[0]
        wf_data_path = result_db.get_wf_data_path(scan_point_id)
        if wf_data_path is not None:
            old_wfn_rel_storage_path = os.path.join(result_db_dir, wf_data_path)
            new_wfn_storage_file_name = global_const.wfn_storage_prefix + \
                                        repr(scan_point_id) + \
                                        "_worker_{}".format(worker_num) + \
                                        global_const.cp2k_wfn_suffix
            new_wfn_rel_storage_path = os.path.join(global_const.wfn_storage_folder,
                                        new_wfn_storage_file_name)
            os.rename(old_wfn_rel_storage_path, new_wfn_rel_storage_path)
            result_db.delete_wf_data_path(scan_point_id)
            result_db.write_wf_data_path(scan_point_id, new_wfn_rel_storage_path)
