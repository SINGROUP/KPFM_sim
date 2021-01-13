# coding: utf-8

import sqlite3
from kpfm_sim_tasks import Descend_tip_task, Tune_bias_task

eps = 1.0e-13


# Class for storing information on separate KPFM simulation tasks to a
# database. Simulation tasks can be planned and fetched for execution
# in a batch job later.
class Task_control_db(object):
    def __init__(self, filename):
        print("Initializing task database")
        self.db_filename = filename
        try:
            self.db_con = self.open_db()
            self.create_tables()
        except sqlite3.Error as e:
            print("An error occurred:", e.args[0])
        finally:
            self.close_db()
        
        
    def __del__(self):
        print("Closing task database connection and deleting database object.")
        self.close_db()


    def __enter__(self):
        self.db_con = self.open_db()
        self.db_con.row_factory = sqlite3.Row
        self.db_con.text_factory = str
        
        
    def __exit__(self, exc_type, exc_value, traceback):
        self.close_db()
        if exc_type is None:
            return True
        else:
            return False


    def open_db(self):
        con = sqlite3.connect(self.db_filename, detect_types = sqlite3.PARSE_DECLTYPES)
        return con


    def close_db(self):
        if self.db_con:
            self.db_con.close()


    def create_tables(self):
        sqlite3.register_adapter(bool, int)
        sqlite3.register_converter("BOOLEAN", lambda v: bool(int(v)))
        cur = self.db_con.cursor()
        cur.execute("CREATE TABLE IF NOT EXISTS active_tasks(id INTEGER PRIMARY KEY, "
                    "slurm_id INTEGER, type TEXT, state TEXT, x REAL, y REAL, s REAL, "
                    "V REAL, param_start REAL, param_end REAL, param_step REAL, "
                    "result_db TEXT, global_res_db TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS completed_tasks(id INTEGER PRIMARY KEY, "
                    "type TEXT, x REAL, y REAL, s REAL, V REAL, "
                    "param_start REAL, param_end REAL, param_step REAL)")
        self.db_con.commit()


    def write_task(self, task):
        cur = self.db_con.cursor()
        cur.execute("INSERT INTO active_tasks(slurm_id, type, state, x, y, s, "
                    "V, param_start, param_end, param_step, result_db, global_res_db) "
                    "VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
                    task.get_db_record())
        self.db_con.commit()


    def write_task_data(self, task_id, slurm_id, slurm_output_path = None,
                            cp2k_output_fn = None, cp2k_restart_fn = None,
                            wf_filename = None):
        cur = self.db_con.cursor()
        if slurm_output_path is not None:
            slurm_output = self._read_file_if_exists(slurm_output_path)
        else:
            slurm_output = None
        if cp2k_output_fn is not None:
            cp2k_output = self._read_file_if_exists(cp2k_output_fn)
        else:
            cp2k_output = None
        if cp2k_restart_fn is not None:
            cp2k_restart = self._read_file_if_exists(cp2k_restart_fn)
        else:
            cp2k_restart = None
        if wf_filename is not None:
            wf_data = self._read_file_if_exists(wf_filename, is_binary=True)
        else:
            wf_data = None
        cur.execute("INSERT INTO task_data VALUES(?,?,?,?,?,?)",
                    (task_id, slurm_id, slurm_output,
                    cp2k_output, cp2k_restart, wf_data))
        self.db_con.commit()
        print("Wrote task data to database.")


    def move_task_to_completed(self, task_id):
        cur = self.db_con.cursor()
        cur.execute("SELECT type, x, y, s, V, param_start, param_end, param_step "
                    "FROM active_tasks WHERE id=?", (task_id,))
        row = cur.fetchone()
        if row is not None:
            cur.execute("INSERT INTO completed_tasks(type, x, y, s, "
                        "param_start, param_end, param_step) "
                        "VALUES(?,?,?,?,?,?,?)",
                        (row["type"], row["x"], row["y"], row["s"], row["param_start"],
                        row["param_end"], row["param_step"]))
            cur.execute("DELETE FROM active_tasks WHERE id=?",
                        (task_id,))
        self.db_con.commit()


    def update_task_state(self, task_id, new_state):
        cur = self.db_con.cursor()
        cur.execute("UPDATE active_tasks SET state=? WHERE id=?",
                    (new_state, task_id))
        self.db_con.commit()


    def update_task_slurm_id(self, task_id, new_slurm_id):
        cur = self.db_con.cursor()
        cur.execute("UPDATE active_tasks SET slurm_id=? WHERE id=?",
                    (new_slurm_id, task_id))
        self.db_con.commit()


    def update_task(self, task_id, task):
        cur = self.db_con.cursor()
        db_update = task.get_db_update()
        db_update.append(task_id)
        cur.execute("UPDATE active_tasks "
                    "SET state=?, x=?, y=?, s=?, V=? "
                    "WHERE id=?",
                    db_update)
        self.db_con.commit()


    def fetch_active_task_id(self, type_constraint=None, state_constraint=None):
        cur = self.db_con.cursor()
        if type_constraint is not None:
            if state_constraint is not None:
                cur.execute("SELECT id FROM active_tasks "
                            "WHERE type=? AND state=?",
                            (type_constraint, state_constraint))
            else:
                cur.execute("SELECT id FROM active_tasks "
                            "WHERE type=?",
                            (type_constraint,))
        else:
            if state_constraint is not None:
                cur.execute("SELECT id FROM active_tasks "
                            "WHERE state=?",
                            (state_constraint,))
            else:
                cur.execute("SELECT id FROM active_tasks")
        row = cur.fetchone()
        if row is not None:
            return row[0]
        else:
            return None


    def reserve_active_task(self, type_constraint=None, state_constraint=None):
        cur = self.db_con.cursor()
        cur.execute("BEGIN EXCLUSIVE")
        task_id = self.fetch_active_task_id(type_constraint, state_constraint)
        if task_id is not None:
            task = self.extract_task(task_id)
            self.update_task_state(task_id, "reserved")
            return task_id, task
        else:
            self.db_con.commit()
            return None, None


    def extract_task(self, task_id):
        cur = self.db_con.cursor()
        cur.execute("SELECT * FROM active_tasks WHERE id=?",
                    (task_id,))
        row = cur.fetchone()
        if row is not None:
            if row["type"] == "descend_tip":
                return Descend_tip_task(row["x"], row["y"], row["s"], row["V"],
                                            row["param_start"], row["param_end"],
                                            row["param_step"], row["result_db"],
                                            row["global_res_db"], row["state"],
                                            row["slurm_id"])
            elif row["type"] == "tune_bias":
                return Tune_bias_task(row["x"], row["y"], row["s"], row["V"],
                                            row["param_start"], row["param_end"],
                                            row["param_step"], row["result_db"],
                                            row["global_res_db"], row["state"],
                                            row["slurm_id"])
            else:
                return None
        else:
            return None


    def extract_restart_files(self, task_id, cp2k_output_fn, cp2k_restart_fn, wf_filename):
        cur = self.db_con.cursor()
        restart_file_exists = False
        cur.execute("SELECT cp2k_output, cp2k_restart, cp2k_wf_data "
                    "FROM task_data WHERE task_id=?",
                    (task_id,))
        row = cur.fetchone()
        if row is not None:
            if row["cp2k_output"] is not None:
                with open(cp2k_output_fn, 'w') as cp2k_output_file:
                    cp2k_output_file.write(row["cp2k_output"])
            if row["cp2k_restart"] is not None:
                with open(cp2k_restart_fn, 'w') as cp2k_restart_file:
                    cp2k_restart_file.write(row["cp2k_restart"])
                    restart_file_exists = True
            if row["cp2k_wf_data"] is not None:
                with open(wf_filename, 'wb') as cp2k_wf_file:
                    cp2k_wf_file.write(row["cp2k_wf_data"])
        self.db_con.commit()
        cur.execute("DELETE FROM task_data WHERE task_id=?",
                    (task_id,))
        self.db_con.commit()
        return restart_file_exists


    def extract_output(self, task_id):
        cur = self.db_con.cursor()
        cur.execute("SELECT cp2k_output "
                    "FROM task_data WHERE task_id=?",
                    (task_id,))
        row = cur.fetchone()
        return row[0]


    def _read_file_if_exists(self, filename, is_binary=False):
        try:
            if not is_binary:
                with open(filename, 'r') as in_file:
                    contents = in_file.read()
            else:
                with open(filename, 'rb') as in_file: 
                    contents = in_file.read()
            return contents
        except IOError:
            return None
