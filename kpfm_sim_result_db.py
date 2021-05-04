# -*- coding: utf-8 -*-
#!/usr/bin/python

import shutil, os
import sqlite3
import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms

eps = 1.0e-13
bigeps = 1.0e-6

debug = True

# Class for storing the results of a KPFM simulation to a database and
# extracting the results from there.
class Result_db(object):
    def __init__(self, filename):
        self.db_con = None
        self.db_filename = filename
        try:
            self.db_con = self.open_db()
            self.create_tables()
        except sqlite3.Error as e:
            print("An error occurred:", e.args[0])
        finally:
            self.close_db()
        print("Initialized result database \'{}\'.".format(filename))


    def __del__(self):
        print("Closed result database connection and deleted database object.")
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
        cur.execute("CREATE TABLE IF NOT EXISTS scan_point(id INTEGER PRIMARY KEY, "
                    "x REAL, y REAL, s REAL, V REAL, energy REAL)")
        cur.execute("CREATE TABLE IF NOT EXISTS wf_data(scan_point_id INTEGER, "
                    "wf_data BLOB)") #TODO: is this table needed anymore?
        cur.execute("CREATE TABLE IF NOT EXISTS wf_data_path(scan_point_id INTEGER, "
                    "wf_path TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS output_file(scan_point_id INTEGER, "
                    "output TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS calc_forces_output(scan_point_id INTEGER, "
                    "output TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS model_part(id INTEGER PRIMARY KEY, "
                    "part TEXT, position_in_part TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS atoms(id INTEGER PRIMARY KEY, "
                    "atom_type TEXT, is_fixed BOOLEAN, model_part_id INTEGER)")
        cur.execute("CREATE TABLE IF NOT EXISTS atomic_geometry(scan_point_id INTEGER, "
                    "atom_id INTEGER, x REAL, y REAL, z REAL, mulliken_charge REAL)")
        cur.execute("CREATE TABLE IF NOT EXISTS unit_cell(scan_point_id INTEGER, "
                    "a REAL, b REAL, c REAL, periodicity_id INTEGER)")
        cur.execute("CREATE TABLE IF NOT EXISTS periodicity(id INTEGER PRIMARY KEY, "
                    "periodic_in_x BOOLEAN, periodic_in_y BOOLEAN, periodic_in_z BOOLEAN)")
        cur.execute("CREATE TABLE IF NOT EXISTS atomic_forces(scan_point_id INTEGER, "
                    "atom_id INTEGER, Fx REAL, Fy REAL, Fz REAL)")
        self.db_con.commit()


    # detach is giving error - don't know why
    def copy_atoms_data(self, from_db_file):
        cur = self.db_con.cursor()
        cur.execute("ATTACH ? AS from_db", (from_db_file,))
        cur.execute("INSERT INTO model_part SELECT * FROM from_db.model_part")
        cur.execute("INSERT INTO atoms SELECT * FROM from_db.atoms")
        cur.execute("DETACH from_db")
        self.db_con.commit()


    def get_scan_point_id(self, x, y, s, V):
        cur = self.db_con.cursor()
        cur.execute("SELECT id FROM scan_point WHERE (x BETWEEN ? AND ?) "
                    "AND (y BETWEEN ? AND ?) AND (s BETWEEN ? AND ?) "
                    "AND (V BETWEEN ? AND ?)",
                    (x-eps, x+eps, y-eps, y+eps, s-eps, s+eps, V-bigeps, V+bigeps))
        row = cur.fetchone()
        if row is None:
            return None # change here maybe !
            #return 1
        else:
            return row[0]


    def get_model_part(self):
        cur = self.db_con.cursor()
        cur.execute("SELECT * FROM model_part")
        model_part = []
        pos_in_part = []
        for row in cur:
            model_part.append(row["part"])
            pos_in_part.append(row["position_in_part"])
        if debug:
            print ("model_part_full",model_part)
            print ('pos_in_part',pos_in_part)
        return model_part, pos_in_part;       
        

    def get_all_scan_point_entries(self):
        cur = self.db_con.cursor()
        cur.execute("SELECT * FROM scan_point")
        scan_point_entries = []
        for row in cur:
            scan_point_entries.append((row["id"], row["x"], row["y"],
                                    row["s"], row["V"], row["energy"]))
        return scan_point_entries


    def get_all_s_scan_points(self, x, y, V=0.0):
        cur = self.db_con.cursor()
        cur.execute("SELECT id, s FROM scan_point WHERE (x BETWEEN ? AND ?) "
                    "AND (y BETWEEN ? AND ?) AND (V BETWEEN ? AND ?) ORDER BY s",
                    (x-eps, x+eps, y-eps, y+eps, V-bigeps, V+bigeps))
        scan_points = []
        for row in cur:
            scan_points.append((row["id"], row["s"]))
        return scan_points


    def get_all_V_scan_points(self, x, y, s):
        cur = self.db_con.cursor()
        cur.execute("SELECT id, V FROM scan_point WHERE (x BETWEEN ? AND ?) "
                    "AND (y BETWEEN ? AND ?) AND (s BETWEEN ? AND ?) ORDER BY V",
                    (x-eps, x+eps, y-eps, y+eps, s-eps, s+eps))
        scan_points = []
        for row in cur:
            scan_points.append((row["id"], row["V"]))
        return scan_points


    def get_s_range_scan_points(self, x, y, s_min, s_max, V=0.0):
        cur = self.db_con.cursor()
        if V is not None:
            cur.execute("SELECT id, s FROM scan_point WHERE (x BETWEEN ? AND ?) "
                        "AND (y BETWEEN ? AND ?) AND (V BETWEEN ? AND ?) "
                        "AND (s BETWEEN ? AND ?) ORDER BY s",
                        (x-eps, x+eps, y-eps, y+eps, V-bigeps, V+bigeps, s_min-eps, s_max+eps))
        else:
            cur.execute("SELECT id, s FROM scan_point WHERE (x BETWEEN ? AND ?) "
                        "AND (y BETWEEN ? AND ?) AND (s BETWEEN ? AND ?) "
                        "GROUP BY s ORDER BY s",
                        (x-eps, x+eps, y-eps, y+eps, s_min-eps, s_max+eps))
        scan_points = []
        for row in cur:
            scan_points.append((row["id"], row["s"]))
        return scan_points
    
    
    def get_larger_s_scan_points(self, x, y, s, V=0.0):
        cur = self.db_con.cursor()
        cur.execute("SELECT id, s FROM scan_point WHERE (x BETWEEN ? AND ?) "
                    "AND (y BETWEEN ? AND ?) AND (V BETWEEN ? AND ?) "
                    "AND s > ? ORDER BY s",
                    (x-eps, x+eps, y-eps, y+eps, V-bigeps, V+bigeps, s))
        scan_points = []
        for row in cur:
            scan_points.append((row["id"], row["s"]))
        return scan_points


    def get_neighbor_s_points(self, x, y, s, V):
        cur = self.db_con.cursor()
        
        cur.execute("SELECT id, s FROM scan_point WHERE (x BETWEEN ? AND ?) "
                    "AND (y BETWEEN ? AND ?) AND (V BETWEEN ? AND ?) "
                    "AND s < ? ORDER BY s DESC LIMIT 1",
                    (x-eps, x+eps, y-eps, y+eps, V-eps, V+eps, s-eps))
        row = cur.fetchone()
        if row is None:
            return None, None
        else:
            lower_point = (row["id"], row["s"])
        
        cur.execute("SELECT id, s FROM scan_point WHERE (x BETWEEN ? AND ?) "
                    "AND (y BETWEEN ? AND ?) AND (V BETWEEN ? AND ?) "
                    "AND s > ? ORDER BY s ASC LIMIT 1",
                    (x-eps, x+eps, y-eps, y+eps, V-bigeps, V+bigeps, s+eps))
        row = cur.fetchone()
        if row is None:
            return None, None
        else:
            higher_point = (row["id"], row["s"])
        
        return lower_point, higher_point


    def get_no_forces_scan_points(self):
        cur = self.db_con.cursor()
        cur.execute("SELECT id, s, V FROM scan_point WHERE id NOT IN "
                    "(SELECT scan_point_id FROM atomic_forces) "
                    "ORDER BY s ASC")
        scan_points = []
        for row in cur:
            scan_points.append((row["id"], row["s"], row["V"]))
        return scan_points


    def get_model_part_id(self, model_part, position_in_part):
        cur = self.db_con.cursor()
        cur.execute("SELECT id FROM model_part WHERE part=? AND position_in_part=?",
                    (model_part, position_in_part))
        row = cur.fetchone()
        if row is None:
            return None
        else:
            return row[0]


    def get_model_part_atom_ids(self, model_part, position_in_part=None):
        cur = self.db_con.cursor()
        if position_in_part is None:
            cur.execute("SELECT atoms.id FROM atoms JOIN model_part "
                        "ON atoms.model_part_id=model_part.id "
                        "WHERE part=?",
                        (model_part,))
        else:
            cur.execute("SELECT atoms.id FROM atoms JOIN model_part "
                        "ON atoms.model_part_id=model_part.id "
                        "WHERE part=? AND position_in_part=?",
                        (model_part, position_in_part))
        atom_inds = []
        for row in cur:
            atom_inds.append(row["id"])
        return atom_inds


    def get_sim_parameter(self, param_key):
        cur = self.db_con.cursor()
        cur.execute("SELECT value FROM sim_parameters WHERE key=?",
                    (param_key,))
        row = cur.fetchone()
        if row is None:
            return None
        else:
            return row[0]


    def get_closest_pot_scan_points(self, s):
        cur = self.db_con.cursor()
        cur.execute("SELECT id, s FROM pot_scan_point WHERE s <= ? "
                    "ORDER BY s DESC LIMIT 1", (s,))
        row = cur.fetchone()
        if row is None:
            lower_point = None
        else:
            lower_point = (row["id"], row["s"])
        cur.execute("SELECT id, s FROM pot_scan_point WHERE s > ? "
                    "ORDER BY s ASC LIMIT 1", (s,))
        row = cur.fetchone()
        if row is None:
            higher_point = None
        else:
            higher_point = (row["id"], row["s"])
        return lower_point, higher_point


    def get_external_potential(self, pot_scan_point_id):
        cur = self.db_con.cursor()
        cur.execute("SELECT * FROM potential WHERE pot_scan_point_id=? "
                    "ORDER BY r, z",
                    (pot_scan_point_id,))
        rs = []
        zs = []
        pot = []
        for row in cur:
            rs.append(row["r"])
            zs.append(row["z"])
            pot.append(row["value"])
        cur.execute("SELECT COUNT(DISTINCT z) FROM potential "
                    "WHERE pot_scan_point_id=?",
                    (pot_scan_point_id,))
        row = cur.fetchone()
        nz = row[0]
        r_array = np.array(rs).reshape((-1, nz))
        z_array = np.array(zs).reshape((-1, nz))
        pot_array = np.array(pot).reshape((-1, nz))
        return r_array[:,1], z_array[1,:], pot_array


    def get_energy(self, scan_point_id):
        cur = self.db_con.cursor()
        cur.execute("SELECT energy FROM scan_point WHERE id=?",
                    (scan_point_id,))
        row = cur.fetchone()
        return row[0]


    def get_atomic_forces(self, scan_point_id, atom_ids=None):
        cur = self.db_con.cursor()
        cur.execute("SELECT atom_id, Fx, Fy, Fz FROM atomic_forces "
                    "WHERE scan_point_id=?",
                    (scan_point_id,))
        forces = []
        if atom_ids is None:
            for row in cur:
                forces.append((row["Fx"], row["Fy"], row["Fz"]))
        else:
            for row in cur:
                if row["atom_id"] in atom_ids:
                    forces.append((row["Fx"], row["Fy"], row["Fz"]))
        if forces:
            return np.array(forces)
        else:
            return None


    def get_wf_data_path(self, scan_point_id):
        cur = self.db_con.cursor()
        cur.execute("SELECT wf_path FROM wf_data_path WHERE scan_point_id=?",
                    (scan_point_id,))
        row = cur.fetchone()
        if row is not None:
            return row["wf_path"]
        else:
            return None


    # (x,y,s) are coordinates of the macroscopic tip.
    # V is the effective bias voltage.
    # energy is the potential energy of the tip at the current scan point.
    def write_scan_point(self, x, y, s, V, energy):
        cur = self.db_con.cursor()
        if self.get_scan_point_id(x, y, s, V) is None:
            cur.execute("INSERT INTO scan_point(x, y, s, V, energy) VALUES(?,?,?,?,?)",
                        (x, y, s, V, energy))
            self.db_con.commit()
            return cur.lastrowid
        else:
            #return None # change here!
            return 1


    # wf_filename points to the file where the state information is located.
    def write_wf_data(self, scan_point_id, wf_filename):
        with open(wf_filename, 'rb') as file_in:
            wf_data = file_in.read()
        cur = self.db_con.cursor()
        cur.execute("INSERT INTO wf_data VALUES(?,?)",
                    (scan_point_id, wf_data))
        self.db_con.commit()


    def write_wf_data_path(self, scan_point_id, wf_path):
        cur = self.db_con.cursor()
        cur.execute("INSERT INTO wf_data_path VALUES(?,?)",
                    (scan_point_id, wf_path))
        self.db_con.commit()


    def write_output_file(self, scan_point_id, output):
        cur = self.db_con.cursor()
        cur.execute("INSERT INTO output_file VALUES(?,?)",
                    (scan_point_id, output))
        self.db_con.commit()


    def write_calc_forces_output(self, scan_point_id, output):
        cur = self.db_con.cursor()
        cur.execute("INSERT INTO calc_forces_output VALUES(?,?)",
                    (scan_point_id, output))
        self.db_con.commit()


    def write_model_part(self, model_part, position_in_part):
        cur = self.db_con.cursor()
        cur.execute("INSERT INTO model_part(part, position_in_part) VALUES(?,?)",
                    (model_part, position_in_part))
        self.db_con.commit()
        return cur.lastrowid


    # atoms_model is a ASE Atoms object
    # is_fixed and belongs_to are lists or rank-1 arrays with elements for each atom
    def write_atoms(self, atoms_model, is_fixed, belongs_to, simplistic = False):
        cur = self.db_con.cursor()
        if simplistic:
            model_part_ids = belongs_to;
        else:
            model_part_ids = []
            for atom_belongs_to in belongs_to:
                model_part_id = self.get_model_part_id(atom_belongs_to[0], atom_belongs_to[1])
                if model_part_id is None:
                    model_part_id = self.write_model_part(atom_belongs_to[0], atom_belongs_to[1])
                model_part_ids.append(model_part_id)
        if debug:
            print ("model_part_ids:",model_part_ids)
        for atom_i, atom in enumerate(atoms_model):
            cur.execute("INSERT INTO atoms VALUES(?,?,?,?)",
                        (atom_i, atom.symbol, is_fixed[atom_i], model_part_ids[atom_i]))
        self.db_con.commit()


    # atoms_model is a ASE Atoms object
    # mulliken_charge is a list or rank-1 array with elements for each atom
    def write_atomic_geo(self, scan_point_id, atoms_model, mulliken_charge):
        cur = self.db_con.cursor()
        for atom_i, atom in enumerate(atoms_model):
            pos = atom.position
            cur.execute("INSERT INTO atomic_geometry VALUES(?,?,?,?,?,?)",
                        (scan_point_id, atom_i, pos[0], pos[1], pos[2],
                        mulliken_charge[atom_i]))
        self.db_con.commit()


    # forces is a rank-2 array where the first index identifies the atom and
    # the second one the force component
    def write_atomic_forces(self, scan_point_id, forces):
        cur = self.db_con.cursor()
        for atom_i in range(forces.shape[0]):
            cur.execute("INSERT INTO atomic_forces VALUES(?,?,?,?,?)",
                        (scan_point_id, atom_i, forces[atom_i, 0],
                        forces[atom_i, 1], forces[atom_i, 2]))
        self.db_con.commit()


    # atoms_model is a ASE Atoms object
    def write_unit_cell(self, scan_point_id, atoms_model):
        cell = atoms_model.get_cell()
        pbc = atoms_model.get_pbc()
        cur = self.db_con.cursor()
        cur.execute("SELECT id FROM periodicity WHERE periodic_in_x=? AND periodic_in_y=? AND periodic_in_z=?",
                    (bool(pbc[0]), bool(pbc[1]), bool(pbc[2])))
        row = cur.fetchone()
        if row is None:
            cur.execute("INSERT INTO periodicity(periodic_in_x, periodic_in_y, periodic_in_z) VALUES(?,?,?)",
                        (bool(pbc[0]), bool(pbc[1]), bool(pbc[2])))
            periodicity_id = cur.lastrowid
            self.db_con.commit()
        else:
            periodicity_id = row[0]
        cur.execute("INSERT INTO unit_cell VALUES(?,?,?,?,?)",
                    (scan_point_id, cell[0,0], cell[1,1], cell[2,2], periodicity_id))
        self.db_con.commit()


    def extract_wf_data(self, scan_point_id, wf_filename, project_path=""):
        cur = self.db_con.cursor()
        cur.execute("SELECT wf_data FROM wf_data WHERE scan_point_id=?",
                    (scan_point_id,))
        row = cur.fetchone()
        if row is not None:
            with open(wf_filename, 'wb') as file_out:
                file_out.write(row[0])
            return True
        else:
            cur.execute("SELECT wf_path FROM wf_data_path WHERE scan_point_id=?",
                        (scan_point_id,))
            row = cur.fetchone()
            if row is not None:
                wf_rel_path = row["wf_path"]
                # For backward compatibility
                if os.path.isabs(wf_rel_path):
                    wf_path = wf_rel_path
                    if debug: 
                        print ("debug: row[wf_path]")
                else:
                    wf_path = os.path.join(project_path, wf_rel_path)
                    if debug: 
                        print ("debug: os.path.join(project_path,wf_rel_path)")
                if debug:
                    print("debug: wf_path",wf_path)
                    print("debug: wf_rel_path",wf_rel_path)
                    print("debug: wf_filename",wf_filename)
                try:
                    shutil.copy(wf_path, wf_filename)
                except IOError:
                    print("***\nwf_path entry for scan point {} points to non-existent file \n***".format(scan_point_id))
                    if debug: 
                        print ("debug: extract_wf_data 1")
                    return False
                return True
            else:
                return False


    def extract_output(self, scan_point_id, output_type="geo_opt"):
        cur = self.db_con.cursor()
        if output_type == "geo_opt":
            cur.execute("SELECT output FROM output_file WHERE scan_point_id=?",
                        (scan_point_id,))
        elif output_type == "forces":
            cur.execute("SELECT output FROM calc_forces_output WHERE scan_point_id=?",
                        (scan_point_id,))
        row = cur.fetchone()
        if row is not None:
            return row[0]
        else:
            return None


    def extract_atoms_object(self, scan_point_id, get_charges = False, get_model = False):
        cur = self.db_con.cursor()
        if debug:
           print ("cur:", cur)
           print ("scan_point_id", scan_point_id)
        if scan_point_id is not None:
            cur.execute("SELECT atom_id, x, y, z, atom_type, is_fixed, mulliken_charge, model_part_id "
                        "FROM atomic_geometry "
                        "JOIN atoms ON atomic_geometry.atom_id=atoms.id "
                        "WHERE scan_point_id=? ORDER BY atom_id", (scan_point_id,))
        else:
            cur.execute("SELECT atom_id, x, y, z, atom_type, is_fixed, mulliken_charge, model_part_id "
                        "FROM atomic_geometry "
                        "JOIN atoms ON atomic_geometry.atom_id=atoms.id "
                        "ORDER BY atom_id") #, (scan_point_id,))
        if debug:
            print ("after 1st execution")
        positions = []
        atom_types = []
        fixed_atom_inds = []
        charges = []
        model_part = []
        is_fixed = []
        for row in cur:
            if False:
                print (row["atom_type"])
                print (row["is_fixed"])
                print (row["model_part_id"])
            positions.append((row["x"], row["y"], row["z"]))
            atom_types.append(row["atom_type"])
            if row["is_fixed"]:
                fixed_atom_inds.append(row["atom_id"])
            charges.append(row["mulliken_charge"])
            if get_model:
                model_part.append(row["model_part_id"])
                is_fixed.append(row["is_fixed"])
        if debug:
            print ("after positions and atom_types_appending:")
            print ("positions:", positions)
            print ("atom_types:", atom_types)
            print ("fixed_atom_inds",fixed_atom_inds)
            print ("is_fixed",is_fixed)
            print ('model_part',model_part)
        if atom_types:
            atoms_model = Atoms(symbols=atom_types, positions=positions, charges=charges)
            if debug:
                 print ("atoms_model:", atoms_model)
        else:
            if debug:
                 print ("it is going to return None")
            return None
        
        if debug:
            print ("before 2nd execution")
        if scan_point_id is not None:
            cur.execute("SELECT a, b, c, periodic_in_x, periodic_in_y, periodic_in_z\
                         FROM unit_cell\
                         JOIN periodicity\
                         ON unit_cell.periodicity_id=periodicity.id\
                         WHERE scan_point_id=?", (scan_point_id,))
        else:
            cur.execute("SELECT a, b, c, periodic_in_x, periodic_in_y, periodic_in_z\
                         FROM unit_cell\
                         JOIN periodicity\
                         ON unit_cell.periodicity_id=periodicity.id")
           
        if debug:
            print ("after 2nd execution")
        row = cur.fetchone()
        atoms_model.set_cell([row["a"], row["b"], row["c"]])
        atoms_model.set_pbc([row["periodic_in_x"], row["periodic_in_y"],
                                row["periodic_in_z"]])
        fix_bulk = FixAtoms(fixed_atom_inds)
        atoms_model.set_constraint(fix_bulk)
        
        if debug:
            print ("atoms_model:", atoms_model)
            print ("atoms_model.get_global_number_of_atoms():", atoms_model.get_global_number_of_atoms())
        if get_model:
            return (atoms_model, model_part, is_fixed), charges
        if get_charges:
            return atoms_model, charges
        else:
            return atoms_model


    def delete_scan_point(self, scan_point_id):
        cur = self.db_con.cursor()
        cur.execute("DELETE FROM scan_point WHERE id=?",
                    (scan_point_id,))
        cur.execute("DELETE FROM wf_data WHERE scan_point_id=?",
                    (scan_point_id,))
        cur.execute("DELETE FROM output_file WHERE scan_point_id=?",
                    (scan_point_id,))
        cur.execute("DELETE FROM calc_forces_output WHERE scan_point_id=?",
                    (scan_point_id,))
        cur.execute("DELETE FROM atomic_geometry WHERE scan_point_id=?",
                    (scan_point_id,))
        cur.execute("DELETE FROM atomic_forces WHERE scan_point_id=?",
                    (scan_point_id,))
        cur.execute("DELETE FROM unit_cell WHERE scan_point_id=?",
                    (scan_point_id,))
        cur.execute("SELECT * FROM wf_data_path WHERE scan_point_id=?",
                    (scan_point_id,))
        row = cur.fetchone()
        if row is not None:
            if debug:
                print("debug: wf_path",wf_path)
            try:
                os.remove(row["wf_path"])
            except OSError:
                print("***\nwf_path entry for scan point {} pointed to non-existent file\n***".format(scan_point_id))
                if debug:
                    print("debug: delete_scan_point 2")
            cur.execute("DELETE FROM wf_data_path WHERE scan_point_id=?",
                        (scan_point_id,))
        self.db_con.commit()


    def delete_wf_data_path(self, scan_point_id):
        cur = self.db_con.cursor()
        cur.execute("DELETE FROM wf_data_path WHERE scan_point_id=?",
                    (scan_point_id,))
        self.db_con.commit()
