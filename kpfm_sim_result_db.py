# -*- coding: utf-8 -*-
#!/usr/bin/python

import shutil, os
import sqlite3
import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io.trajectory import Trajectory

from scipy.interpolate import UnivariateSpline # only for post_processing6;15M6;15m

# important constants: #

eV_to_J = 1.602176565e-19
au_to_eV = 27.211
bohr_to_m = 5.2917721092e-11
au_to_N = au_to_eV*eV_to_J/bohr_to_m

smooth_factor = 1.0e-7

eps = 1.0e-13
bigeps = 1.0e-6

nd = 6 # number of digits to round

debug = False

# -- at the end other functions used for handling db when copying, or extracting and post-processing data ... --- #

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
        cur.execute("CREATE TABLE IF NOT EXISTS pot_data_path(scan_point_id INTEGER, "
                    "pot_path TEXT)")
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

    def get_pot_data_path(self, scan_point_id, point=None, V= 0. ):
        cur = self.db_con.cursor()
        cur.execute("SELECT pot_path FROM pot_data_path WHERE scan_point_id=?",
                    (scan_point_id,))
        row = cur.fetchone()
        if row is not None:
            return row["pot_path"]
        else:
            print("D: path not found for this")
            print("D: V, point",V, point)
            if abs(V) > eps:
                now_id_tmp = self.get_scan_point_id(point[0], point[1], point[2], 0.0)
                print("D: now_id_tmp",now_id_tmp)
                cur = self.db_con.cursor()
                cur.execute("SELECT pot_path FROM pot_data_path WHERE scan_point_id=?",
                    (now_id_tmp,))
                row = cur.fetchone()
                print ("D: row",row)
                if row is not None:
                    return row["pot_path"]
                else:
                    print("STILL NO PATH FOUND")
                    return None
            else:
                print("!!! Path not found at all, or something went wrong!!!")
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

    def write_pot_data_path(self, scan_point_id, pot_path):
        cur = self.db_con.cursor()
        cur.execute("INSERT INTO pot_data_path VALUES(?,?)",
                    (scan_point_id, pot_path))
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
            print ("debug: model_part_ids:",model_part_ids)
            print ("debug: len(is_fixed)",len(is_fixed))
            print ("debug: len(model_part_ids)",len(model_part_ids))
        for atom_i, atom in enumerate(atoms_model):
            cur.execute("INSERT INTO atoms VALUES(?,?,?,?)",
                        (atom_i, atom.symbol, is_fixed[atom_i], model_part_ids[atom_i]))
        self.db_con.commit()


    # atoms_model is a ASE Atoms object
    # mulliken_charge is a list or rank-1 array with elements for each atom
    def write_atomic_geo(self, scan_point_id, atoms_model, mulliken_charge):
        cur = self.db_con.cursor()
        if mulliken_charge is not None:
            for atom_i, atom in enumerate(atoms_model):
                pos = atom.position
                cur.execute("INSERT INTO atomic_geometry VALUES(?,?,?,?,?,?)",
                            (scan_point_id, atom_i, pos[0], pos[1], pos[2],
                            mulliken_charge[atom_i]))
        else: # charges = None #
            for atom_i, atom in enumerate(atoms_model):
                pos = atom.position
                tmp_val = None
                cur.execute("INSERT INTO atomic_geometry VALUES(?,?,?,?,?,?)",
                            (scan_point_id, atom_i, pos[0], pos[1], pos[2],
                            tmp_val)) # putting there NULL
        
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
        if debug:
            print ("D atoms_model",atoms_model)
            print ("D fixed_atoms_inds",fixed_atom_inds)
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

# ******************************************************************************************************************* #
# *************   Now the functions for database handling, when something  .... ************************************* #
# ******************************************************************************************************************* #
#

def prepare_db_for_task(global_res_db_file, result_db_file, task_db_file):
    '''
    prepare_db_for_task(global_res_db_file, result_db_file, task_db_file)
    adjust the result db file and the task db file with the all/last results from the global results file
    only the important parts (scan points and geometry) copied
    '''
    from_db = Result_db(global_res_db_file)
    to_db = Result_db(result_db_file)
    control_db = Result_db(task_db_file)
    gm = True
    with from_db:
        scan_points = from_db.get_all_scan_point_entries()
        #if scan_points == None:
        #    print("The global input file was not found; \n DEBUG: from_db",from_db,"\n DEBUG: global_res_db_file",global_res_db_file, "\n" )
        with to_db :
            with control_db :
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
                            global_res_db_file, to_id, result_db_file))
                        tmp, charges = from_db.extract_atoms_object(from_id, get_charges=True, get_model=gm);
                        pot_data_path = from_db.get_pot_data_path(from_id)
                        if gm:
                            atoms = tmp[0]
                            model_part = tmp[1]
                            is_fixed = tmp[2]
                            full_model, pos_in_part  = from_db.get_model_part()
                            if debug:
                                print ("debug: model_part",model_part)
                                print ('debug: is_fixed', is_fixed)
                                print ('debug: full_model',full_model);
                                print ('debug: pos_in_part',pos_in_part);
                                print ('debug: pot_data_path',pot_data_path)
                            to_db.write_atoms(atoms, is_fixed, model_part, simplistic = True)
                            for i in range(len(full_model)):
                                to_db.write_model_part(full_model[i],pos_in_part[i])
                            gm = False;
                        else:
                            atoms = tmp
                        #print("debug: atoms",atoms)
                        #print("debug: atoms.positions", atoms.positions )
                        to_db.write_atomic_geo(to_id, atoms, charges)
                        to_db.write_unit_cell(to_id, atoms)
                        if pot_data_path is not None:
                            to_db.write_pot_data_path(to_id,pot_data_path)
                    ####
                    if control_db.get_scan_point_id(x, y, s, V) is None:
                        control_id = control_db.write_scan_point(x, y, s, V, energy)
                        print("Copying scan point {} from {} to scan point {} in {}".format(from_id,
                            global_res_db_file, control_id, task_db_file))
    print()
    print("results and tasks db files updated")

def prepare_db_for_small(global_res_db_file, result_db_file,pxy=None):
    '''
    prepare_db_for_small(global_res_db_file, result_db_file,pxy=None)
    adjust the result db file with the all/last results from the global results file
    only the important parts (scan points and geometry) copied (if not created before)
    pxy - None or 2 floats /x/ and /y/ for points for which you want to calculate for only one /z/ scan
    also - it gives total numbers of calculated geometries, and numbers of relevant() geometries that already does not have a potential
    !!! ASUMPTION: there are the same 
    '''
    from_db = Result_db(global_res_db_file)
    to_db = Result_db(result_db_file)
    gm = True # get_full_model #
    ids=[]; pot_ids1=[]; pot_ids2=[] # ids - all the (relevant) IDS for which one can calculate potentials ; pot_ids - IDS for which the potential is calculated #
    with from_db:
        scan_points = from_db.get_all_scan_point_entries()
        with to_db :
            for scan_point in scan_points:
                from_id = scan_point[0]
                x = scan_point[1]
                y = scan_point[2]
                s = scan_point[3]
                V = scan_point[4]
                energy = scan_point[5]
                if abs(V) < eps:
                    if pxy is None:
                        ids.append(from_id)
                    else: # pxy = [x,y]
                        if (abs(x - pxy[0])<eps) and (abs(y - pxy[1])<eps):
                            ids.append(from_id)
                to_scan_point = to_db.get_scan_point_id(x, y, s, V)
                if to_scan_point is None: #to_db.get_scan_point_id(x, y, s, V) is None:
                    to_id = to_db.write_scan_point(x, y, s, V, energy)
                    print("Copying scan point {} from {} to scan point {} in {}".format(from_id,
                        global_res_db_file, to_id, result_db_file))
                    gm = gm if to_id < 2 else False ; #not to get and write full model, if you are copying to already partially written file#
                    # print ("debug: gm",gm)
                    tmp, charges = from_db.extract_atoms_object(from_id, get_charges=True, get_model=gm);
                    pot_data_path = from_db.get_pot_data_path(from_id)
                    if gm:
                        atoms = tmp[0]
                        model_part = tmp[1]
                        is_fixed = tmp[2]
                        full_model, pos_in_part  = from_db.get_model_part()
                        if debug:
                            print ("debug: model_part",model_part)
                            print ('debug: is_fixed', is_fixed)
                            print ('debug: full_model',full_model);
                            print ('debug: pos_in_part',pos_in_part);
                            print ('debug: pot_data_path',pot_data_path)
                        to_db.write_atoms(atoms, is_fixed, model_part, simplistic = True)
                        for i in range(len(full_model)):
                            to_db.write_model_part(full_model[i],pos_in_part[i])
                        gm = False;
                    else:
                        atoms = tmp
                    to_db.write_atomic_geo(to_id, atoms, charges)
                    to_db.write_unit_cell(to_id, atoms)
                    if pot_data_path is not None:
                        to_db.write_pot_data_path(to_id,pot_data_path)
                        pot_ids1.append(from_id) 
                else: #->to_db.get_scan_point_id(x,y,s,V) is not None:
                    # print ("debug: to_scan_point",to_scan_point) ; exit()
                    if to_db.get_pot_data_path(to_scan_point) is not None :
                        pot_ids2.append(from_id)
                        print ("potential for ",to_scan_point," is already calculated")
                    ####
    print()
    print("result db file updated; calculating the relevant non-calculated ids")
    ids = np.array(ids, dtype = np.int32)
    pot_ids = np.union1d(np.array(pot_ids1, dtype = np.int32),np.array(pot_ids2, dtype = np.int32))
    mask = np.isin(ids,pot_ids, invert=True); no_pot_ids = ids[mask]
    if True:
        print("debug:  no_pot_id", no_pot_ids)
    return no_pot_ids ; # list of numbers of (relevant) geometries, that doesn't have potential calculated #

def copy_db_ft(from_db_file, to_db_file, w1_path=None, wo_path=None, wf_ext=None):
    '''
    copy_db_ft()
    '''
    #w1_path = w1_path if w1_path is not "None" else None
    #wo_path = wo_path if wo_path is not "None" else None
    if debug:
        print ("debug: w1_path, wo_path",w1_path, wo_path)
    from_db = Result_db(from_db_file)
    to_db = Result_db(to_db_file)
    gm = True ## - just to see if the model and fixed part are supposed to be copied #
    bwfc = True; ecl = [] # # for error message, if problem with wfn copying ##
    with from_db:
        scan_points = from_db.get_all_scan_point_entries()
        #if scan_points == None:
        #    print("The global input file was not found; \n DEBUG: from_db",from_db,"\n DEBUG: global_res_db_file",global_res_db_file, "\n" )
        with to_db :
            for scan_point in scan_points:
                from_id = scan_point[0]
                x = scan_point[1]
                y = scan_point[2]
                s = scan_point[3]
                V = scan_point[4]; Vtmp = V if V < 1000.0 else 1000.0 ; #just to adjust reruns with 1000, when the number grows with every new /s/ point #
                energy = scan_point[5]
                if to_db.get_scan_point_id(x, y, s, Vtmp) is None:
                    if debug:
                        print ("D: V, Vtmp", V , Vtmp)
                    to_id = to_db.write_scan_point(x, y, s, Vtmp, energy)
                    print("Copying scan point {} from {} to scan point {} in {}".format(from_id,
                            from_db_file, to_id, to_db_file))
                    gm_tmp = True if (gm and from_id < 2) else False
                    tmp, charges = from_db.extract_atoms_object(from_id, get_charges=True, get_model=gm_tmp);
                    forces = from_db.get_atomic_forces(from_id)
                    output = from_db.extract_output(from_id)
                    calc_forces_output = from_db.extract_output(from_id, "forces")
                    wf_path = from_db.get_wf_data_path(from_id)
                    pot_data_path = from_db.get_pot_data_path(from_id)
                    if gm_tmp :
                        atoms = tmp[0]
                        model_part = tmp[1]
                        is_fixed = tmp[2]
                        full_model, pos_in_part  = from_db.get_model_part()
                        if debug:
                            print ("debug: model_part",model_part)
                            print ('debug: is_fixed', is_fixed)
                            print ('debug: full_model',full_model);
                            print ('debug: pos_in_part',pos_in_part);
                            print ('debug: pot_data_path',pot_data_path)
                        to_db.write_atoms(atoms, is_fixed, model_part, simplistic = True)
                        for i in range(len(full_model)):
                            to_db.write_model_part(full_model[i],pos_in_part[i])
                        gm = False;
                    else:
                        atoms = tmp
                    if debug:
                        print("debug: atoms",atoms) ; print("debug: atoms.positions", atoms.positions )
                    to_db.write_atomic_geo(to_id, atoms, charges)
                    to_db.write_unit_cell(to_id, atoms)
                    to_db.write_output_file(to_id, output)
                    if pot_data_path is not None:
                        to_db.write_pot_data_path(to_id,pot_data_path)
                    if (wf_path is not None) and (w1_path is not None) and (wo_path is not None) and (wf_ext is not None):
                        wf_path1 = wo_pth +"/scan_point-"+str(to_id)+"-RESTART."+wf_ext
                        f_path   = w1_pth +"/"+ wf_path
                        if debug:
                            print("wf_ext",wf_ext)
                            print("wf_path",wf_path)
                            print("wf_path1",wf_path1)
                            print("f_path",f_path)
                        try:
                            shu.copyfile(f_path,wf_path1)
                            print ("the wave-function file:",f_path,"was coppied to",wf_path1)
                        except:
                            bwfc = False
                            tmp  = "PROBLEM: CANNOT COPY - the wave-function file: "+f_path+" cannot be coppied to: "+wf_path1
                            print ( tmp )
                            ecl.append( tmp )
                        to_db.write_wf_data_path(to_id, wf_path1)
    print()
    print("results and tasks db files updated")
    #
    return bwfc, ecl ; # if there are some problems during wf files copying #

# --- functions for post-processing : --- #
def force_from_energy_fd(s, energies):
    forces = []
    for i in range(1, len(s)-1):
        force = eV_to_J*1.0e10*(energies[i-1]-energies[i+1])/(s[i+1]-s[i-1])
        forces.append(force)
    return forces


def derivate_force(s, forces):
    dforces = []
    for i in range(1, len(s)-1):
        dforce = (forces[i-1]-forces[i+1])/(s[i+1]-s[i-1])
        dforces.append(dforce)
    return dforces


def calc_force_curve_from_energy(result_db, x, y, V,  ignore_first = False):
    # debug = True
    energies = []
    s = []
    with result_db:
        scan_points = result_db.get_all_s_scan_points(x, y, V)
        for point in scan_points:
            energies.append(result_db.get_energy(point[0])*au_to_eV)
            s.append(round(point[1],nd))
    if ignore_first:
        s=s[:-1]
        energies=energies[:-1]
    if debug:
        print("D: energies,l",energies,len(energies))
    energies[-1] = energies[-2]
    energy_array = np.array(energies)
    s_array = np.array(s)
    if debug:
        print("D: energy_array,l",energy_array,len(energy_array))
        print("D: s_array,l",s_array,len(s_array))
    energy_spl = UnivariateSpline(s_array, energy_array, k=3, s=smooth_factor*len(s_array))
    s_interp = np.linspace(s_array[0], s_array[-1], num=100)
    energy_array_interp = energy_spl(s_interp)
    energy_interp_data = np.column_stack((s_interp, energy_array_interp))
    forces = -eV_to_J*1.0e10*energy_spl(s_array, nu=1)
    return s_array, energy_array, energy_interp_data, forces


def calc_force_curve_from_energy_fd(result_db, x, y, V):
    energies = []
    s = []
    with result_db:
        scan_points = result_db.get_all_s_scan_points(x, y, V)
        for point in scan_points:
            energies.append(result_db.get_energy(point[0])*au_to_eV)
            s.append(point[1])
    forces = force_from_energy_fd(s, energies)
    return s, energies, forces


def calc_force_curve(result_db, x, y, V):
    forces = []
    ss = []
    with result_db:
        scan_points = result_db.get_all_s_scan_points(x, y, V)
        tip_top_atoms = result_db.get_model_part_atom_ids("tip", "top")
        tip_center_atoms = result_db.get_model_part_atom_ids("tip", "center")
        tip_apex_atom = result_db.get_model_part_atom_ids("tip", "apex")
        sample_top_atoms = result_db.get_model_part_atom_ids("sample", "top")
        sample_center_atoms = result_db.get_model_part_atom_ids("sample", "center")
        sample_bottom_atoms = result_db.get_model_part_atom_ids("sample", "bottom")
        for point in scan_points:
            scan_point_id = point[0]
            s = point[1]
            atomic_forces = result_db.get_atomic_forces(scan_point_id, tip_top_atoms)
            if debug:
                print ("debug: atomic_forces", atomic_forces)
                print ("debug: len(atomic_forces)", scan_point_id, len(atomic_forces))
            if atomic_forces is not None:
                ss.append(s)
                forces.append(np.sum(atomic_forces[:,1])*au_to_N)
    s_array = np.array(ss)
    force_array = np.array(forces)
    return s_array, force_array


def extract_geometry_traj(result_db, traj_file, x, y, V):
    traj = Trajectory(traj_file, "w")
    with result_db:
        scan_points = result_db.get_all_s_scan_points(x, y, V)
        latoms = []
        for point in scan_points:
            atoms = result_db.extract_atoms_object(point[0])
            latoms.append(atoms)
            #traj.write(atoms)
            #write(traj_file, atoms)
        for atoms in latoms[::-1]:
            traj.write(atoms)
    del latoms;
    traj.close()

def extract_charges_descent(result_db,  x, y, V):
    #traj = Trajectory(traj_file, "w")
    with result_db:
        scan_points = result_db.get_all_s_scan_points(x, y, V)
<<<<<<< HEAD
        print("D scan_points",scan_points)
=======
        if debug:
            print("D scan_points",scan_points)
>>>>>>> 3038211cc4f31f47a2c47a9131218cd53fc6ef61
        lcharges = []
        for point in scan_points:
            atoms, charges  = result_db.extract_atoms_object(point[0], get_charges = True)
            lcharges.append(charges)
        lcharges = np.array(lcharges[::-1])
        if debug:
            print ("D full charges.shape",lcharges.shape)
    return lcharges;

# The END ??? #

