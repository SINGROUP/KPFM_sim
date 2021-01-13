# coding: utf-8
#!/usr/bin/python

import sys
import sqlite3

m_to_angst = 1.0e10
eps = 1.0e-13


def create_macro_parameters(pot_con, result_con):
    pot_cur = pot_con.cursor()
    result_cur = result_con.cursor()
    pot_cur.execute("SELECT r_tip, t_sample, eps_sample FROM capacitance\
                    GROUP BY r_tip, t_sample, eps_sample")
    for row in pot_cur:
        result_cur.execute("INSERT INTO macro_parameters(r_tip, t_sample, eps_sample)\
                            VALUES(?,?,?)",
                            (row["r_tip"], row["t_sample"], row["eps_sample"]))
    result_con.commit()


def create_macro_scan_points(pot_con, result_con):
    pot_cur = pot_con.cursor()
    result_cur_point = result_con.cursor()
    result_cur_param = result_con.cursor()
    result_cur_param.execute("SELECT * FROM macro_parameters")
    for param_row in result_cur_param:
        pot_cur.execute("SELECT s FROM capacitance\
                        WHERE r_tip=? AND t_sample=? AND eps_sample=?",
                        (param_row["r_tip"], param_row["t_sample"],
                        param_row["eps_sample"]))
        for point_row in pot_cur:
            result_cur_point.execute("INSERT INTO c_scan_point(param_id, s)\
                                    VALUES(?,?)",
                                    (param_row["id"], point_row["s"]*m_to_angst))
        pot_cur.execute("SELECT s FROM potential\
                        WHERE r_tip=? AND t_sample=? AND eps_sample=?\
                        GROUP BY s",
                        (param_row["r_tip"], param_row["t_sample"],
                        param_row["eps_sample"]))
        for point_row in pot_cur:
            result_cur_point.execute("INSERT INTO pot_scan_point(param_id, s)\
                                    VALUES(?,?)",
                                    (param_row["id"], point_row["s"]*m_to_angst))
    result_con.commit()


def copy_capacitance(pot_con, result_con):
    pot_cur = pot_con.cursor()
    result_cur = result_con.cursor()
    result_cur_point = result_con.cursor()
    
    result_cur_point.execute("SELECT * FROM c_scan_point JOIN macro_parameters\
                            ON c_scan_point.param_id=macro_parameters.id")
    for point_row in result_cur_point:
        pot_cur.execute("SELECT value FROM capacitance WHERE r_tip=? AND\
                        t_sample=? AND eps_sample=? AND (s BETWEEN ? AND ?)",
                        (point_row["r_tip"], point_row["t_sample"],
                        point_row["eps_sample"], point_row["s"]/m_to_angst-eps,
                        point_row["s"]/m_to_angst+eps))
        c_row = pot_cur.fetchone()
        result_cur.execute("INSERT INTO capacitance VALUES(?,?)",
                            (point_row["id"], c_row["value"]))
    result_con.commit()


def copy_potential(pot_con, result_con):
    pot_cur = pot_con.cursor()
    result_cur = result_con.cursor()
    result_cur_point = result_con.cursor()
    
    result_cur_point.execute("SELECT * FROM pot_scan_point JOIN macro_parameters\
                            ON pot_scan_point.param_id=macro_parameters.id")
    for point_row in result_cur_point:
        pot_cur.execute("SELECT r, z, value FROM potential WHERE r_tip=? AND\
                        t_sample=? AND eps_sample=? AND (s BETWEEN ? AND ?)",
                        (point_row["r_tip"], point_row["t_sample"],
                        point_row["eps_sample"], point_row["s"]/m_to_angst-eps,
                        point_row["s"]/m_to_angst+eps))
        for pot_row in pot_cur:
            result_cur.execute("INSERT INTO potential VALUES(?,?,?,?)",
                                (point_row["id"], pot_row["r"]*m_to_angst,
                                pot_row["z"]*m_to_angst, pot_row["value"]))
    result_con.commit()


if len(sys.argv) == 3:
    pot_db_file = sys.argv[1]
    result_db_file = sys.argv[2]
else:
    sys.exit("Usage: python copy_pot_to_result_db.py pot_db_file result_db_file")

with sqlite3.connect(result_db_file) as result_db_con:
    result_db_con.row_factory = sqlite3.Row
    cur = result_db_con.cursor()
    
    cur.execute("DROP TABLE IF EXISTS macro_parameters")
    cur.execute("DROP TABLE IF EXISTS c_scan_point")
    cur.execute("DROP TABLE IF EXISTS pot_scan_point")
    cur.execute("DROP TABLE IF EXISTS potential")
    cur.execute("DROP TABLE IF EXISTS capacitance")
    result_db_con.commit()
    cur.execute("CREATE TABLE macro_parameters(id INTEGER PRIMARY KEY,\
                r_tip REAL, t_sample REAL, eps_sample REAL)")
    cur.execute("CREATE TABLE c_scan_point(id INTEGER PRIMARY KEY,\
                param_id INTEGER, s REAL)")
    cur.execute("CREATE TABLE pot_scan_point(id INTEGER PRIMARY KEY,\
                param_id INTEGER, s REAL)")
    cur.execute("CREATE TABLE potential(pot_scan_point_id INTEGER,\
                r REAL, z REAL, value REAL)")
    cur.execute("CREATE TABLE capacitance(c_scan_point_id INTEGER, value REAL)")
    result_db_con.commit()
    
    with sqlite3.connect(pot_db_file) as pot_db_con:
        pot_db_con.row_factory = sqlite3.Row
        create_macro_parameters(pot_db_con, result_db_con)
        print("Macro parameters created")
        create_macro_scan_points(pot_db_con, result_db_con)
        print("Macro scan points created")
        copy_capacitance(pot_db_con, result_db_con)
        print("Capacitance data copied")
        copy_potential(pot_db_con, result_db_con)
        print("Potential data copied")
        
    #cur.execute("ATTACH DATABASE '{}' AS 'pot_db'".format(pot_db_file))
    #cur.execute("INSERT INTO potential SELECT * FROM pot_db.potential")
    #cur.execute("INSERT INTO capacitance SELECT * FROM pot_db.capacitance")
    #result_db_con.commit()
