# coding: utf-8
#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from ase.io.trajectory import PickleTrajectory

from kpfm_sim_result_db import Result_db

eV_to_J = 1.602176565e-19
au_to_eV = 27.211
bohr_to_m = 5.2917721092e-11
au_to_N = au_to_eV*eV_to_J/bohr_to_m


def force_from_energy(s, energies):
    force = eV_to_J*1.0e10*(energies[0]-energies[2])/(s[2]-s[0])
    return force


def calc_force_curve_from_energy(result_db, x, y, s):
    forces = []
    Vs = []
    with result_db:
        scan_points = result_db.get_all_V_scan_points(x, y, s)
        for point in scan_points:
            V = point[1]
            lower_point, higher_point = result_db.get_neighbor_s_points(x, y, s, V)
            if lower_point is None:
                continue
            ss = [lower_point[1], s, higher_point[1]]
            energies = [result_db.get_energy(lower_point[0])*au_to_eV,
                        result_db.get_energy(point[0])*au_to_eV,
                        result_db.get_energy(higher_point[0])*au_to_eV]
            forces.append(force_from_energy(ss, energies))
            Vs.append(V)
    return Vs, forces


def calc_force_curve(result_db, x, y, s):
    forces = []
    Vs = []
    with result_db:
        scan_points = result_db.get_all_V_scan_points(x, y, s)
        tip_top_atoms = result_db.get_model_part_atom_ids("tip", "top")
        for point in scan_points:
            scan_point_id = point[0]
            V = point[1]
            atomic_forces = result_db.get_atomic_forces(scan_point_id, tip_top_atoms)
            if atomic_forces is not None:
                Vs.append(V)
                forces.append(np.sum(atomic_forces[:,1])*au_to_N)
    return Vs, forces


def extract_energies(result_db, x, y, s):
    energies = []
    Vs = []
    with result_db:
        scan_points = result_db.get_all_V_scan_points(x, y, s)
        for point in scan_points:
            Vs.append(point[1])
            energies.append(result_db.get_energy(point[0])*au_to_eV)
    return Vs, energies


def extract_geometry_traj(result_db, traj_file, x, y, s):
    traj = PickleTrajectory(traj_file, "w")
    with result_db:
        scan_points = result_db.get_all_V_scan_points(x, y, s)
        for point in scan_points:
            atoms = result_db.extract_atoms_object(point[0])
            traj.write(atoms)
    traj.close()


if len(sys.argv) == 5:
    result_db_file = sys.argv[1]
    x_tip = float(sys.argv[2])
    y_tip = float(sys.argv[3])
    s_tip = float(sys.argv[4])
else:
    sys.exit("Usage: python extract_tune_bias_data.py <result_db_file> <x_tip> <y_tip> <s_tip>")

result_db = Result_db(result_db_file)

extract_geometry_traj(result_db, "geometry_x{}_y{}_s{}.traj".format(x_tip, y_tip, s_tip),
                        x_tip, y_tip, s_tip)

Vs_1, energies = extract_energies(result_db, x_tip, y_tip, s_tip)
energy_array = np.array([Vs_1, energies])

Vs_1, forces_1 = calc_force_curve_from_energy(result_db, x_tip, y_tip, s_tip)
force_array_1 = np.array([Vs_1, forces_1])

Vs_2, forces_2 = calc_force_curve(result_db, x_tip, y_tip, s_tip)
force_array_2 = np.array([Vs_2, forces_2])

np.savetxt("energy_x{}_y{}_s{}.txt".format(x_tip, y_tip, s_tip), np.transpose(energy_array))
np.savetxt("force_x{}_y{}_s{}_dE.txt".format(x_tip, y_tip, s_tip), np.transpose(force_array_1))
np.savetxt("force_x{}_y{}_s{}.txt".format(x_tip, y_tip, s_tip), np.transpose(force_array_2))

line_fit = np.polyfit(Vs_1, forces_1, 1)
line = np.poly1d(line_fit)
force_line = line(Vs_1)
plt.plot(force_array_1[0,:], 1.0e12*force_array_1[1,:], 'rs', markersize=8.0)
plt.plot(Vs_1, 1.0e12*force_line, 'b-', linewidth=2.0)
plt.xlabel("Bias voltage (V)", size=16)
plt.ylabel("Force (pN)", size=16)
plt.legend(["data points", "linear fit"])
plt.show()
