# coding: utf-8
#!/usr/bin/python

import sys
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from ase.io.trajectory import Trajectory

from kpfm_sim_result_db import Result_db

eV_to_J = 1.602176565e-19
au_to_eV = 27.211
bohr_to_m = 5.2917721092e-11
au_to_N = au_to_eV*eV_to_J/bohr_to_m

smooth_factor = 1.0e-7

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


def calc_force_curve_from_energy(result_db, x, y, V):
    energies = []
    s = []
    with result_db:
        scan_points = result_db.get_all_s_scan_points(x, y, V)
        for point in scan_points:
            energies.append(result_db.get_energy(point[0])*au_to_eV)
            s.append(point[1])
    energies[-1] = energies[-2]
    energy_array = np.array(energies)
    s_array = np.array(s)
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
            print ("debug: atomic_forces", atomic_forces)
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
        for point in scan_points:
            atoms = result_db.extract_atoms_object(point[0])
            traj.write(atoms)
            #write(traj_file, atoms)
    traj.close()


if len(sys.argv) == 4:
    result_db_file = sys.argv[1]
    x_tip = float(sys.argv[2])
    y_tip = float(sys.argv[3])
else:
    sys.exit("Usage: python extract_descend_tip_data.py result_db_file x_tip y_tip")

result_db = Result_db(result_db_file)

extract_geometry_traj(result_db, "geometry_x{}_y{}.traj".format(x_tip, y_tip), x_tip, y_tip, 0.0)

ss_1, energies, energy_interp_data, forces_1 = calc_force_curve_from_energy(result_db, x_tip, y_tip, 0.0)
energy_array = np.column_stack((ss_1, energies))
force_array_1 = np.column_stack((ss_1, forces_1))

ss_2, forces_2 = calc_force_curve(result_db, x_tip, y_tip, 0.0)
force_array_2 = np.column_stack((ss_2, forces_2))

np.savetxt("energy_x{}_y{}.txt".format(x_tip, y_tip), energy_array)
np.savetxt("force_x{}_y{}_dE.txt".format(x_tip, y_tip), force_array_1)
np.savetxt("force_x{}_y{}.txt".format(x_tip, y_tip), force_array_2)

#dforces_1 = derivate_force(ss_1[1:-1], forces_1)
#dforces_2 = derivate_force(ss_2, forces_2)

plt.plot(ss_1, energies, 'bs-', energy_interp_data[:,0], energy_interp_data[:,1], 'ro-')
plt.legend(['data points', 'smoothing spline'])
plt.xlabel("Macroscopic tip-sample distance, s (Å)", size=16)
plt.ylabel("Energy (eV)", size=16)
plt.show()
plt.plot(ss_1, 1.0e12*forces_1, 'bs-', ss_2, 1.0e12*forces_2, 'ro-')
plt.legend(['derivative of energy', 'sum of atomic forces'])
plt.xlabel("Macroscopic tip-sample distance, s (Å)", size=16)
plt.ylabel("Force (pN)", size=16)
plt.show()
#plt.plot(ss_1[2:-2], dforces_1, 'bs-', ss_2[1:-1], dforces_2, 'ro-')
#plt.show()
