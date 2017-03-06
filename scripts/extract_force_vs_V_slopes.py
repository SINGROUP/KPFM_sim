# coding: utf-8
#!/usr/bin/python

import sys
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from ase.io.trajectory import PickleTrajectory

from kpfm_sim_result_db import Result_db

eps = 1.0e-13
eV_to_J = 1.602176565e-19
au_to_eV = 27.211
bohr_to_m = 5.2917721092e-11
au_to_N = au_to_eV*eV_to_J/bohr_to_m

weight_E = 1.0e4
weight_F = 1.0e12
smooth_factor_E = 1.0
smooth_factor_F = 1.0
use_atomic_forces = True

def calc_force_curve_from_energy(result_db, x, y, s_min, s_max, V):
    energies = []
    ss = []
    with result_db:
        scan_points = result_db.get_s_range_scan_points(x, y, s_min, s_max, V)
        for point in scan_points:
            energies.append(result_db.get_energy(point[0])*au_to_eV)
            ss.append(point[1])
    energy_array = np.array(energies)
    s_array = np.array(ss)
    weights = weight_E*np.ones(len(s_array))
    energy_spl = UnivariateSpline(s_array, energy_array, w=weights, k=3, s=smooth_factor_E*len(s_array))
    s_interp = np.linspace(s_array[0], s_array[-1], num=100)
    energy_array_interp = energy_spl(s_interp)
    energy_interp_data = np.column_stack((s_interp, energy_array_interp))
    forces = -eV_to_J*1.0e10*energy_spl(s_array, nu=1)
    return s_array, energy_array, energy_interp_data, forces


def calc_force_curve(result_db, x, y, s_min, s_max, V):
    forces = []
    ss = []
    with result_db:
        scan_points = result_db.get_s_range_scan_points(x, y, s_min, s_max, V)
        tip_top_atoms = result_db.get_model_part_atom_ids("tip", "top")
        tip_center_atoms = result_db.get_model_part_atom_ids("tip", "center")
        tip_apex_atom = result_db.get_model_part_atom_ids("tip", "apex")
        sample_top_atoms = result_db.get_model_part_atom_ids("sample", "top")
        sample_center_atoms = result_db.get_model_part_atom_ids("sample", "center")
        sample_bottom_atoms = result_db.get_model_part_atom_ids("sample", "bottom")
        for point in scan_points:
            scan_point_id = point[0]
            s = point[1]
            atomic_forces = result_db.get_atomic_forces(scan_point_id, sample_top_atoms+sample_center_atoms+sample_bottom_atoms)
            if atomic_forces is not None:
                ss.append(s)
                forces.append(-np.sum(atomic_forces[:,1])*au_to_N)
    s_array = np.array(ss)
    force_array = np.array(forces)
    return s_array, force_array


def smoothen_force_curve(ss, forces):
    weights = weight_F*np.ones(len(ss))
    force_spl = UnivariateSpline(ss, forces, w=weights, k=3, s=smooth_factor_F*len(ss))
    smooth_forces = force_spl(ss, nu=0)
    return ss, smooth_forces


def calc_slopes(Vs, force_array):
    ns = force_array.shape[1]
    slopes = np.zeros(ns)
    slope_variances = np.zeros(ns)
    for i in range(ns):
        if Vs.shape[0] < 2:
            raise Exception("Not enough V data points to make a linear fit (at least 2 needed)")
        elif Vs.shape[0] == 2:
            line_fit = np.polyfit(Vs, force_array[:,i], 1, cov=False)
            slopes[i] = line_fit
        else:
            line_fit, cov = np.polyfit(Vs, force_array[:,i], 2, cov=True)
            slopes[i] = line_fit[1]
            slope_variances[i] = cov[0,0]
    return slopes, slope_variances


def get_V_values(result_db, x, y, s_min, s_max):
    V_list = []
    with result_db:
        s_scan_points = result_db.get_s_range_scan_points(x, y, s_min, s_max, V=None)
        scan_point = s_scan_points[0]
        s = scan_point[1]
        V_scan_points = result_db.get_all_V_scan_points(x, y, s)
    for scan_point in V_scan_points:
        V = scan_point[1]
        V_list.append(V)
    V_array = np.array(V_list)
    return V_array


if len(sys.argv) == 6 or len(sys.argv) == 9:
    result_db_file = sys.argv[1]
    x_tip = float(sys.argv[2])
    y_tip = float(sys.argv[3])
    s_min = float(sys.argv[4])
    s_max = float(sys.argv[5])
    V_min = None
    V_max = None
    V_step = None
    if len(sys.argv) == 9:
        V_min = float(sys.argv[6])
        V_max = float(sys.argv[7])
        V_step = float(sys.argv[8])
else:
    sys.exit("Usage: python {} result_db_file x_tip y_tip s_min s_max [V_min V_max V_step]".format(sys.argv[0]))

result_db = Result_db(result_db_file)

if V_min is None:
    Vs = get_V_values(result_db, x_tip, y_tip, s_min, s_max)
    print Vs
else:
    Vs = np.arange(V_min, V_max+eps, V_step)
    print Vs

forces_list = []
smooth_forces_list = []
forces_grad_list = []
energies_list = []
energy_interp_list = []
for V in Vs:
    ss, energies, energy_interp_data, forces_grad = calc_force_curve_from_energy(result_db, x_tip, y_tip, s_min, s_max, V)
    forces_grad_list.append(forces_grad)
    energies_list.append(energies)
    energy_interp_list.append(energy_interp_data[:,1])
    if use_atomic_forces:
        ss_2, forces = calc_force_curve(result_db, x_tip, y_tip, s_min, s_max, V)
        forces[-1] = 0.0
        forces_list.append(forces)
        ss_3, smooth_forces = smoothen_force_curve(ss_2, forces)
        smooth_forces_list.append(smooth_forces)
force_grad_array = np.array(forces_grad_list)
if use_atomic_forces:
    force_array = np.array(forces_list)
    smooth_force_array = np.array(smooth_forces_list)
energy_array = np.array(energies_list)
energy_interp_array = np.array(energy_interp_list)
ss_interp = energy_interp_data[:,0]

if use_atomic_forces:
    slopes, slope_variances = calc_slopes(Vs, force_array)
    slopes_smooth, slope_variances_smooth = calc_slopes(Vs, smooth_force_array)
slopes_grad, slope_variances_grad = calc_slopes(Vs, force_grad_array)

np.savetxt("slopes_x{}_y{}.txt".format(x_tip, y_tip), np.column_stack((ss, slopes_grad)))

plt.figure()
if use_atomic_forces:
    plt.plot(ss, 1.0e12*slopes, 's-', label="atomic forces")
    plt.plot(ss, 1.0e12*slopes_smooth, 'o-', label="a forces smoothed")
plt.plot(ss, 1.0e12*slopes_grad, 'x-', label="energy gradient")
plt.xlabel(u"Macroscopic tip-sample distance, s (Å)", size=16)
plt.ylabel(u"Force slope (pN/V)", size=16)
plt.legend()

plt.figure(figsize=(10,6))
for i,V in enumerate(Vs):
    plt.plot(ss, energy_array[i,:], 's', label="{:.1f} V".format(V))
plt.gca().set_color_cycle(None)
for i,V in enumerate(Vs):    
    plt.plot(ss_interp, energy_interp_array[i,:], '-', label="{:.1f} V, fit".format(V))
plt.xlabel(u"Macroscopic tip-sample distance, s (Å)", size=16)
plt.ylabel(u"Energy (eV)", size=16)
plt.subplots_adjust(right=0.8)
plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)

plt.figure(figsize=(10,6))
if use_atomic_forces:
    for i,V in enumerate(Vs):
        if abs(V-int(V)) < 1.0e-8:
            plt.plot(ss, 1.0e12*force_array[i,:], 's-', label="{:.1f} V".format(V))
        else:
            plt.plot(ss, 1.0e12*force_array[i,:], 's-', label="{:.1f} V".format(V))
    plt.gca().set_color_cycle(None)
    '''
    for i,V in enumerate(Vs):
        if abs(V-int(V)) < 1.0e-8:
            plt.plot(ss, 1.0e12*smooth_force_array[i,:], 's-', label="{:.1f} V, smooth".format(V))
        else:
            plt.plot(ss, 1.0e12*smooth_force_array[i,:], 's-', label="{:.1f} V, smooth".format(V))
    plt.gca().set_color_cycle(None)
    '''

for i,V in enumerate(Vs):
    if abs(V-int(V)) < 1.0e-8:
        plt.plot(ss, 1.0e12*force_grad_array[i,:], 'o--', label="{:.1f} V, grad".format(V))
    else:
        plt.plot(ss, 1.0e12*force_grad_array[i,:], 'o--', label="{:.1f} V, grad".format(V))

plt.title(u"Force vs. s for different biases", size=16)
plt.xlabel(u"Macroscopic tip-sample distance, s (Å)", size=16)
plt.ylabel(u"Force (pN)", size=16)
plt.subplots_adjust(right=0.75)
plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)

plt.figure(figsize=(10,6))
for i,s in enumerate(ss):
    plt.plot(Vs, 1.0e12*force_grad_array[:,i], 's', label=u"{:.2f} Å".format(s))
plt.gca().set_color_cycle(None)
for i,s in enumerate(ss):
    force_vs_V_fit, force_vs_V_cov = np.polyfit(Vs, force_grad_array[:,i], 2, cov=True)
    force_vs_V_poly = np.poly1d(force_vs_V_fit)
    poly_Vs = np.linspace(Vs[0], Vs[-1], 100)
    force_vs_V_poly_values = force_vs_V_poly(poly_Vs)
    plt.plot(poly_Vs, 1.0e12*force_vs_V_poly_values, '-', label=u"{:.2f} Å, fit".format(s))
plt.title(u"Force vs. V for different s, with quadratic fit", size=16)
plt.xlabel(u"Bias voltage (V)", size=16)
plt.ylabel(u"Force (pN)", size=16)
plt.subplots_adjust(right=0.75)
plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)

plt.show()
