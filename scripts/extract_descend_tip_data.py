# coding: utf-8
#!/usr/bin/python

import sys
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
#from ase.io.trajectory import Trajectory

from kpfm_sim_result_db import *

from optparse import OptionParser


erm = "Usage: python copy_scan_points_and_wfn_files.py -i from_db_file -o to_db_file (-w <worker_path>  <results_wfn_path> ) (--kpts)"
debug = True
parser = OptionParser()
parser.add_option('-i', '--input_db_file', action ="store", type="string", nargs=1,
                    help='1 file expected <input_db_file>', default = "None")
parser.add_option('-p', '--point', action ="store", type="float", nargs=2,
                    help='1 point = 2 floats expected', default = (float('nan'),float('nan')) )
parser.add_option('-n', '--no_forces', action='store_false', # not mandatory 
                    help="do not get, plot or save the atomic forces", default = True)
parser.add_option('-v', '--visualize', action='store_true', # not mandatory 
                    help="show the calculated curves", default = False)
(options,args) = parser.parse_args()

result_db_file = options.input_db_file
x_tip          = options.point[0]
y_tip          = options.point[1]
AtForces       = options.no_forces
Visu           = options.visualize # normally you do not want to visualize the curves if the amount of points is high #

if (result_db_file == "None") or (x_tip == float('nan') )  or (y_tip == float('nan')):
    sys.exit("Usage: python extract_descend_tip_data.py -i result_db_file -p x_tip y_tip (--no_forces) (--visualize)")

result_db = Result_db(result_db_file)

extract_geometry_traj(result_db, "geometry_x{}_y{}.traj".format(x_tip, y_tip), x_tip, y_tip, 0.0)

ss_1, energies, energy_interp_data, forces_1 = calc_force_curve_from_energy(result_db, x_tip, y_tip, 0.0)
energy_array = np.column_stack((ss_1[::-1], energies[::-1])) # all the data descending now - otherwise it makes the problem,if the scan is not complete !!! #
force_array_1 = np.column_stack((ss_1[::-1], forces_1[::-1]))

if AtForces :
    ss_2, forces_2 = calc_force_curve(result_db, x_tip, y_tip, 0.0)
    force_array_2 = np.column_stack((ss_2[::-1], forces_2[::-1]))

np.savetxt("energy_x{}_y{}.txt".format(x_tip, y_tip), energy_array)
np.savetxt("force_x{}_y{}_dE.txt".format(x_tip, y_tip), force_array_1)
if AtForces :
    np.savetxt("force_x{}_y{}.txt".format(x_tip, y_tip), force_array_2)

#dforces_1 = derivate_force(ss_1[1:-1], forces_1)
#dforces_2 = derivate_force(ss_2, forces_2)

markersize = 2

plt.plot(ss_1, energies, 'bs-', energy_interp_data[:,0], energy_interp_data[:,1], 'ro-', markersize = markersize)
plt.legend(['data points', 'smoothing spline'])
plt.xlabel("Macroscopic tip-sample distance, s (Å)", size=16)
plt.ylabel("Energy (eV)", size=16)
plt.savefig("EvsD_x{}_y{}.png".format(x_tip, y_tip))
if Visu:
    plt.show()
else:
    plt.close()

if AtForces:
    plt.plot(ss_1, 1.0e12*forces_1, 'bs-', ss_2, 1.0e12*forces_2, 'ro-', markersize = markersize)
else: # only derivated forces
    plt.plot(ss_1, 1.0e12*forces_1, 'bs-', markersize = markersize)
plt.legend(['derivative of energy', 'sum of atomic forces'])
plt.xlabel("Macroscopic tip-sample distance, s (Å)", size=16)
plt.ylabel("Force (pN)", size=16)
plt.savefig("FvsD_x{}_y{}.png".format(x_tip, y_tip))
if Visu:
    plt.show()
else:
    plt.close()

#plt.plot(ss_1[2:-2], dforces_1, 'bs-', ss_2[1:-1], dforces_2, 'ro-')
#plt.show()

