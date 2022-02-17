# coding: utf-8
#!/usr/bin/python

import sys
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
#from ase.io.trajectory import Trajectory

from kpfm_sim_result_db import *

from optparse import OptionParser


erm = "Usage: python extract_descent_tip_data.py -i result_db_file -p x_tip y_tip (--no_forces) (--visualize) (--voltage V)"
debug = True
parser = OptionParser()
parser.add_option('-i', '--input_db_file', action ="store", type="string", nargs=1,
                    help='1 file expected <input_db_file>', default = "None")
parser.add_option('-p', '--point', action ="store", type="float", nargs=2,
                    help='1 point = 2 floats expected', default = (float('nan'),float('nan')) )
parser.add_option('-V', '--voltage', action ="store", type="float", nargs=1,
                    help='optional: voltage 1 float expected -- default 0.0', default = 0.0 )
(options,args) = parser.parse_args()

result_db_file = options.input_db_file
x_tip          = options.point[0]
y_tip          = options.point[1]
V              = options.voltage

if (result_db_file == "None") or (x_tip == float('nan') )  or (y_tip == float('nan')):
    sys.exit("Usage: python extract_descend_tip_data.py -i result_db_file -p x_tip y_tip (--no_forces) (--visualize) (--voltage V)")

result_db = Result_db(result_db_file)

lcharges = extract_charges_descent(result_db,  x_tip, y_tip, V)

fname = "charges_x{}_y{}_V{}.txt".format(x_tip,y_tip,V)

if debug:
    print ("D: fname",fname)
np.savetxt(fname,lcharges)

