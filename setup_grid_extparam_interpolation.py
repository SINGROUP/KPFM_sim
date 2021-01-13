from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
  name = 'Linear interpolation of data on grid with respect to an external parameter',
  ext_modules = cythonize("grid_linear_interpolation_extparam.pyx"),
  include_dirs=[numpy.get_include()]
)
