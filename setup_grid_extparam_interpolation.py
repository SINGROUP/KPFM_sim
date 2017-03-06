from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'Linear interpolation of data on grid with respect to an external parameter',
  ext_modules = cythonize("grid_linear_interpolation_extparam.pyx"),
)
