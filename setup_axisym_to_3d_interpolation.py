from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
  name = 'Bilinear interpolation from axisymmetric data to 3D grid',
  ext_modules = cythonize("axisym_to_3d_grid_bilinear_interpolation.pyx"),
  include_dirs=[numpy.get_include()]
)
