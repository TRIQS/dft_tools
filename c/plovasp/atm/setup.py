
from distutils.core import setup, Extension
import numpy

c_atm_dos_mod = Extension('c_atm_dos', sources=['dos_tetra3d.c', 'argsort.c'],
    include_dirs=[numpy.get_include()])

setup(ext_modules=[c_atm_dos_mod])

