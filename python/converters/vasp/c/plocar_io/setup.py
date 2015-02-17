
from distutils.core import setup, Extension
import numpy

c_plocar_io_mod = Extension('c_plocar_io', sources=['c_plocar_io.c'],
    include_dirs=[numpy.get_include()])

setup(ext_modules=[c_plocar_io_mod])

