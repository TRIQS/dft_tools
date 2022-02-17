# Generated automatically using the command :
# c++2py.py -m atm -o atm --moduledoc "Analytical Tetrahedron Method for DOS" ../../../../c++/triqs_dft_tools/converters/vasp/dos_tetra3d.hpp
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "atm", doc = "Analytical Tetrahedron Method for calculating DOS", app_name = "atm")

# All the triqs C++/Python modules

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("triqs_dft_tools/converters/vasp/dos_tetra3d.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/cpp2py_converters/arrays.hpp>
using namespace triqs::arrays;
""")

module.add_function ("array<double,2> dos_tetra_weights_3d (array_view<double,1> eigk, double en, array_view<long,2> itt)", doc = """DOS of a band by analytical tetrahedron method\n\n   Returns corner weights for all tetrahedra for a given band and real energy.""")

module.generate_code()
