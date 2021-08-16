from triqs_dft_tools.converters.vasp import *
Converter = VaspConverter(filename = 'nio', proj_or_hk = 'hk')
Converter.convert_dft_input()
