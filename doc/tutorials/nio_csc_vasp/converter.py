from triqs_dft_tools.converters import VaspConverter
import triqs_dft_tools.converters.plovasp.converter as plo_converter

# Generate and store PLOs
plo_converter.generate_and_output_as_text('plo.cfg', vasp_dir='./')

# run the converter
Converter = VaspConverter(filename = 'vasp', proj_or_hk = 'hk')
Converter.convert_dft_input()
