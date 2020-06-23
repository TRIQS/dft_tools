.. _refconverters:

Converters
==========

Wien2k Converter
----------------
.. autoclass:: triqs_dft_tools.converters.wien2k.Wien2kConverter
   :members:
   :special-members:
   :show-inheritance:

H(k) Converter
--------------
.. autoclass:: triqs_dft_tools.converters.hk.HkConverter
   :members:
   :special-members:

Wannier90 Converter
-------------------
.. autoclass:: triqs_dft_tools.converters.wannier90.Wannier90Converter
   :members:
   :special-members:

PLOVASP
----------
.. _refPLOVASP:

PLOVASP reference, the classes / functions are sorted the way the converter uses them.

.. automodule:: triqs_dft_tools.converters.plovasp.converter
  :members: generate_and_output_as_text

.. automodule:: triqs_dft_tools.converters.plovasp.inpconf
  :members: ConfigParameters

.. automodule:: triqs_dft_tools.converters.plovasp.vaspio
  :members: VaspData, Plocar, Poscar, Kpoints, Eigenval, Doscar, read_symmcar

.. automodule:: triqs_dft_tools.converters.plovasp.elstruct
  :members: ElectronicStructure

.. automodule:: triqs_dft_tools.converters.plovasp.plotools
  :members:

.. automodule:: triqs_dft_tools.converters.plovasp.proj_shell
  :members:

.. automodule:: triqs_dft_tools.converters.plovasp.proj_group
  :members:


VASP Converter
-------------------
.. _refVASPconverter:
.. autoclass:: triqs_dft_tools.converters.vasp.VaspConverter
   :members:
   :special-members:


Converter Tools
---------------
.. autoclass:: triqs_dft_tools.converters.converter_tools.ConverterTools
   :members:
   :special-members:
