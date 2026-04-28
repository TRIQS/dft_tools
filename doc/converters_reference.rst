.. _converters_reference:

DFT Converter Modules
=====================

The DFT converter modules are maintained in the standalone
`dftkit <https://triqs.github.io/dftkit/latest>`_ package.
For the full API reference, see the
`dftkit Python API documentation <https://triqs.github.io/dftkit/latest/_autosummary/triqs_dftkit.html>`_.

The converters can still be imported from ``triqs_dft_tools.converters`` for
backward compatibility:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Converter
     - Backward-compatible import
     - dftkit module
   * - Wien2k
     - ``from triqs_dft_tools.converters import Wien2kConverter``
     - :py:class:`triqs_dftkit.wien2k.Converter`
   * - VASP
     - ``from triqs_dft_tools.converters import VaspConverter``
     - :py:class:`triqs_dftkit.vasp.Converter`
   * - Elk
     - ``from triqs_dft_tools.converters import ElkConverter``
     - :py:class:`triqs_dftkit.elk.Converter`
   * - Wannier90
     - ``from triqs_dft_tools.converters import Wannier90Converter``
     - :py:class:`triqs_dftkit.wannier90.Converter`
   * - H(k)
     - ``from triqs_dft_tools.converters import HkConverter``
     - :py:class:`triqs_dftkit.hk.Converter`
