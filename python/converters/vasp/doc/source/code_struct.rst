Code Structure
##############

.. toctree::
   
   vaspio
   plotools
   converter

The program consists of three main parts:
 * :doc:`Import of data from VASP files </vaspio>`
 * Processing of projectors according to an input config-file
 * Conversion of the DFT data to TRIQS h5-file

Import of data from VASP files is implemented in `vaspio.py`. The data
is read from VASP files and then stored in objects in raw format
(i.e. practically no processing is done at this stage).
These objects are then combined into a dictionary which can be easily
passed to other routines.

The basic workflow is prescribed as follows:
 * raw data is read from VASP files and passed to the main part
 * in the main part the input config-file is read and the projectors are selected and process accordingly
 * the processed data is stored into output text-files
 * when the TRIQS-converter is requested the data is read from text-files and written into a h5-file in an appropriate format
