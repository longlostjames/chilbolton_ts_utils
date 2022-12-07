API reference
=============


High-level functions
--------------------

CloudnetPy's high-level functions provide a simple mechanism to process
cloud remote sensing measurements into Cloudnet products. A full processing
goes in steps. Each step produces a file which used as an input for the
next step.

Raw data conversion
...................

Different Cloudnet instruments provide raw data in various formats (netCDF, binary, text)
that first need to be converted into homogeneous Cloudnet netCDF files
containing harmonized units and other metadata. This initial processing step
is necessary to ensure that the subsequent processing steps work with
all supported instrument combinations.

.. automodule:: wivern_chilbolton_utils
  :members:

.. autofunction:: wivern_chilbolton_utils.convert_camra_ts_l0a2l0b

.. autofunction:: wivern_chilbolton_utils.convert_camra_ts_l0b2l1

.. autofunction:: wivern_chilbolton_utils.convert_galileo_ts_l0b2l1



Visualizing results
...................

CloudnetPy offers an easy-to-use plotting interface:

.. autofunction:: plotting.generate_figure

There is also possibility to compare CloundetPy files with the
Matlab-processed legacy files
(tagged "legacy" in the `Cloudnet data portal <https://cloudnet.fmi.fi>`_):

.. autofunction:: plotting.compare_files
