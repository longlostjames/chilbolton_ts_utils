========
Overview
========

WIVERN-2 Chilbolton Campaign Dataset
------------------------------------

This documentation describes the dataset collected using ground-based radars
at Chilbolton Obervatory in the UK for the Doppler Wind Radar Science
Performance Study (WIVERN-2), funded by the European Space Agency (ESA) under
Contract no. 4000130864/20/NL/CT.

Level 1 calibrated data in NetCDF-4 format (see :ref:`file-format`) are provided
from the following three radars:

* 3GHz CAMRa radar (ncas-radar-camra-1)
* 35GHz Copernicus radar (ncas-radar-ka-band-1)
* 94GHz Galileo radar (ncas-radar-w-band-1)

The data are arranged into a series of one or more rain events per day of
observation, with any extended non-rainy periods being omitted.  The data may be
read by a variety of standard packages including MATLAB, IDL and R, as well as by
tailored programs written in Python, C or Fortran.  For full
documentation of the NetCDF standard see https://www.unidata.ucar.edu/software/netcdf/.
Documentation is provided for the
:ref:`Python routines<raw-data-conversion>`
used in the processing chain leading to these Level 1 data.

A set of height-time quicklooks is also provided to enable the user to identify
periods of particular interest. An example is given below.

.. figure:: _static/example_quicklook_w-band.png
	   :width: 700 px
	   :align: center


:ref:`Python routines<quicklook-generation>` are provided for the production
of such quicklooks.  Specifically, they use the time-series data to generate
plots of radar equivalent reflectivity and of pulse-pair estimates of the
Doppler velocity.  No velocity unfolding is performed.



..
  .. note::

    Near real-time Cloudnet data can be accessed at https://cloudnet.fmi.fi.


..
  See also:
