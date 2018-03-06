#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Simple standalone tests of netcdf_aux module
# 1) Inject time frames into data set newly open data set
# 2) Append time frames into existing data set
# export PYTHONPATH="../..":${PYTHONPATH}
# -----------------------------------------------------------------------------

from numpy import *
from datetime import *
import netcdf_aux
import netCDF4  as netcdf
ncf = netcdf.Dataset("test4.nc", "w")

lon   = linspace(2, 3, 4)
lat   = linspace(55, 56, 5)
depth = linspace(0,10,3)
time  = None # leave time dimension unlimited 

info = {'lon': lon,
        'lat': lat,
        'time': time,
        'time_unit': "seconds",
        'time_offset': date(2017,3,5),
        'variable_name':"sPHYT",
        'units':'something',
        'long_name':'daily surface phytoplankton'}


data0 = random.random((1, len(lon), len(lat))) # have 1 as initial time dimension
netcdf_aux.write_lonlatdata_in_COARDS_format(ncf, data0, info, data_layout="txy")

# --- write a specific time frame to open set
data  = random.random((len(lon), len(lat)))
info['time'] = 28
netcdf_aux.write_lonlatdata_in_COARDS_format(ncf, data, info, data_layout="xy", time_frame_number=3)
ncf.close()


# --- inject time frame in existing set
ncf = netcdf.Dataset("test4.nc", "a")
info['time'] = 12
netcdf_aux.write_lonlatdata_in_COARDS_format(ncf, data, info, data_layout="xy", time_frame_number=2)
ncf.close()
