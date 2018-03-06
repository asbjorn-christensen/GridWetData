#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Simple standalone tests of netcdf_aux module
# 1) Inject time frames into data set with known time limits
# 2) Append time frames into existing data set
# export PYTHONPATH="../..":${PYTHONPATH}
# -----------------------------------------------------------------------------

from numpy import *
from datetime import *
import netcdf_aux
import netCDF4  as netcdf
ncf = netcdf.Dataset("test5.nc", "w")

lon   = linspace(2, 3, 4)
lat   = linspace(55, 56, 5)
time  = range(1000, 3900, 900) # time grid known


info = {'lon': lon,
        'lat': lat,
        'time': time,
        'time_unit': "seconds",
        'time_offset': date(2017,3,5),
        'variable_name':"sPHYT",
        'units':'something',
        'long_name':'daily surface phytoplankton'}


data0 = zeros((len(time), len(lon), len(lat)))  # initialization on specific known time grid
netcdf_aux.write_lonlatdata_in_COARDS_format(ncf, data0, info, data_layout="txy")

# --- write a specific time frames to open set
for it in range(len(time)):
    data  = random.random((len(lon), len(lat)))
    #info['time'] = 28   # you may overwrite time by setting info['time']
    netcdf_aux.write_lonlatdata_in_COARDS_format(ncf, data, info, data_layout="xy", time_frame_number=it)
ncf.close()

