#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Simple standalone tests of netcdf_aux module
# 1) Dump a single xy and xyz frame to file, no time reference
#
# export PYTHONPATH="../..":${PYTHONPATH}
# -----------------------------------------------------------------------------
from numpy import *
from datetime import *
import netcdf_aux
import netCDF4  as netcdf
ncf = netcdf.Dataset("test.nc", "w")

lon   = linspace(2, 3, 4)
lat   = linspace(55, 56, 5)
depth = linspace(0,10,3)
time  = range(1000, 3900, 900)

info = {'lon': lon,
        'lat': lat,
        'variable_name':"sPHYT",
        'units':'something',
        'long_name':'daily surface phytoplankton'}

# -------- just XY data --------
data  = random.random( (len(lon), len(lat)) )
ncf = netcdf.Dataset("test1.nc", "w")
netcdf_aux.write_lonlatdata_in_COARDS_format(ncf, data, info, data_layout="xy")
ncf.close()

# -------- XYZ data --------
info['depth'] = depth
data = random.random( (len(lon), len(lat), len(depth)) )
ncf = netcdf.Dataset("test2.nc", "w")
netcdf_aux.write_lonlatdata_in_COARDS_format(ncf, data, info, data_layout="xyz")
ncf.close()

# -------- TXYZ data --------
info['time']        = time
info['time_unit']   = "seconds"        # optional
info['time_offset'] = date(2017,3,5)   # optional

data = random.random( (len(time), len(lon), len(lat), len(depth)) )
ncf = netcdf.Dataset("test3.nc", "w")
netcdf_aux.write_lonlatdata_in_COARDS_format(ncf, data, info, data_layout="txyz")
ncf.close()
