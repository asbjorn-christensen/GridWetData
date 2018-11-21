#!/usr/bin/env python
#######################################################################################
#
# Extract bottom currents (lowest wet grid cell) from hydrographic subset by looping
# with GridWetData.HBM.DataManager.loop_over(...)
#
# Creates netcdf data variables (ubott, vbott) in output file out.nc
# 
# If GridWetData is not in default search path, then (in bash) apply:
#    export PYTHONPATH=<PATH_TO_GridWetData>:${PYTHONPATH}
#######################################################################################

from GridWetData import *                 # import the GridWetData environment
from GridWetData.HBM import DataManager, HBMGrid_3D, GridData_3D_EastStaggerd, GridData_3D_SouthStaggerd

dmg = DataManager("/home/data/GUDP-VIND_test/preprocessed_data", HBMGrid_3D)   # could be sys.argv[1]

import  netCDF4 as netcdf

ncuv = netcdf.Dataset("out.nc", "w")                                           # could be sys.argv[2]
times = []
i = 0
for (tim, data3d) in dmg.loop_over(GridData_3D_EastStaggerd, "u", "ns"):
        times.append(tim)
        dt = tim-times[0]
        bLayer = data3d.get_bottom_layer()
        bLayer.write_data_as_netCDF(ncuv, index=i, dataParam = "ubott")
        i += 1
        print "u:", i, tim
times = []
i = 0
for (tim, data3d) in dmg.loop_over(GridData_3D_SouthStaggerd, "v", "ns"):
        times.append(tim)
        dt = tim-times[0]
        bLayer = data3d.get_bottom_layer()
        bLayer.write_data_as_netCDF(ncuv, index=i, dataParam = "vbott")
        i += 1
        print "v:", i, tim
