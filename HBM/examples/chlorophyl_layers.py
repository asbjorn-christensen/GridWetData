#!/usr/bin/env python
#######################################################################################
#
# Extract misc chlorophyl attributes (strat front, bottom, surface, vertical average)
# by looping with GridWetData.HBM.DataManager.loop_over(...)
# 
# If GridWetData is not in default search path, then (in bash) apply:
#    export PYTHONPATH=<PATH_TO_GridWetData>:${PYTHONPATH}
#######################################################################################

from GridWetData import *                 # import the GridWetData environment
from GridWetData.HBM import DataManager, HBMGrid_3D

dmg = DataManager("/home/data/GUDP-VIND_test/preprocessed_data", HBMGrid_3D)  # could be sys.argv[1]

import  netCDF4 as netcdf


ncChlorofyl = netcdf.Dataset("out.nc", "w")                                   # could be sys.argv[2]
times = []
i = 0
for (tim, data3d) in dmg.loop_over(GridData_3D, "chl", "ns"):
        times.append(tim)
        dt = tim-times[0]
        fLayer  = derived_layers.StratificationFront(data3d)
        bLayer = data3d.get_bottom_layer()
        sLayer = data3d.get_surface_layer()
        aLayer = data3d.get_vertical_average()
        fLayer.write_data_as_netCDF(ncChlorofyl, index=i, dataParam = "front")
        bLayer.write_data_as_netCDF(ncChlorofyl, index=i, dataParam = "buttom")
        sLayer.write_data_as_netCDF(ncChlorofyl, index=i, dataParam = "surface")
        aLayer.write_data_as_netCDF(ncChlorofyl, index=i, dataParam = "average")
        i += 1
