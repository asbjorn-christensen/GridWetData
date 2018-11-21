#!/usr/bin/env python
#######################################################################################
#
# Generate a time series of a set of time frames by looping over data files 
# from GridWetData.HBM.DataManager
# 
# If GridWetData is not in default search path, then (in bash) apply:
#    export PYTHONPATH=<PATH_TO_GridWetData>:${PYTHONPATH}
#######################################################################################

from GridWetData import *                 # import the GridWetData environment
from GridWetData.HBM import DataManager, HBMGrid_3D

dmg = DataManager("/home/data/GUDP-VIND_test/preprocessed_data", HBMGrid_3D)  # could be sys.argv[1]

import  netCDF4 as netcdf
ncTempElevation = netcdf.Dataset("out.nc", "w")                               # could be sys.argv[2]   
times  = []
i = 0
for (thash, fname, tim) in zip(*dmg.datasets["z"]["ns"]):
        times.append(tim)
        #print tim, fname
        dt = tim-times[0]
        zindex = GridData_2D(dmg.grids["ns"].surf_grid, fname)
        #print zindex.data.min(), zindex.data.max()
        zindex.write_data_as_netCDF(ncTempElevation, index=i, dataParam="elevation")
        i += 1
