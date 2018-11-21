#!/usr/bin/env python
#######################################################################################
#
# Reexport wet point data format into regular grid data format 
# Use DataManager inventory list to pick file name
# 
# If GridWetData is not in default search path, then (in bash) apply:
#    export PYTHONPATH=<PATH_TO_GridWetData>:${PYTHONPATH}
#######################################################################################

from GridWetData import *                 # import the GridWetData environment
from GridWetData.HBM import DataManager, HBMGrid_3D

dmg = DataManager("/home/data/GUDP-VIND_test/preprocessed_data", HBMGrid_3D)  # could be sys.argv[1]

fname = dmg.datasets["z"]["ns"][1][0]                                         # pick wet point data file
zindex = GridData_2D(dmg.grids["ns"].surf_grid, fname)
zindex.write_data_as_netCDF("out.nc")                                         # could be sys.argv[2]
