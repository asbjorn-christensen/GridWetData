#!/usr/bin/env python
###############################################################################
#   Essential tutorial demonstrating the GridWetData API
#
#   Each demonstration section is independent
#   The examples assume data and grid definitions in a directory "DMI_data",
#   where the kode is executed
###############################################################################
import os
import sys; sys.path[1:1] = [os.pardir]   # allow to see this directory from above, if this script is run in source directory

from GridWetData import *                 # import the GridWetData environment

###############################################################################
#
#   1)  Single point interpolations of data somewhere at sometime
#
#       As example, we use ocean temperature
#       First argument of data_manager.DataManager is the directory with
#       data and grid descriptors; second argument (here DMIGrid_3D) is the grid class that
#       should be used to represent grids
#       TemperatureData_3DwithTime is a grid data sub class that can resolve
#       which specific data files should be loaded (when needed)
#       position is the position, where we assess the temperature, as
#       any indexable type as (lon[deg E], lat[deg N], meters below the sea surface)
#       Temperature is assessed using the call hook of the TemperatureData_3DwithTime
#       instance: tt(position, when), where when is a datetime instance
#       Necessary data for interpoaltion is first loaded, when the time argument is received, 
#       to avoid loading excessive amount of data in the data base
#       The second interpolation call illustrates the caching property of
#       GridData_3DwithTime sub classes (like TemperatureData_3DwithTime).
#       GridData_3DwithTime determines that cache data is sufficient to
#       serve the requested second interpolation, and does not load more data from the data base
#
###############################################################################


dmg = data_manager.DataManager("DMI_data", DMIGrid_3D)  # interface to data base
tt  = TemperatureData_3DwithTime(dmg)                   # interface to temperature

pos = [4.2083, 56.025, 1.2]  #  lon[deg E], lat[deg N], meters below the sea surface
print tt(pos, datetime(2015, 05, 20, 00, 16, 00)), " - should give 8.55091339007" # interpolation using call hook
print tt(pos, datetime(2015, 05, 20, 00, 36, 00)), " - should give 8.54531238808" # interpolation using call hook


###############################################################################
#
#   2)  Vertical scan down the water column of data somewhere at sometime
#       by manually selected data frames (i.e. not using the data manager)
#
#       This time, we'll not use the data manager, but do it manually to see
#       what goes on behind the scene. As example, we'll try to assess the
#       salinity near the sea bed
#       First line instantiates a DMIGrid_3D instance from a grid descriptor file
#       Second line updates grid topography corresponding to a specific sea level elevation data frame
#       Third line loads a salinity data set from a specific file (you are responsible
#       for the correspondance between data and dynamic sea level when circumventing
#       the data manager and *_3DwithTime classes)
#       Fifth line assesses full water depth at selected position (corresponding to
#       loaded sea level elevation)
#       Line 6-9 scans down the water column, interpolating salinity along the way
###############################################################################

g3D = DMIGrid_3D(       os.path.join("DMI_data", "ns_grid.nc"))  # load this grid
g3D.update_sea_level(   os.path.join("DMI_data", "z_ns_2015_05_20_00_15_00.nc")) # manually select frame
salt = GridData_3D(g3D, os.path.join("DMI_data", "s_ns_2015_05_20_00_15_00.nc")) # manually select frame

pos = [4.2083, 56.025, 0] #  lon[deg E], lat[deg N], meters below the sea surface
cwd = g3D.interpolate_wdepth(pos) # assess water full depth at position pos
for depth in arange(0, cwd, 0.01):
    pos[2] = depth
    print depth, salt(pos) # only a position argument for GridData_3D interpolation call hook


###############################################################################
#
#   3)  Surface layer / bottom layer / column average extraction of data on native grid
#
#       This time, we'll use the DataManager to create a GridData_3D
#       by interpolating two time frames. snapshot_at() interpolates
#       both data and sea level interpolation between time frames 
###############################################################################

dmg    = data_manager.DataManager("DMI_data", DMIGrid_3D)   # interface to data base
tt     = TemperatureData_3DwithTime(dmg)                    # interface to temperature
tframe = tt.snapshot_at(datetime(2015, 05, 20, 00, 16, 00)) # GridData_3D instance 

tsurf = tframe.get_surface_layer()
tbott = tframe.get_bottom_layer()
tavg  = tframe.get_vertical_average()

print "temperature in water column             :", tframe.data[100,200,:]
print "cell thickness in water column          :", tframe.grid.cellw[100,200]
print "surf cell temperature from layer        :", tsurf.data[100,200]
print "bottom cell temperature from layer      :", tbott.data[100,200]
print "vertical average temperature from layer :", tavg.data[100,200]

###############################################################################
#
#   4)  Interpolation of data at custom provided mesh at sometime   
#
#   GridData_3D and GridData_3DwithTime call hooks accepts a list of positions
#   and returns corresponding list of interplated values
###############################################################################

dmg    = data_manager.DataManager("DMI_data", DMIGrid_3D)   # interface to data base
tt     = TemperatureData_3DwithTime(dmg)

my_mesh = []  # create example mesh
for x in arange(4, 6, 0.01):
    my_mesh.append([x, 56.025, 1.33])
    
data_at_my_mesh = tt(my_mesh, datetime(2015, 05, 20, 00, 16, 00))
print data_at_my_mesh

###############################################################################
#
#   5)  Derived data: thermal front index layer on native grid
#
#   The module derived_properties contain functionality to
#   generate derived data from raw physical data sets
###############################################################################



dmg     = data_manager.DataManager("DMI_data", DMIGrid_3D)   # interface to data base
temp    = TemperatureData_3DwithTime(dmg).snapshot_at(datetime(2015, 05, 20, 00, 16, 00))     

sindex  = derived_layers.StratificationFront(temp) # GridData_2D instance
sindex.write_data("sindex.nc")                     # save data to file - netcdf format in this case

