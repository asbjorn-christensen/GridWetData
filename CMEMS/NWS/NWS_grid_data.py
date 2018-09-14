#!/usr/bin/env python

from GridWetData.CMEMS.CMEMS_grid_data import *
from NWS_grids                         import *

_verbose          = False  # control debugging output

############################################################################
#   Shallow subclassing of CMEMS_grid_data classes, applying
#   NWS_Grid_2D/NWS_Grid_3D as default, if a grid is not provided explicitly
############################################################################
#    
class NWS_GridData_3DwithTime(CMEMS_GridData_3DwithTime):
    #  -----------------------------------------------------
    ## constructor
    #  Invoke super class constructor with grid = NWS_Grid_3D as default
    #  @param self      The object pointer
    #  @param fname     netcdf data set file name  
    #  @param varname   variable name to look for (optional)
    #                   If absent (and class attributre prop is not set) the data set
    #                   is scanned for known variable names
    #  @param gridclass Grid instance (optional). Default: instantiate a NWS_grid
    def __init__(self, fname, varname=None, grid=None):
        if grid is None:
            CMEMS_GridData_3DwithTime.__init__(self, fname, varname, NWS_Grid_3D(fname))
        else:    
            CMEMS_GridData_3DwithTime.__init__(self, fname, varname, grid)


#  ====================================================================================
## Super class for interpolation of vector data in space+time by first interpolating space, then time
#  specific for 3D + interaction with DataManager 
# <b> Roles </b> 
#    \li keep track of data cache, decide whether to reload
#    \li interact with data manager
#    \li call hook is redirectd to interpolate, so instances are callable as obj(where, when)
#
# compositional relation to GridVector_3D
# Attributes set by sub classes:
#       prop: data property types vector (set by sub classes)
#    
# Inherit update_cache from super class
#
class NWS_GridVector_3DwithTime(CMEMS_GridVector_3DwithTime):
    #  -----------------------------------------------------
    ## constructor
    #  Invoke super class constructor with grid = NWS_Grid_3D as default
    #  @param self       The object pointer
    #  @param fname      netcdf data set file name  
    #  @param varnames   variable names to look for (vector of strings) (optional)
    #                    If absent, class attributre prop must be set
    #                    No scanning for known variable names
    #                    for each token "zero_data", a zero-valued field is inserted in that position
    #                    (allow to e.g. uprank a horizontal field to 3D)
    #  @param gridclass  Grid instance (optional). Default: instantiate a CMEMS_grid
    def __init__(self, fname, varname=None, grid=None):
        if grid is None:
            CMEMS_GridVector_3DwithTime.__init__(self, fname, varname, NWS_Grid_3D(fname))
        else:    
            CMEMS_GridVector_3DwithTime.__init__(self, fname, varname, grid)
        
# class NWS_GridVector_2DwithTime(CMEMS_GridVector_withTime):  deferred until needed

#  ====================================================================================
## Super class for interpolation of 2D data in space+time by first interpolating space, then time
#  

class NWS_GridData_2DwithTime(CMEMS_GridData_2DwithTime):
    #  -----------------------------------------------------
    ## constructor
    #  Invoke super class constructor with grid = NWS_Grid_2D as default
    #  @param self      The object pointer
    #  @param fname     netcdf data set file name  
    #  @param varname   variable name to look for (optional)
    #                   If absent (and class attributre prop is not set) the data set
    #                   is scanned for known variable names
    def __init__(self, fname, varname=None, grid=None, get_wetmask=False):
        if grid is None:
            CMEMS_GridData_2DwithTime.__init__(self, fname, varname=varname,
                                                     grid=NWS_Grid_2D(fname),
                                                     get_wetmask=get_wetmask)
        else:    
            CMEMS_GridData_2DwithTime.__init__(self, fname, varname=varname,
                                                    grid=grid,
                                                     get_wetmask=get_wetmask)
        
        
    

##########################################################################################
#    Specific elementary property classes
##########################################################################################

# The stagger method will trigger that load_frame restagger data to
# cell centered data after load. If class does not have this method, no
# restaggering is performed


## GridData_2D sub classes 
class NWS_SeaFloorTemperature_2DwithTime(NWS_GridData_2DwithTime):
    ## data property type 
    prop = "sotemper"   # set class attribute, passed to when soliciting data

class NWS_MixedLayerThickness_2DwithTime(NWS_GridData_2DwithTime):
    ## data property type 
    prop = "karamld"   # set class attribute, passed to when soliciting data

class NWS_SeaSurfaceHeight_2DwithTime(NWS_GridData_2DwithTime):
    ## data property type 
    prop = "sossheig"   # set class attribute, passed to when soliciting data


    
    
## GridData_3D sub classes for temperature data in both space and time
class NWS_Temperature_3DwithTime(NWS_GridData_3DwithTime):
    ## data property type 
    prop = "votemper"   # set class attribute, passed to when soliciting data

class NWS_Salinity_3DwithTime(NWS_GridData_3DwithTime):
    ## data property type 
    prop = "vosaline"   # set class attribute, passed to when soliciting data
                      
class NWS_VolumeBeamAttenuation_3DwithTime(NWS_GridData_3DwithTime):
    ## data property type 
    prop = "attn"   # set class attribute, passed to when soliciting data                        
                
class NWS_Chlorophyll_3DwithTime(NWS_GridData_3DwithTime):
    ## data property type 
    prop = "CHL"   # set class attribute, passed to when soliciting data                        
                      
class NWS_DissolvedOxygen_3DwithTime(NWS_GridData_3DwithTime):
    ## data property type 
    prop = "O2o"   # set class attribute, passed to when soliciting data                        
                   
class NWS_Nitrate_3DwithTime(NWS_GridData_3DwithTime):
    ## data property type 
    prop =  "N3n" # set class attribute, passed to when soliciting data 

class NWS_Phosphate_3DwithTime(NWS_GridData_3DwithTime):
    ## data property type 
    prop =  "N1p"  # set class attribute, passed to when soliciting data                        
                    
class NWS_Phytoplankton_3DwithTime(NWS_GridData_3DwithTime):
    ## data property type 
    prop = "PhytoC"   # set class attribute, passed to when soliciting data                        

class NWS_PrimaryProductivity_3DwithTime(NWS_GridData_3DwithTime):
    ## data property type 
    prop = "netPP"   # set class attribute, passed to when soliciting data                        



## GridVector_3D sub classes for temperature data in both space and time    
class NWS_HorizontalCurrents_3DwithTime(NWS_GridVector_3DwithTime):
    ## data property type 
    prop = ("vozocrtx", "vomecrty")   # set class attribute, passed to when soliciting data               

class NWS_Currents_vertical0_3DwithTime(NWS_GridVector_3DwithTime):
    ## vertical component currently set to zero
    prop = ("vozocrtx", "vomecrty", "zero_data")   # set class attribute, passed to when soliciting data
    

                        
#################### self test #########################
#
#  mid North Sea test point: ix,iy = 100,150  <->  lon,lat = 4.2083, 56.025
#
if __name__ == "__main__":
    temp  = NWS_Temperature_3DwithTime("data/MetO-NWS-PHYS-hi-TEM.nc")
    temp0 = temp.load_frame(0)
    pos   = [4.2083, 56.025, 0]
    ## --- test vertical get_vertex_indices
    maxwd = temp.grid.interpolate_wdepth(pos)
    for wd in arange(0, maxwd, 0.01):
        pos[2] = wd
        #print wd, temp0(pos)
    #sys.exit()
    ## import data_manager
    ## dmg = data_manager.DataManager("../DMI_data", HBMGrid_3D)
    ## tt = TemperatureData_3DwithTime(dmg)
    ## pos = [4.2083, 56.025, 1.2]
    ## tt.interpolate(pos, datetime(2015, 05, 20, 00, 16, 00))
    ## print tt
    ## tt.interpolate(pos, datetime(2015, 05, 20, 00, 26, 00))
    ## print tt
    ## tt interpolate(pos, datetime(2015, 05, 20, 00, 36, 00))
    ## print tt
    pos  = [4.2083, 56.025, 0]
    now  = now  = datetime(2017, 6, 29, 5, 0)
    for i in range(120):
        print i, temp(pos, now)
        now = now + timedelta(seconds=60)
