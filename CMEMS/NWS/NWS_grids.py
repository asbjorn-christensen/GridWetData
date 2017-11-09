#!/usr/bin/env python

import sys; sys.path[1:1] = ["../../.."]; print "fix sys.path hack"
from GridWetData.grids import *

import NWS_dataformat

# ==============================================================================
#  Shared grid classes for these CMEMS model output:
#
#    NORTHWESTSHELF_ANALYSIS_FORECAST_PHYS_004_001_b 
#    NORTHWESTSHELF_REANALYSIS_PHYS_004_009
#    NORTHWESTSHELF_ANALYSIS_FORECAST_BIO_004_002_b
#    NORTHWESTSHELF_REANALYSIS_BIO_004_011
#
#  NWS data does not currently have dynamic sealevel elevation




# ==============================================================================
## 2D surface wet point NWS grid
class NWSGrid_2D(LonLatGrid):
    #  -----------------------------------------------------
    ## constructor from file data
    #  @param self    The object pointer
    #  @param fname   file name of an NWS data netcdf file (contains grid info)
    #
    def __init__(self, fname):
        ncfile     = NWS_dataformat.NWSDataSet(fname)
        ## file name corresponding to loaded data set
        self.fname = fname
        nx, ny, lon0, lat0, dlon, dlat = ncfile._extract_2D_grid_params()
        ncfile.close()
        #
        LonLatGrid.__init__(self,  nx, ny, lon0, lat0, dlon, dlat)
        #
    

#  =====================================================================
## 3D  NWS grid
#  compositional relation to the 2D surface NWS grid
#  NWS grid data are interpolated to specific depths at generation time,
#  so vertical 
class NWSGrid_3D(LonLatZGrid):
    #  -----------------------------------------------------
    ## constructor from file data
    #  @param self    The object pointer
    #  @param fname   file name of an NWS data netcdf file (contains grid info)
    #  @param diagvar Name of diagnostic variable, to deduce topography from
    #
    def __init__(self, fname):
        ncfile     = NWS_dataformat.NWSDataSet(fname)
        ## file name corresponding to loaded data set
        self.fname = fname
        nx, ny, lon0, lat0, dlon, dlat = ncfile._extract_2D_grid_params()
        cellw0                         = ncfile._extract_3D_cellw0()
        ncfile.close()
        #
        LonLatZGrid.__init__(self,  nx, ny, lon0, lat0, dlon, dlat, cellw0)
    
    
        

        
#################### self test #########################
#
#  mid North Sea test point: ix,iy = 210,240
#
if __name__ == "__main__":
    _verbose = True
    #g3D = NWSGrid_3D("../../../tmp/MetO-NWS-PHYS-hi-TEM.nc", "votemper")
    g3D = NWSGrid_3D("../../../tmp/MetO-NWS-PHYS-hi-TEM.nc")
    for x in arange(-1, 12, 0.001):
        try:
            wd = g3D.interpolate_wdepth((x,54.13))
            print x,wd
        except:
            pass
    
                    
         
