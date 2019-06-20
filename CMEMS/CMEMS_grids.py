#!/usr/bin/env python

from GridWetData.grids import *
import CMEMS_dataformat

# ==============================================================================
#  Shared grid classes for these CMEMS model output:
#   




# ==============================================================================
## 2D surface/at-depth  CMEMS grid
# sub class must specify valid names for topography resolution setting attribute resolve_topo_varnames_2D
class CMEMS_Grid_2D(LonLatGrid):
    #  ------------------------------------------------------------------------------
    resolve_topo_varnames_2D = None   # set class attribute (sub class must specify valid names for topography resolution)
    #  ------------------------------------------------------------------------------
    ## constructor from file data
    #  @param self    The object pointer
    #  @param fname   file name of an CMEMS data netcdf file (contains grid info)
    #
    def __init__(self, fname):
        ncfile     = CMEMS_dataformat.CMEMS_DataSet(fname)  # includes invoking _extract_2D_grid_params
        ## file name corresponding to loaded data set
        self.fname = fname
        #
        LonLatGrid.__init__(self, ncfile.nx, ncfile.ny, ncfile.lon0, ncfile.lat0, ncfile.dlon, ncfile.dlat)
        ncfile.close()
        #
        

#  =====================================================================
## 3D  CMEMS grid 
#  compositional relation to the 2D surface CMEMS grid
#  CMEMS grid data are interpolated to specific depths at generation time,
#  
#  sub class must specify valid names for topography resolution setting attribute resolve_topo_varnames_3D
class CMEMS_Grid_3D(LonLatZGrid):
    #  ------------------------------------------------------------------------------
    resolve_topo_varnames_3D = None     # set class attribute (sub class must specify valid names for topography resolution)
    #  ------------------------------------------------------------------------------
    ## constructor from file data
    #  @param self    The object pointer
    #  @param fname   file name of an CMEMS data netcdf file (contains grid info)
    #  @param diagvar Name of diagnostic variable, to deduce topography from (optional), otherwise apply resolve_topo_varnames_3D
    #
    def __init__(self, fname, diagvar=None):
        ncfile     = CMEMS_dataformat.CMEMS_DataSet(fname)   # invokes _extract_2D_grid_params
        ## file name corresponding to loaded data set
        self.fname = fname
        if diagvar is None:
            cellw0  = ncfile._extract_3D_cellw0(self.resolve_topo_varnames_3D) 
        else:
            cellw0  = ncfile._extract_3D_cellw0(diagvar)
        #
        LonLatZGrid.__init__(self, ncfile.nx, ncfile.ny, ncfile.lon0, ncfile.lat0, ncfile.dlon, ncfile.dlat, cellw0)    
        ncfile.close()
    
        

        
#################### self test #########################
#
#  mid North Sea test point: ix,iy = 210,240
#
if __name__ == "__main__":
    _verbose = True
    g2D = CMEMS_Grid_2D("data/MetO-NWS-PHYS-hi-TEM.nc")
    print "nx=%d ny=%d" % (g2D.nx, g2D.ny)
    #
    g3D = CMEMS_Grid_3D("data/MetO-NWS-PHYS-hi-TEM.nc", "votemper")
    for x in arange(-1, 12, 0.001):
        try:
            wd = g3D.interpolate_wdepth((x,54.13))
            print x,wd
        except:
            pass
    
                    
         
