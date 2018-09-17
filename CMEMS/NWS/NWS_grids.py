#!/usr/bin/env python

from GridWetData.CMEMS.CMEMS_grids import *


NWS_variablenames_3D = {"temperature"             : "votemper",
                        "salinity"                : "vosaline",
                        "east_current"            : "vozocrtx",
                        "north_current"           : "vomecrty",
                        "volume_beam_attenuation" : "attn",
                        "chlorophyll"             : "CHL",
                        "dissolved_oxygen"        : "O2o",
                        "nitrate"                 : "N3n",
                        "phosphate"               : "N1p",
                        "phytoplankton"           : "PhytoC",
                        "primary_productivity"    : "netPP"} 

NWS_variablenames_2D = {"sea_floor_temperature" : "sotemper", 
                        "mixed_layer_thickness" : "karamld", 
                        "sea_surface_height"    : "sossheig"}


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
## 2D surface/at-depth  NWS grid
class NWS_Grid_2D(CMEMS_Grid_2D):
    #  ----------------------------------------------------------------------------------------
    resolve_topo_varnames_2D = NWS_variablenames_2D  # valid variables to resolve 2D topography
    #  ----------------------------------------------------------------------------------------
    

#  =====================================================================
## 3D  NWS grid

class NWS_Grid_3D(CMEMS_Grid_3D):
    #  ----------------------------------------------------------------------------------------
    resolve_topo_varnames_3D = NWS_variablenames_3D  # valid variables to resolve 3D topography
    #  ----------------------------------------------------------------------------------------
    
    
        
if __name__ == "__main__":
    g2D = NWS_Grid_2D("data/MetO-NWS-PHYS-hi-TEM.nc")
    print "nx=%d ny=%d" % (g2D.nx, g2D.ny)
    g3D = NWS_Grid_3D("data/MetO-NWS-PHYS-hi-TEM.nc")
    print "nx=%d ny=%d nz=%d" % (g3D.nx, g3D.ny, g3D.nz)
