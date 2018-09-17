#!/usr/bin/env python

from GridWetData.CMEMS.CMEMS_grids import *


#   float mlotst(time, lat, lon)        standard_name = "ocean_mixed_layer_thickness_defined_by_sigma_theta" ;
#	float siconc(time, lat, lon)        standard_name = "sea_ice_area_fraction" ;
#	float sob(time, lat, lon)           standard_name = "sea_water_salinity" ;
#	float sla(time, lat, lon)           standard_name = "sea_surface_height_above_sea_level" ;
#   float sithick(time, lat, lon)       standard_name = "sea_ice_thickness" ;
#	float bottomT(time, lat, lon)       standard_name = "sea_water_potential_temperature_at_sea_floor" ;

#	float vo(time, depth, lat, lon)     standard_name = "northward_sea_water_velocity" ;
#	float thetao(time, depth, lat, lon) standard_name = "sea_water_potential_temperature" ;
#	float uo(time, depth, lat, lon)     standard_name = "eastward_sea_water_velocity" ;
#	float so(time, depth, lat, lon)     standard_name = "sea_water_salinity" ;
    
BAL_variablenames_3D = {"temperature"             : "thetao",
                        "salinity"                : "so",
                        "east_current"            : "uo",
                        "north_current"           : "vo"} 

BAL_variablenames_2D = {"mixed_layer_thickness" : "mlotst", 
                        "sea_ice_area_fraction" : "siconc", 
                        "sea_water_salinity"    : "sob",
                        "sea_surface_height"    : "sla",
                        "sea_ice_thickness"     : "sithick",
                        "sea_water_potential_temperature_at_sea_floor"    : "bottomT"
                        }


# ==============================================================================


# ==============================================================================
## 2D surface/at-depth  BAL grid
class BAL_Grid_2D(CMEMS_Grid_2D):
    #  ----------------------------------------------------------------------------------------
    resolve_topo_varnames_2D = BAL_variablenames_2D  # valid variables to resolve 2D topography
    #  ----------------------------------------------------------------------------------------
    

#  =====================================================================
## 3D  BAL grid

class BAL_Grid_3D(CMEMS_Grid_3D):
    #  ----------------------------------------------------------------------------------------
    resolve_topo_varnames_3D = BAL_variablenames_3D  # valid variables to resolve 3D topography
    #  ----------------------------------------------------------------------------------------
    
    
        
if __name__ == "__main__":
    g2D = BAL_Grid_2D("data/dataset-bal-analysis-forecast-phy-hourly.nc")
    print "nx=%d ny=%d" % (g2D.nx, g2D.ny)
    g3D = BAL_Grid_3D("data/dataset-bal-analysis-forecast-phy-hourly.nc")
    print "nx=%d ny=%d nz=%d" % (g3D.nx, g3D.ny, g3D.nz)
