#!/usr/bin/env python

import sys; sys.path[1:1] = ["../../.."]; print "fix sys.path hack"
from GridWetData.grids import *
import re

################################################################################
#               Raw CMEMS NWS data query utilities
################################################################################
#  
#  Since grid information is contained in any data file, NWSGrids does not have
#  a separate grid descriptor, but grid information is extracetd from any data file
#  The dictionary below is a mapping of the known data variable names in the CMEMS NWS
#  data set. It is a user convenience used for inferring the topography from a
#  data set so that the user does not have to provide an explicit variable name
#  that should be used for probing the topography, but can just point to a data file


_verbose = False

NWS_variablenames_3D = {"temperature"             : "votemper",
                        "salinity"                : "vosaline",
                        "east_current"            : "vozocrtx",
                        "north_Current"           : "vomecrty",
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


    
_time_offset = datetime(2014,1,1)

# for parsing the time:units attribute

daysec_time_offset_designation = re.compile("(?P<tunit>\w+)[ ]+since[ ]+(?P<toffset>\d{4,4}-\d{2,2}-\d{2,2} \d{2,2}:\d{2,2}:\d{2,2})")
day_time_offset_designation    = re.compile("(?P<tunit>\w+)[ ]+since[ ]+(?P<toffset>\d{4,4}-\d{2,2}-\d{2,2})")                                
seconds_per_time_interval      = {"seconds": 1,
                                  "hours"  : 3600,
                                  "days"   : 86400}


def _to_datetime(thash):
    #  ------------------------------------------------------------
    ## convert timehash scalar to datetime object
    #  ------------------------------------------------------------
    return _time_offset + timedelta(seconds=int(thash))



def _cc2lw(cc):
    #  ------------------------------------------------------------
    ## Convert array cc of cell centers to array lw of layer widths
    #
    #  @param cc            array of cell center depths
    #  @return              corresponding array of layer withs
    #
    #  Surface is at z=0. Notice that it is not in general possible
    #  to convert an array of cell centers to an implied list of
    #  layer withs. The algortihm below always gives a result, but
    #  in general, some layers may be negative for arbitrary input
    #  ------------------------------------------------------------
    n  = len(cc)
    lw = ones(n, float)
    lw[0] = 2*cc[0]
    for l in range(1,n):
        lw[l] = 2*(cc[l]-cc[l-1] - 0.5*lw[l-1])
    assert all(lw >= 0) # test solution is meaningfull. Allow skin layers with lw = 0
    return lw


# ==============================================================================
#    Dress up basic netcdf data set class with additional handy methods
#
#    Use compositional relation, not heritance, because all user-set sub class
#    attributes within the netcdf.Dataset instance scope are written to file as
#    netCDF attributes, in confict with data being opened in reading mode
#    Use auto relay to netcdf.Dataset, if method/attribute not present, as pseudo heritance
# 
#    base class: http://unidata.github.io/netcdf4-python/
# ==============================================================================
class NWSDataSet:
    def __init__(self, *args, **kwargs):
        self.ncf = netcdf.Dataset(*args, **kwargs)
        self._update_time_hash()

    def close(self):
        self.ncf.close()

    def has_variable(self, vname): return self.ncf.variables.has_key(vname)

        
    def _update_time_hash(self):
        # -----------------------------------------------------
        # Parse time:units attribute and create a time hash : _timehash 
        # based on the variable time for data frames in the set
        # _timehash has units seconds since _time_offset (module constant)
        # 
        # Time designation is slightly different bewteen data products; examples :
        #   MetO-NWS-BIO-dm-DOXY.h:		   time:units = "seconds since 2017-06-27 00:00:00" ;
        #   MetO-NWS-PHYS-hi-TEM.h:	  	   time:units = "seconds since 2015-07-03 00:00:00" ;
        #   MetO-NWS-REAN-PHYS-daily-TEM.h: time:units = "days since 1985-01-01" ;
        #
        # It may be considered, if time frames with day annotation should be assigned to
        # noon that day, rather than mid night
        # Assumes time frames are stored in time order and that a
        # time point is only represented once
        # -----------------------------------------------------
        unitattr = self.ncf.variables["time"].units

        onlyday = day_time_offset_designation.match(unitattr)
        secday  = daysec_time_offset_designation.match(unitattr)
        if  secday: # test first, because onlyday will match when secday does
            toff = datetime.strptime(secday.group('toffset'), "%Y-%m-%d %H:%M:%S")
            tsca = seconds_per_time_interval[ secday.group('tunit') ]
        elif onlyday:
            toff = datetime.strptime(onlyday.group('toffset'), "%Y-%m-%d") # implicitly 00:00:00
            tsca = seconds_per_time_interval[ onlyday.group('tunit') ]
        else: # do not capture the situation with matches for both sunit and dunit ...
            msg  = "In data set %s: attribute time:units has unknown format: %s" % (self.ncf.file, unitattr)
            raise exceptions.ValueError(msg)
        #
        dtsec0 = int(0.5 + (toff - _time_offset).total_seconds())
        t = self.ncf.variables["time"][:]
        if t.dtype == float:
                t = asarray(0.5 + t, int) # 0.5 : protect against down rounding
        self.timehash = tsca*t + dtsec0   # int vector

            
    def get_time_interpolation_bracket(self, when):
        # -----------------------------------------------------
        # Pick time frames for time interpolation and compute interpolation weights
        # Return None, if when is not interior (or at an end point)
        # Return also a bracket in the case where when is coincident with a time point
        # to avoid this special case at interplation
        #
        # if i = searchsorted(a, v) then the bracket for v into a is (i-1, i)
        # -----------------------------------------------------
        th = int((when - _time_offset).total_seconds())
        if th < self.timehash[0] or th > self.timehash[-1]:
            return None
        else:
            i1 = searchsorted(self.timehash, th) # right delimiter, i1 >= 0
            i1 = min(len(self.timehash)-1, i1)   # 0 <= i1 <= len(self.timehash)-1
            i0 = max(0, i1-1)                    # 0 <= i0 <= len(self.timehash)-1
            dt = self.timehash[i1]-self.timehash[i0]
            if dt > 0:
                w1 = 1.0*(th - self.timehash[i0])/dt # force float division
                w0 = 1.0 - w1
            else:
                w0 = w1 = 0.5  # splitting arbitrary
            return (i0,w0), (i1,w1)

        
    def _extract_2D_grid_params(self):
        #  -----------------------------------------------------
        ## Extract basic 2D descriptors from open netCDF file ncf
        #
        #  @param self       open netCDF file object
        #  @return           nx, ny, lon0, lat0, dlon, dlat (basic grid descriptors)
        #  -----------------------------------------------------
        nx         = self.ncf.dimensions["lon"].size # integer 
        ny         = self.ncf.dimensions["lat"].size # integer 
        lonvec     = self.ncf.variables["lon"][:]
        latvec     = self.ncf.variables["lat"][:]
        dlon       = lonvec[1]-lonvec[0]
        dlat       = latvec[1]-latvec[0]
        assert dlon > 0
        assert dlat > 0
        return nx, ny, lonvec[0], latvec[0], dlon, dlat

        
    def _extract_3D_cellw0(self, diagvarname=NWS_variablenames_3D):
        #  ------------------------------------------------------------------------
        ## Extract cell widths cellw0(nx,ny,nz) from open netCDF file ncf
        #  @param self            open netCDF file
        #  @param diagvarname    name of variable from which to detect topography (optional)
        #  @return               cellw0(nx,ny,nz)  (array of cell widths, -1 for dry cells)
        #   
        #  Assess topography by probing for mask==True at data variable diagvarname in ncf.
        #  If optional parameter diagvarname is not present, it defaults to an entry
        #  in a pre-coded catalogue of variable names to look for (NWS_variablenames).
        #  If lookup fails an exception is raised 
        #
        #  Data variables are apparently loaded as masked arrays by
        #  the netcdf4 API. Therefore detect topography by testing mask==True (dry)
        #  ------------------------------------------------------------------------
        nx, ny, lon0, lat0, dlon, dlat = self._extract_2D_grid_params()
        nz         = self.ncf.dimensions["depth"].size # integer
        #
        # ------ identify topography probing variable ------
        #
        if (diagvarname == NWS_variablenames_3D):   # search pre-coded catalogue
            for vnam in NWS_variablenames_3D.values():
                if vnam in self.ncf.variables.keys(): # found a match in pre-coded catalogue
                    break
            else:
                raise exceptions.LookupError("_extract_3D_cellw0: netCDF file %s does not contain a known data variable" % self.ncf.filepath())
        else:
            vnam = diagvarname
        var           = self.ncf.variables[vnam][0,:]     # topography test variable, asked array, shape == (time, depth, lat, lon)
        var           = swapaxes(var, 0, 2)          # now shape = (lon, lat, depth)
        #
        ccz           = self.ncf.variables["depth"][:]    # layer centers
        cellw0        = zeros((nx,ny,nz), float)
        cellw0[:,:,:] = _cc2lw(ccz)                  # convert cell center positions to layer widths
        cellw0        = where(var.mask, -1, cellw0)  # dry cells assigned -1
        return cellw0

    def get_time_frame(self, vname, itime=0):
        ## Extract frame itime of variable vname and swap axes, corresponding to GWD data layout
        var   = self.ncf.variables[vname][itime,:]
        var   = swapaxes(var, 0, 2) # cast to GWD layout
        if _verbose:
            print "get_time_frame: exporting frame %d" % itime
        return var.data             # strip mask
        
    def get_time(self, itime=0):
        ## return datetime corresponding to frame itime 
        return _to_datetime(self.timehash[itime])
    
    def get_number_of_frames(self, vname):
        return len(self.ncf.variables[vname])
        
        
#################### self test #########################
#
#  mid North Sea test point: ix,iy = 210,240
#
if __name__ == "__main__":
    _verbose = True
    raw = NWSDataSet("../../../tmp/MetO-NWS-PHYS-hi-TEM.nc")
    v = raw.get_time_frame("votemper", 3)
    print v[210,240,:]
    #print raw.get_time_interpolation_bracket(datetime(2015,7,3,23,59))
    #NWSGrid_3D("../../../tmp/MetO-NWS-PHYS-hi-TEM.nc", "votemper")
    #g3D = NWSGrid_3D("../../../tmp/MetO-NWS-PHYS-hi-TEM.nc")
    
    
                    
         
