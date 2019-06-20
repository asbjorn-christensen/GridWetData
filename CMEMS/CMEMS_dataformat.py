#!/usr/bin/env python


from GridWetData.grids import *
import re

################################################################################
#               Raw CMEMS data query utilities
################################################################################
#  
#  Since grid information is contained in any data file, CMEMS_Grids does not load
#  a separate grid descriptor, but extracts grid information from any data file
#  Sub classes of CMEMS_Grids provide lookup variable name tables, corresponding to
#  specific data sets 
#
#  horizontal axes must be called either ('lon' and 'lat') OR ('longitude' and 'latitude')
#  not a hybrid combination; dimension and axis name must be identical
#  Data set contains only one set of horizontal axes, which are considered singletons
#
#  numpy.flip requires numpy version >= 1.12
################################################################################
_lon_names = ('lon', 'longitude')
_lat_names = ('lat', 'latitude')

_verbose = False

    
_time_offset = datetime(2014,1,1)  # module internal - not related to actual offset of time axis

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
class CMEMS_DataSet:
    def __init__(self, *args, **kwargs):
        self.ncf = netcdf.Dataset(*args, **kwargs)
        self._update_time_hash()
        self._extract_2D_grid_params()  # make sure flip_lon, flip_lat attrs are set
        
    def close(self):
        self.ncf.close()

    def has_variable(self, vname): return self.ncf.variables.has_key(vname)

        
    def _update_time_hash(self):
        # -----------------------------------------------------
        # Parse time:units attribute and create a time hash : timehash 
        # based on the variable time for data frames in the set
        # timehash has units seconds since _time_offset (module constant)
        # 
        # Time designation is slightly different bewteen data products; examples :
        #   MetO-NWS-BIO-dm-DOXY.h:		    time:units = "seconds since 2017-06-27 00:00:00" ;
        #   MetO-NWS-PHYS-hi-TEM.h:	  	    time:units = "seconds since 2015-07-03 00:00:00" ;
        #   MetO-NWS-REAN-PHYS-daily-TEM.h: time:units = "days since 1985-01-01" ;
        #   wave_BS_2010_1h.h:              time:units = "hours since 1900-01-01 00:00:00.0" 
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
        dtsec0 = int(0.5 + (toff - _time_offset).total_seconds()) # seconds from internal offset to data time offset 
        t = asarray(self.ncf.variables["time"][:], int64)         # int32 is not sufficient when converted to seconds
        self.timehash = asarray(0.5 + tsca*t + dtsec0, int64)     # round to long int64 vector
        #
        
            
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
        ## Extract basic 2D descriptors from open netCDF file ncf (setting them as attributes on self)
        #  
        #  capture and flag if lon/lat axes in data are reversed (flip_lon, flip_lat)
        #  @param self       open netCDF file object
        #  @return           None
        #  -----------------------------------------------------
        for axname in _lon_names:
            if self.ncf.dimensions.has_key(axname):
                self.nx         = self.ncf.dimensions[axname].size # integer
                lonvec     = self.ncf.variables[axname][:]    # dimension and axis name must be identical
                if lonvec[1] > lonvec[0]:                     # normal axis orientation
                    self.dlon       = lonvec[1]-lonvec[0]
                    self.lon0       = lonvec[0]
                    self.flip_lon   = False
                else:                                         # reversed axis orientation
                    self.dlon       = lonvec[0]-lonvec[1]
                    self.lon0       = lonvec[-1]
                    self.flip_lon   = True
                break
        else:
            raise exceptions.InputError("lon axis not recognized")
        #
        for axname in _lat_names:
            if self.ncf.dimensions.has_key(axname):
                self.ny         = self.ncf.dimensions[axname].size # integer
                latvec     = self.ncf.variables[axname][:]    # dimension and axis name must be identical
                if latvec[1] > latvec[0]:                     # normal axis orientation
                    self.dlat       = latvec[1]-latvec[0]
                    self.lat0       = latvec[0]
                    self.flip_lat   = False
                else:                                         # reversed axis orientation
                    self.dlat       = latvec[0]-latvec[1]
                    self.lat0       = latvec[-1]
                    self.flip_lat   = True
                break
        else:
            raise exceptions.InputError("lat axis not recognized")
        #

        
    def _extract_3D_cellw0(self, diagvarname):
        #  ------------------------------------------------------------------------
        ## Extract cell widths cellw0(nx,ny,nz) from open netCDF file ncf
        #  @param self           open netCDF file
        #  @param diagvarname    name of variable / variable name table, from which to detect topography
        #  @return               cellw0(nx,ny,nz)  (array of cell widths, -1 for dry cells)
        #
        #  assume cellw0 is time independent (pick time frame 0)
        #  Assess topography by probing for mask==True at data variable / table entry diagvarname in ncf.
        #  if diagvarname == None is passed, a KeyError will be raised
        #  If lookup fails an exception is raised 
        #
        #  Data variables are apparently loaded as masked arrays by
        #  the netcdf4 API. Therefore detect topography by testing mask==True (dry)
        #  ------------------------------------------------------------------------
        self.nz  = self.ncf.dimensions["depth"].size # integer
        #
        # ------ identify topography probing variable ------
        #
        if isinstance(diagvarname, dict): # search a pre-coded catalogue - pick first match
            for vnam in diagvarname.values():
                if vnam in self.ncf.variables.keys(): # found a match in pre-coded catalogue
                    break
            else:
                msg = "_extract_3D_cellw0: netCDF file %s does not contain a known data variable" % self.ncf.filepath()
                raise exceptions.LookupError(msg)
        else: # a specific variable name was provided
            vnam = diagvarname
            if not self.ncf.variables.has_key(vnam):
                msg = "_extract_3D_cellw0: netCDF file %s does not have variable %s" % (self.ncf.filepath(), str(vnam))
                raise exceptions.LookupError(msg)
            
        var           = self.ncf.variables[vnam][0,:]     # topography test variable, asked array, shape == (time, depth, lat, lon)
        var           = swapaxes(var, 0, 2)          # now shape = (lon, lat, depth)
        #
        ccz           = self.ncf.variables["depth"][:]    # layer centers
        cellw0        = zeros((self.nx, self.ny, self.nz), float)
        cellw0[:,:,:] = _cc2lw(ccz)                  # convert cell center positions to layer widths
        cellw0        = where(var.mask, -1, cellw0)  # dry cells assigned -1
        return cellw0

    def get_time_frame(self, vname, itime=0):
        ## Extract frame itime of variable vname and swap axes, corresponding to GWD data layout
        #  3D: (time, depth, lat, lon) ->  (lon, lat, depth)
        #  2D: (time, lat, lon)        ->  (lon, lat)
        #  netcdf4 module internally handles integer affine decompression, if data stored as integers
        # ------------------------------------------------------------------------------------------
        ncvar   = self.ncf.variables[vname]
        maskvar = ncvar[itime,:]
        var     = maskvar.data          # strip mask
        if  (len(ncvar.dimensions) == 4) and \
            (ncvar.dimensions[2] in _lat_names) and \
            (ncvar.dimensions[3] in _lon_names): # 3D case
            var   = swapaxes(var, 0, 2)                         # cast to 3D GWD layout
        elif  (len(ncvar.dimensions) == 3) and \
              (ncvar.dimensions[1] in _lat_names) and \
              (ncvar.dimensions[2] in _lon_names):  # 2D case
            var   = swapaxes(var, 0, 1)                         # cast to 2D GWD layout
        else:
            msg = "get_time_frame: unexpected axes / axes ordering : %s" % str(ncvar.dimensions)
            raise exception.ValueError(msg)
        # --- flip axes, if required, to obtain standard layout ---
        #
        if self.flip_lon:
            var = flip(var, axis=0)  # flip requires numpy >= 1.12
        if self.flip_lat:
            var = flip(var, axis=1)  # flip requires numpy >= 1.12
        #
        if _verbose:
            print "get_time_frame: exporting frame %d" % itime
        return var            
        
    def get_time(self, itime=0):
        ## return datetime corresponding to frame itime 
        return _to_datetime(self.timehash[itime])
    
    def get_number_of_frames(self, vname):
        return len(self.ncf.variables[vname])

    def get_wetmask(self, vname):
        ## extract wetmask from variable vname
        #  assume wetmask is time independent (pick time frame 0)
        #  3D: (time, depth, lat, lon) ->  (lon, lat, depth)
        #  2D: (time, lat, lon)        ->  (lon, lat)
        ncvar   = self.ncf.variables[vname]
        wetmask = where(ncvar[0,:].mask, 0, 1)   # pick time frame 0
        if ncvar.dimensions == ('time', 'depth', 'lat', 'lon'): # 3D case
            wetmask   = swapaxes(wetmask, 0, 2)                         # cast to 3D GWD layout
        elif ncvar.dimensions == ('time', 'lat', 'lon'):        # 2D case
            wetmask   = swapaxes(wetmask, 0, 1)                         # cast to 2D GWD layout
        else:
            msg = "get_wetmask: unexpected axes / axes ordering : %s" % str(ncvar.dimensions)
            raise exception.ValueError(msg)
        # --- flip axes, if required, to obtain standard layout ---
        if self.flip_lon:
            wetmask = flip(wetmask, axis=0)  # flip requires numpy >= 1.12
        if self.flip_lat:
            wetmask = flip(wetmask, axis=1)  # flip requires numpy >= 1.12
        #
        return wetmask
        
#################### self test #########################
#
#  mid North Sea test point: ix,iy = 210,240
#
if __name__ == "__main__":
    _verbose = True
    raw = CMEMS_DataSet("data/MetO-NWS-PHYS-hi-TEM.nc")
    v = raw.get_time_frame("votemper", 3)
    print v[210,240,:]
    print raw.get_wetmask("votemper")[210,240,:]
    
    
    
                    
         
