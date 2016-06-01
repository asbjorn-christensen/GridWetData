#!/usr/bin/env python
# -*- coding: utf-8 -*-  
##########################################################################################
#        Data spectator example, deriving from grids / grid_data
#
#        Input data not shipped along with this module
#        git clone https://github.com/asbjorn-christensen/GridWetData
#        ncdump -h /home/data/ECMWF/NorthSea_waves/swh_NS_2004-2015.nc
#
##########################################################################################

from grids import * 
from grid_data import GridData_2D, GridData_withTime
from datetime import *

def interpolation_bracket(x, xlist):
    # constant extrapolation at left/right endpoints
    nt = len(xlist) 
    i1  = searchsorted(xlist, x)
    i0  = i1-1
    i1 = min(nt-1, max(i1, 0))
    i0 = min(nt-1, max(i0, 0))
    dx = xlist[i1] - xlist[i0] 
    if dx > 0:
        w0 = (xlist[i1]-x)/dx
    elif i0 == 0:
        w0 = 1.0   # left end
    elif i1 == nt-1:
        w0 = 0.0   # right end
    else:
        w0 = 0.5   # degenerate xlist with one element
    w1 = 1.0 - w0    
    return i0,w0,i1,w1


def in_days(dt):
    return 1.0*dt.days + dt.seconds/8.64e4 + dt.microseconds/8.64e10


# ===============================================================================
##  Provride interpolation for a specific 2D+time data set for surface waves
# 
#   Example spectator class, deriving from GridData_withTime
#   This class must provide update_cache and __init__
#   
# ===============================================================================
class WaveMovie(GridData_withTime):  #   
    
    ## -----------------------------------------------------
    ## constructor
    #  @param self   The object pointer.
    #  @param fname  file name to netCDF data set
    ## -----------------------------------------------------
    def __init__(self, fname):
        self.ncfile     = netcdf.NetCDFFile(fname, mmap=False) # keep open
        self.swh        = self.ncfile.variables["swh"]
        self.swh_scale  = self.ncfile.variables["swh"].scale_factor
        self.swh_offset = self.ncfile.variables["swh"].add_offset
        self.swh_fill   = self.ncfile.variables["swh"]._FillValue
        # ------------ time ------------
        self.time       = self.ncfile.variables["time"][:] 
        self.time_offset = datetime(1900,1,1) # time:units = "hours since 1900-01-01 00:00:0.0" ;
        # ------------ grids ------------
        lon  = self.ncfile.variables["longitude"][:]
        lat  = self.ncfile.variables["latitude"][:]
        dlon = lon[1] - lon[0]
        dlat = lat[1] - lat[0]
        self.grid      = LonLatGrid(len(lon), len(lat), lon[0], lat[0], dlon, dlat)
        self.when_last = None # empty cache
        #
    ## -----------------------------------------------------
    ## deconstructor
    #  @param self   The object pointer.
    #  close the netCDF file associated with the object
    ## -----------------------------------------------------    
    def __del__(self):
        if hasattr(self, "ncfile"): self.ncfile.close()

    ## -----------------------------------------------------
    ## Fetch a specific frame from data set - internal subroutine
    #  @param self   The object pointer.
    #  @param ifr    which frame to fetch from data set
    ## -----------------------------------------------------        
    def get_data_frame(self, ifr):
        # avoid transforming entire array
        return self.swh_offset + self.swh[ifr]*self.swh_scale

    ## -----------------------------------------------------
    ## Update cache corresponding to time when
    #  This method prepares time interpolation by instantiating grid data objects (gdata0, gdata1)
    #  at corresponding to neighboring time points. Attributes (w0,w1) are weight that should be used for
    #  linear time interpolation.
    #  Try to avoid I/O overhead if necessary, by analyzing whether needed data is already available in cache
    #  @param self   The object pointer.
    #  @param when   datetime object for time, which cache should correspond to
    ## -----------------------------------------------------   
    def update_cache(self, when):  
        dt    = when - self.time_offset # timedelta instance: Only days, seconds and microseconds are stored internally
        hours = 24.0*dt.days + dt.seconds/3.6e3 + dt.microseconds/3.6e9
        i0,w0,i1,w1 = interpolation_bracket(hours, self.time)
        ## interpolation weight corresponding to lower time point
        self.w0 = w0 # may have changed, even though no reload is needed
        ## interpolation weight corresponding to upper time point
        self.w1 = w1 # may have changed, even though no reload is needed
        # ------- check, if we already have needed data in cache
        if hasattr(self, "gdata0") and hasattr(self, "gdata1"):
            old0 = self.gdata0
            old1 = self.gdata1
            # -- update self.gdata0
            if   i0 == old0.index:
                self.gdata0 = old0
            elif i0 == old1.index:
                self.gdata0 = old1
            else:
                self.gdata0 = GridData_2D(self.grid, self.get_data_frame(i0))
                self.gdata0.index = i0
            # -- update self.gdata1
            if   i1 == old0.index:
                self.gdata1 = old0
            elif i1 == old1.index:
                self.gdata1 = old1
            else:
                self.gdata1 = GridData_2D(self.grid, self.get_data_frame(i1))
                self.gdata1.index = i1       
        else: # ------- first load, no cache
            ## GridData_2D instance corresponding to lower time point
            self.gdata0 = GridData_2D(self.grid, self.get_data_frame(i0))
            self.gdata0.index = i0
            ## GridData_2D instance corresponding to upper time point
            self.gdata1 = GridData_2D(self.grid, self.get_data_frame(i1))
            self.gdata1.index = i1
        #
        self.when_last = when   
        #

#################### self test #########################
#
#  mid North Sea test point: ix,iy = 100,150
#
if __name__ == "__main__":
    _verbose = True
    wmo = WaveMovie("/home/data/ECMWF/NorthSea_waves/swh_NS_2004-2015.nc")  # Input data not shipped along 
    start_time = datetime(2006,11,20, 5,23)
    end_time   = datetime(2006,11,22, 5,23)
    here       = (4.02, 55.07)
    dt         = timedelta(seconds=900)
    # -------------- test point --------------
    # print "SWH @ (%f,%f) %s = %f" % (here[0],here[1], str(start_time), wmo(here, start_time))
    # -------------- time scan --------------
    #test_time = start_time
    #while test_time < end_time:
    #    print in_days(test_time-start_time),  wmo(here, test_time)
    #    test_time = test_time + dt
    # -------------- space scan --------------
    for lon in arange(4.0, 6.0, 0.01):
        here       = (lon, 55.07)
        print lon, wmo(here, start_time)
    #
