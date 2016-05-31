#!/usr/bin/env python
# -*- coding: utf-8 -*-  
##########################################################################################
#   Interface to DTU10 model to assess astronomical components of sea-surface elevation
#   and derived geostrophical surface currents
#   Requires module jdcal (https://pypi.python.org/pypi/jdcal) and module dtu10 
#   Linux setup procedure for DTU10: f2py -c -m dtu10API ortran_sources/perth3.f
#   f2py is part of SciPy
##########################################################################################

from grids import *
from grid_data import has_rank_2
from datetime import *
from subprocess import Popen, PIPE
import jdcal
from constants import EarthMeanRadius, deg2rad, g, omega
import dtu10API

# ===============================================================================
##  Provide a GridData_withTime like interface to the online model DTU10
#   TODO: capture range errors
#         dry point attempts
#         include summertime (DTU10 refers to UTC)
# ===============================================================================

class DTU10:  
    ## --------------------------------------------------------------------
    ## Space+time interpolation front end (GridData_withTime like)
    #  @param self  The object pointer
    #  @param where single position / array of any shape of positions
    #  @param when  datetime instance corresponding to when data should be interpolated
    #  @return single/array of interpolated values
    def interpolate(self, where, when):
        if  has_rank_2(where):  # assume array-like
            buff = zeros(len(where), float) 
            for i,pos in enumerate(where):
                buff[i] = self._interpolate_single(pos, when)
            return buff
        else:
            return self._interpolate_single(where, when)
        
    __call__       = interpolate
    _numdifflength = 100 # meters, for numerical differentiation
    
    ## --------------------------------------------------------------------
    ## Space+time interpolation at specific point and time
    #  @param self  The object pointer
    #  @param where single position 
    #  @param when  datetime instance corresponding to when data should be interpolated
    #  @return single/array of interpolated values (SSE in meters)
    #  Invoke perth3 subroutine in module dtu10API:
    #  perth3( DLAT, DLON, TIME, TIDE, ISDATA )
    #  Arguments -
    #     name      type  I/O               description
    #     ----      ----  ---               -----------
    #     DLAT       D     I    north latitude (in degrees) for desired
    #                           location.
    #
    #     DLON       D     I    east longitude (in degrees).
    #
    #     TIME       D     I    desired UTC time, in (decimal) Modified Julian Date
    #                           e.g., 1 Jan 1990 noon = 47892.5
    #
    #     TIDE       D     O    computed tidal height.  The units will be
    #                           the same as the 'amplitude' units on the
    #                           input tidal grids (usually cm).
    #
    #     ISDATA     L     O    logical denoting whether tide data exist at
    #                           desired location.  If FALSE, then TIDE is
    #                           not modified.
    def _interpolate_single(self, where, when):
        TIDE       = array(1.0, float)
        ISDATA     = array(False, bool)
        ttup       = when.timetuple()[:3]
        jday       = jdcal.gcal2jd(*ttup)[1]
        hfrac      = (when.hour + (when.minute + (when.second + when.microsecond*1e-6)/60.0)/60.0)/24.0
        lookuptime = jday+hfrac
        dummy = dtu10API.perth3( where[1], where[0], lookuptime, TIDE, ISDATA )
        print 
        if bool(ISDATA):
            return float(TIDE)/100.0 # convert to meters
        else:
            raise GridRangeViolation # perth3 does not resolve further
        
    ## --------------------------------------------------------------------
    ## Geostrophic surface currents at a specific point and time, obtained by finite differences
    #  @param self  The object pointer
    #  @param where single position / array of any shap of positions
    #  @param when  datetime instance corresponding to when data should be interpolated
    #  @param step  [meters] optional length for finite step differentiation
    #  By standard formula, see e.g. Robert H. Stewart, Introduction to Physical Oceanography, 2008, p.155
    def geostrophic_current(self, where, when, step = _numdifflength):
        LAT = where[1]
        LON = where[0]
        f   = 2*omega*sin(LAT*deg2rad) # Coriolis frequency at this latitude
        ttup       = when.timetuple()[:3]
        jday       = jdcal.gcal2jd(*ttup)[1]
        hfrac      = (when.hour + (when.minute + (when.second + when.microsecond*1e-6)/60.0)/60.0)/24.0
        lookuptime = jday+hfrac
        dlon   = 0.5*step/EarthMeanRadius/cos(LAT*deg2rad) # assume spherical earth
        dlat   = 0.5*step/EarthMeanRadius                  # assume spherical earth
        TIDE0  = array(1.0, float)
        TIDE1  = array(1.0, float)
        ISDATA = array(False, bool)
        #
        dummy = dtu10API.perth3( LAT, LON+dlon, lookuptime, TIDE1, ISDATA )
        if not bool(ISDATA):
            raise GridRangeViolation # perth3 does not resolve further
        dummy = dtu10API.perth3( LAT, LON-dlon, lookuptime, TIDE0, ISDATA )
        if not bool(ISDATA):
            raise GridRangeViolation # perth3 does not resolve further
        dtide_dx = float((TIDE1-TIDE0)/step) # centered difference
        #
        dummy = dtu10API.perth3( LAT+dlat, LON, lookuptime, TIDE1, ISDATA )
        if not bool(ISDATA):
            raise GridRangeViolation # perth3 does not resolve further
        dummy = dtu10API.perth3( LAT-dlat, LON, lookuptime, TIDE0, ISDATA )
        if not bool(ISDATA):
            raise GridRangeViolation # perth3 does not resolve further
        dtide_dy = float((TIDE1-TIDE0)/step) # centered difference
        #
        return g*array([-dtide_dy, dtide_dx], float)/f
       
    

#################### self test #########################
#
#  mid North Sea test point: ix,iy = 100,150
#
if __name__ == "__main__":
    _verbose = True
    sse = DTU10()
    start_time = datetime(2006,11,20, 5,23)
    end_time   = datetime(2006,11,22, 5,23)
    here       = (4.02, 55.07)
    dt         = timedelta(seconds=900)
    # -------------- test point --------------
    #print "SWH @ (%f,%f) %s = %f meters" % (here[0],here[1], str(start_time), sse(here, start_time))
    #print "SWH @ (%f,%f) %s = %f meters" % (here[0],here[1], str(end_time),   sse(here, end_time))
    # -------------- time scan --------------
    #test_time = start_time
    #while test_time < end_time:
    #    print in_days(test_time-start_time),  sse(here, test_time)
    #    test_time = test_time + dt
    # -------------- space scan --------------
    #for lon in arange(4.0, 6.0, 0.01):
    #    here       = (lon, 55.07)
    #    print lon, sse(here, start_time)
    #
    # -------------- geostrophical current time scan -------------
    test_time = start_time
    while test_time < end_time:
        print "%f %f" % tuple(sse.geostrophic_current(here, test_time))
        test_time = test_time + dt
