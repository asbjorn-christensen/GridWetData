#!/usr/bin/env python

from NWS_grids             import *
from GridWetData.grid_data import *

_verbose          = False  # control debugging output

#  ====================================================================================
## Super class for interpolation in space+time by first interpolating space, then time
#  specific for 3D + interaction with DataManager
# <b> Roles </b> 
#    \li keep track of data cache, decide whether to reload
#    \li interact with data manager
#    \li call hook is redirectd to interpolate, so instances are callable as obj(where, when)
#
# compositional relation to GridData_3D
# Attributes set by sub classes:
#       prop: data property type (set by sub classes)
#    
#    
class NWS_GridData_3DwithTime(GridData_3DwithTime):
    #  -----------------------------------------------------
    ## constructor
    #  Data is first loaded at interpolation time, when time argument is available
    #  @param self      The object pointer
    #  @param fname     netcdf data set file name  
    #  @param varname   variable name to look for (optional)
    #                   If absent (and class attributre prop is not set) the data set
    #                   is scanned for known variable names
    def __init__(self, fname, varname=None):
        self.grid   = NWSGrid_3D(fname)
        self.fname  = fname
        self.ncfile = NWS_dataformat.NWSDataSet(fname)
        # resolve property attribute prop
        if varname is None:
            if hasattr(self, "prop"):
                assert self.ncfile.has_variable(self.prop) # check solicited data is present
            else: # look for known variable names - pick first match
                for prop in NWS_dataformat.NWS_variablenames_3D.values():
                    if self.ncfile.has_variable(prop):
                        self.prop = prop
                        break
                else: 
                    msg  = "Unable to auto-locate data variable in data set %s" % fname
                    raise exceptions.ValueError(msg)
        else: # an explicit varname has bben provided
            self.prop = varname  
            assert self.ncfile.has_variable(varname) # check solicited data is present
            
        ## datetime instance corresponding to current cache content (None for empty cache)
        self.when_last = None # empty cache
        
    def __del__(self):
        try:
            self.ncfile.close()
        except:
            pass
      
        
    ## --------------------------------------------------------------------
    ## Generate a informative string representation of object
    #  @param self The object pointer
    #  @return string representation of instance
    #
    def __str__(self):
        iam = "GridData_3DwithTime instance for prop=%s for gtag=%s\n" % (self.prop, self.gtag)
        if self.when_last is None:
            iam = iam + "cache empty"
        else:
            iam = iam + "cache loaded for t=%s\n" % str(self.when_last)
            iam = iam + "left  bracket : w=%f file=%s\n" % (self.w0, self.gdata0.file_name)
            iam = iam + "right bracket : w=%f file=%s"   % (self.w1, self.gdata1.file_name)
        return iam

    ## Trigger load of a specific time frame as GridData_3D instance and tag instance with frame number
    #  @param self       The object pointer
    #  @param ifr        time frame to load
    #  @return           GridData_3D instance with loaded data, with attributes time_frame == ifr and
    #                    time corresponding to datetime of time frame ifr
    #
    def load_frame(self, ifr):
        if _verbose: print "load_frame: loading frame %d" % ifr
        data = self.ncfile.get_time_frame(self.prop, ifr)
        g3D  = GridData_3D(self.grid, data) # do not copy grid, which is static
        g3D.time_frame = ifr
        g3D.time       = self.ncfile.get_time(ifr)
        return g3D
        
    # --------------------------------------------------------------
    ## Load data corresponding to datetime instance argument when, if necessary, to cache
    #  In many practical cases, data need not to be reloaded, weights should just be changed or pointers switched
    #  Cache identity is decided time frame numbers 
    #  @param self     The object pointer
    #  @param when     datetime instance, which cache should be correspond to
    def update_cache(self, when):   
        bracket = self.ncfile.get_time_interpolation_bracket(when)
        if bracket is None:
            msg = "update_cache: when = %s is not interior for data set %s" % (str(when), self.fname)
        else:
            (i0,w0), (i1,w1) = bracket
        ## interpolation weight corresponding to lower time point
        self.w0 = w0 # may have changed, even though no reload is needed
        ## interpolation weight corresponding to upper time point
        self.w1 = w1 # may have changed, even though no reload is needed
        # ------- check, if we already have needed data in cache
        if hasattr(self, "gdata0") and hasattr(self, "gdata1"):
            old0 = self.gdata0
            old1 = self.gdata1
            # -- update self.gdata0
            if   i0 == old0.time_frame: 
                self.gdata0 = old0
            elif i0 == old1.time_frame:
                self.gdata0 = old1
            else:
                self.gdata0 = self.load_frame(i0)
            # -- update self.gdata1
            if   i1 == old0.time_frame: 
                self.gdata1 = old0
            elif i1 == old1.time_frame:
                self.gdata1 = old1
            else:
                self.gdata1 = self.load_frame(i1)
        else: # ------- first load, no cache
            ## GridData_3D instance corresponding to lower time point
            self.gdata0 = self.load_frame(i0)
            ## GridData_3D instance corresponding to upper time point
            self.gdata1 = self.load_frame(i1) 
        self.when_last = when

        

##########################################################################################
#    Specific elementary property classes
##########################################################################################

# The stagger method will trigger that load_frame restagger data to
# cell centered data after load. If class does not have this method, no
# restaggering is performed
   
## GridData_3D sub class for temperature data in both space and time
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
                                       


# TODO currents ()                        
#      "east_current"            : "vozocrtx",
#      "north_Current"           : "vomecrty",

                        
#################### self test #########################
#
#  mid North Sea test point: ix,iy = 100,150  <->  lon,lat = 4.2083, 56.025
#
if __name__ == "__main__":
    temp  = NWS_Temperature_3DwithTime("../../../tmp/MetO-NWS-PHYS-hi-TEM.nc")
    temp0 = temp.load_frame(0)
    pos   = [4.2083, 56.025, 0]
    ## --- test vertical get_vertex_indices
    maxwd = temp.grid.interpolate_wdepth(pos)
    for wd in arange(0, maxwd, 0.01):
        pos[2] = wd
        #print wd, temp0(pos)
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
    now  = datetime(2015, 7, 3, 2, 23)
    for i in range(120):
        print i, temp(pos, now)
        now = now + timedelta(seconds=60)
