#!/usr/bin/env python

from NWS_grids             import *
from GridWetData.grid_data import *

_verbose          = False  # control debugging output


#  ====================================================================================
## Provider Super class for interpolation in space+time by first interpolating space, then time
#  common for 2D/3D and scalar/vector data
# <b> Roles </b> 
#    \li keep track of data cache, decide whether to reload
#    \li interact with data manager
#    \li call hook is redirectd to interpolate, so instances are callable as obj(where, when)
#  
  
class NWS_GridData_withTime(GridData_withTime):
    #  -----------------------------------------------------
    def __del__(self):
        try:
            self.ncfdata.close()
        except:
            pass
        
    def close(self):
        self.ncfdata.close()
        
    ## Load data corresponding to datetime instance argument when, if necessary, to cache
    #  In many practical cases, data need not to be reloaded, weights should just be changed or pointers switched
    #  Cache identity is decided time frame numbers 
    #  @param self     The object pointer
    #  @param when     datetime instance, which cache should be correspond to
    def update_cache(self, when):   
        bracket = self.ncfdata.get_time_interpolation_bracket(when)
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
class NWS_GridData_3DwithTime(NWS_GridData_withTime, GridData_3DwithTime):
    #  -----------------------------------------------------
    ## constructor
    #  Data is first loaded at interpolation time, when time argument is available
    #  spectator class to view a netcdf data set with time dimension
    #  @param self      The object pointer
    #  @param fname     netcdf data set file name  
    #  @param varname   variable name to look for (optional)
    #                   If absent (and class attributre prop is not set) the data set
    #                   is scanned for known variable names
    def __init__(self, fname, varname=None):
        self.fname  = fname
        self.ncfdata = NWS_dataformat.NWSDataSet(fname)
        self.grid   = NWSGrid_3D(fname)
        # resolve property attribute prop
        if varname is None:
            if hasattr(self, "prop"):
                assert self.ncfdata.has_variable(self.prop) # check solicited data is present
            else: # look for known variable names - pick first match
                for prop in NWS_dataformat.NWS_variablenames_3D.values():
                    if self.ncfdata.has_variable(prop):
                        self.prop = prop
                        break
                else: 
                    msg  = "Unable to auto-locate data variable in data set %s" % fname
                    raise exceptions.ValueError(msg)
        else: # an explicit varname has been provided
            self.prop = varname  
            assert self.ncfdata.has_variable(varname) # check solicited data is present
            
        ## datetime instance corresponding to current cache content (None for empty cache)
        self.when_last = None # empty interpolation cache
        
    ## --------------------------------------------------------------------
    ## Generate a informative string representation of object
    #  @param self The object pointer
    #  @return string representation of instance
    #
    def __str__(self):
        iam = "NWS_GridData_3DwithTime instance for prop=%s\n" % self.prop
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
        data = self.ncfdata.get_time_frame(self.prop, ifr)
        g3D  = GridData_3D(self.grid, data) # do not copy grid, which is static
        g3D.time_frame = ifr
        g3D.time       = self.ncfdata.get_time(ifr)
        return g3D

    ## Inquire number of time frames in current data set for this variable
    #  @param self       The object pointer
    #
    def get_number_of_frames(self): 
        return self.ncfdata.get_number_of_frames(self.prop)
        
    

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
class NWS_GridVector_3DwithTime(NWS_GridData_withTime):
    #  -----------------------------------------------------
    ## constructor
    #  Data is first loaded at interpolation time, when time argument is available
    #  spectator class to view a netcdf data set with time dimension
    #  @param self       The object pointer
    #  @param fname      netcdf data set file name  
    #  @param varnames   variable names to look for (vector of strings) (optional)
    #                    If absent, class attributre prop must be set
    #                    No scanning for known variable names
    #                    for each token "zero_data", a zero-valued field is inserted in that position
    #                    (allow to e.g. uprank a horizontal field to 3D)
    def __init__(self, fname, varnames=None):
        self.fname   = fname
        self.ncfdata = NWS_dataformat.NWSDataSet(fname)
        self.grid    = NWSGrid_3D(fname)
        # resolve property attribute prop
        if varnames is None:
            if not hasattr(self, "prop"):
                msg  = "Unable to identify which variables should be loaded from data set %s" % fname
                raise exceptions.ValueError(msg)
        else: # an explicit varname has been provided
            self.prop = varnames
        for vnam in self.prop: # check solicited data is present
            assert self.ncfdata.has_variable(vnam) 
            
        ## datetime instance corresponding to current cache content (None for empty cache)
        self.when_last = None # empty interpolation cache
        
        
    ## --------------------------------------------------------------------
    ## Generate a informative string representation of object
    #  @param self The object pointer
    #  @return string representation of instance
    #
    def __str__(self):
        iam = "NWS_GridVector_3DwithTime instance for propa=%s\n" % "".join(self.prop)
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
        # --- compose data in vector according to prop as a plain list
        data = []
        for prop in self.prop:
            if prop is "zero_data": # allow to e.g. uprank a horizontal field to 3D
                data.append(zeros(self.grid.nx, self.grid.ny, self.grid.nz), float)  # grid is NWSGrid_3D
            else:
                data.append(self.ncfdata.get_time_frame(prop, ifr))
        g3D  = GridVector_3D(self.grid, data) # do not copy grid, which is static
        g3D.time_frame = ifr
        g3D.time       = self.ncfdata.get_time(ifr)
        return g3D

    ## Inquire number of time frames in current data set for this variable
    #  @param self       The object pointer
    #
    #  since all variables have same time dimension, just probe first variable 
    def get_number_of_frames(self): 
        return self.ncfdata.get_number_of_frames(self.prop[0]) # probe first variable 
        
# class NWS_GridVector_2DwithTime(NWS_GridData_withTime):  deferred until needed

#  ====================================================================================
## Super class for interpolation of 2D data in space+time by first interpolating space, then time
#  

class NWS_GridData_2DwithTime(NWS_GridData_withTime):
    #  -----------------------------------------------------
    ## constructor
    #  Data is first loaded at interpolation time, when time argument is available
    #  spectator class to view a netcdf data set with time dimension
    #  @param self      The object pointer
    #  @param fname     netcdf data set file name  
    #  @param varname   variable name to look for (optional)
    #                   If absent (and class attributre prop is not set) the data set
    #                   is scanned for known variable names
    def __init__(self, fname, varname=None, get_wetmask=False):
        self.fname  = fname
        self.ncfdata = NWS_dataformat.NWSDataSet(fname)
        self.grid   = NWSGrid_2D(fname)
        # resolve property attribute prop
        if varname is None:
            if hasattr(self, "prop"):
                assert self.ncfdata.has_variable(self.prop) # check solicited data is present
            else: # look for known variable names - pick first match
                for prop in NWS_dataformat.NWS_variablenames_2D.values():
                    if self.ncfdata.has_variable(prop):
                        self.prop = prop
                        break
                else: 
                    msg  = "Unable to auto-locate data variable in data set %s" % fname
                    raise exceptions.ValueError(msg)
        else: # an explicit varname has been provided
            self.prop = varname  
            assert self.ncfdata.has_variable(varname) # check solicited data is present

        if get_wetmask:
            self.grid.wetmask = self.ncfdata.get_wetmask(self.prop)
            
        ## datetime instance corresponding to current cache content (None for empty cache)
        self.when_last = None # empty interpolation cache
        
    ## --------------------------------------------------------------------
    ## Generate a informative string representation of object
    #  @param self The object pointer
    #  @return string representation of instance
    #
    def __str__(self):
        iam = "NWS_GridData_2DwithTime instance for prop=%s\n" % self.prop
        if self.when_last is None:
            iam = iam + "cache empty"
        else:
            iam = iam + "cache loaded for t=%s\n" % str(self.when_last)
            iam = iam + "left  bracket : w=%f file=%s\n" % (self.w0, self.gdata0.file_name)
            iam = iam + "right bracket : w=%f file=%s"   % (self.w1, self.gdata1.file_name)
        return iam

    ## Trigger load of a specific time frame as GridData_2D instance and tag instance with frame number
    #  @param self       The object pointer
    #  @param ifr        time frame to load
    #  @return           GridData_2D instance with loaded data, with attributes time_frame == ifr and
    #                    time corresponding to datetime of time frame ifr
    #
    def load_frame(self, ifr):
        if _verbose: print "load_frame: loading frame %d" % ifr
        data = self.ncfdata.get_time_frame(self.prop, ifr)
        g2D  = GridData_2D(self.grid, data) # do not copy grid, which is static
        g2D.time_frame = ifr
        g2D.time       = self.ncfdata.get_time(ifr)
        return g2D

    ## Inquire number of time frames in current data set for this variable
    #  @param self       The object pointer
    #
    def get_number_of_frames(self): 
        return self.ncfdata.get_number_of_frames(self.prop)
    

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
