#!/usr/bin/env python

from HBM_grids             import *
from HBM_data_manager      import * 
from GridWetData.grid_data import *



_default_grid_tag = "ns"   # spatial region / resolution tag associated with data
_verbose          = False  # control debugging output



#  ====================================================================================
#  Narrow sub classes to GridData_3D for stepping through staggered data
#  Experimental solution that may be liftet up to base folder at some point

#  ====================================================================================
class GridData_3D_EastStaggerd(GridData_3D):
    #  -----------------------------------------------------
    ## constructor
    #  Handle data staggered on Eastern side of each cell
    #  @param self       The object pointer
    #  @param grid       grid object
    #  @param inputarg   filename / 3D array
    #  @param padval     pad value to define Western edge (optional)
    def __init__(self, grid, inputarg, padval=0.0):
        if isinstance(inputarg, basestring): # assume inputarg is a file name    
            data      = grid.load_data_frame(inputarg)
            self.file_name = inputarg
        else:  # do not set file_name attribute
            data  = inputarg # assume it is a suitable data instance

        shifted = zeros(data.shape, float)
        shifted[1:,:,:] = 0.5*(data[:-1,:,:] + data[1:,:,:]) # shift data to center point
        shifted[0, :,:] = 0.5*(padval        + data[ 1,:,:]) # use pad value to extrapolate to Western edge
        GridData_3D.__init__(self, grid, shifted)

class GridData_3D_SouthStaggerd(GridData_3D):
    #  -----------------------------------------------------
    ## constructor
    #  Handle data staggered on Southern side of each cell
    #  @param self       The object pointer
    #  @param grid       grid object
    #  @param inputarg   filename / 3D array
    #  @param padval     pad value to define Northern edge (optional)
    def __init__(self, grid, inputarg, padval=0.0):
        if isinstance(inputarg, basestring): # assume inputarg is a file name    
            data      = grid.load_data_frame(inputarg)
            self.file_name = inputarg
        else:  # do not set file_name attribute
            data  = inputarg # assume it is a suitable data instance
    
        shifted = zeros(data.shape, float)
        shifted[:,:-1,:] = 0.5*(data[:,:-1,:] + data[:,1:,:])  # shift data to center point
        shifted[:,-1, :] = 0.5*(padval        + data[:,-1,:])  # use pad value to extrapolate to Northern edge
        GridData_3D.__init__(self, grid, shifted)


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
class GridData_3DwithTime_DataManager(GridData_3DwithTime):
    #  -----------------------------------------------------
    ## constructor
    #  Data is first loaded at interpolation time, when time argument is available
    #  @param self   The object pointer
    #  @param dmg    data manager instance
    #  @param gtag   grid tag (spatial area + resolution). Default is _default_grid_tag
    #
    def __init__(self, dmg, gtag = _default_grid_tag):
        ## data manager instance (Role: map time to file names)
        self.dmg       = dmg
        ## datetime instance corresponding to current cache content (None for empty cache)
        self.when_last = None # empty cache
        ## the grid (spatial area / resolution) tag (set at instantiation)
        self.gtag      = gtag
        
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

    ## Trigger load of a specific time frame as GridData_3D instance and instance grid copy corresponding to that sea level
    #  It is assumed that provided data sets in (fname_data,fname_z) are consistent
    #  @param self       The object pointer
    #  @param fname_data file name of grid data for GridData_3D instantiation
    #  @param fname_z    file name of sea level elevation for grid update
    #  @return GridData_3D instance with loaded data
    #
    def load_frame(self, fname_data, fname_z):
        if _verbose: print "load_frame: %s" % fname_data
        newgrid = copy.copy(self.dmg.grids[self.gtag]) # shallow copy of ancestor grid 
        newgrid.update_sea_level(fname_z)
        g3D = GridData_3D(newgrid, fname_data)
        # stagger class method allows to recenter data
        if hasattr(self, "stagger"):
            self.stagger(g3D.data) # inplace recenter
        return g3D
        
    # --------------------------------------------------------------
    ## Load data corresponding to datetime instance argument when, if necessary, to cache
    #  In many practical cases, data need not to be reloaded, weights should just be changed or pointers switched
    #  Cache identity is decided by filenames of data content 
    #  Cache content is not validated beyond file name comparison
    #  @param self     The object pointer
    #  @param when     datetime instance, which cache should be correspond to
    def update_cache(self, when):    # verbose load
        (dname0, zname0, w0), (dname1, zname1, w1) = self.dmg.get_time_interpolation_files(self.prop, self.gtag, when)
        ## interpolation weight corresponding to lower time point
        self.w0 = w0 # may have changed, even though no reload is needed
        ## interpolation weight corresponding to upper time point
        self.w1 = w1 # may have changed, even though no reload is needed
        # ------- check, if we already have needed data in cache
        if hasattr(self, "gdata0") and hasattr(self, "gdata1"):
            old0 = self.gdata0
            old1 = self.gdata1
            # -- update self.gdata0
            if   dname0 == old0.file_name: 
                self.gdata0 = old0
            elif dname0 == old1.file_name:
                self.gdata0 = old1
            else:
                self.gdata0 = self.load_frame(dname0, zname0)
            # -- update self.gdata1
            if   dname1 == old0.file_name: 
                self.gdata1 = old0
            elif dname1 == old1.file_name:
                self.gdata1 = old1
            else:
                self.gdata1 = self.load_frame(dname1, zname1)
        else: # ------- first load, no cache
            ## GridData_3D instance corresponding to lower time point
            self.gdata0 = self.load_frame(dname0, zname0)
            ## GridData_3D instance corresponding to upper time point
            self.gdata1 = self.load_frame(dname1, zname1) 
        self.when_last = when




        

##########################################################################################
#    Specific elementary property classes
##########################################################################################

# The stagger method will trigger that load_frame restagger data to
# cell centered data after load. If class does not have this method, no
# restaggering is performed
   
## GridData_3D sub class for temperature data in both space and time
class TemperatureData_3DwithTime(GridData_3DwithTime_DataManager):
    ## data property type 
    prop = "t"   # set class attribute, passed to data manager when soliciting data

class SalinityData_3DwithTime(GridData_3DwithTime_DataManager):
    ## data property type 
    prop = "s"   # set class attribute, passed to data manager when soliciting data

## GridData_3D sub class for u currents in both space and time
#  currently make a global restaggering of data, so data corresponds to cell-centered data
class WECurrentData_3DwithTime(GridData_3DwithTime_DataManager):
    ## data property type 
    prop    = "u"           # set class attribute, passed to data manager when soliciting data
    def stagger(self, x):   # data staggered to east side of cell - assume padded with zeros
        x[:-1,:,:] = 0.5*(x[:-1,:,:] + x[1:,:,:]) # inplace

## GridData_3D sub class for u currents in both space and time
#  currently make a global restaggering of data, so data corresponds to cell-centered data
class NSCurrentData_3DwithTime(GridData_3DwithTime_DataManager):
    ## data property type 
    prop    = "v"           # set class attribute, passed to data manager when soliciting data
    def stagger(self, x):   # data staggered to south side of cell - assume padded with zeros
        x[:,:-1,:] = 0.5*(x[:,:-1,:] + x[:,1:,:]) # inplace

## Composition of two (WE,NS) CurrentData_3DwithTime instances
class HorizontalCurrentData_3DwithTime:
    def __init__(self, *args, **kwargs):
        self.u = WECurrentData_3DwithTime(*args, **kwargs)
        self.v = NSCurrentData_3DwithTime(*args, **kwargs)
    def interpolate(self, where, when):
        ures = self.u(where, when)
        vres = self.v(where, when)
        if  has_rank_2(where):
            return array(zip(ures,vres))
        else:
            return array([ures,vres])
    __call__ = interpolate



#################### self test #########################
#
#  mid North Sea test point: ix,iy = 100,150  <->  lon,lat = 4.2083, 56.025
#
if __name__ == "__main__":
    g3D = HBMGrid_3D(os.path.join("../DMI_data", "ns_grid.nc"))
    pick = os.path.join("../DMI_data", "%s_ns_2015_05_20_00_15_00.nc")
    g3D.update_sea_level(pick % "z")
    temp = GridData_3D(g3D, pick % "t")
    pos = [4.2083, 56.025, 0]
    
    ## for i in range(3):
    ##     pos[2] = g3D.ccdepth[100,150,i]
    ##     print g3D.ccdepth[100,150,i], temp.data[100,150,i]
    ## print g3D.get_vertex_indices(pos)
    ## --- test vertical get_vertex_indices
    ## maxwd = g3D.interpolate_wdepth(pos)
    ## for wd in arange(0, maxwd, 0.01):
    ##     pos[2] = wd
    ##     (ix,iy,iz,sx,sy,sz) = g3D.get_vertex_indices(pos)
    ##     print wd, iz,sz
    ##
    ## --- test vertical get_vertex_indices
    maxwd = g3D.interpolate_wdepth(pos)
    for wd in arange(0, maxwd, 0.01):
        pos[2] = wd
        print wd, temp(pos)
    ## import data_manager
    ## dmg = data_manager.DataManager("../DMI_data", HBMGrid_3D)
    ## tt = TemperatureData_3DwithTime(dmg)
    ## pos = [4.2083, 56.025, 1.2]
    ## tt.interpolate(pos, datetime(2015, 05, 20, 00, 16, 00))
    ## print tt
    ## tt.interpolate(pos, datetime(2015, 05, 20, 00, 26, 00))
    ## print tt
    ## tt.interpolate(pos, datetime(2015, 05, 20, 00, 36, 00))
    ## print tt
    
