#!/usr/bin/env python
# -*- coding: utf-8 -*-  
##########################################################################################
#        @package grid_data
#        Light weight generic data grid classes, which are independed on grid specifics
#
##########################################################################################
import data_manager
from grids import * 

_verbose          = False  # module printing level
_default_grid_tag = "ns"   # spatial region / resolution tag associated with data


## informal, but accepted, test for array 
def has_rank_2(obj):
    return hasattr(obj, "__len__") and hasattr(obj[0], "__len__") and not hasattr(obj[0][0], "__len__")


#  ====================================================================================
## Super class for interpolation in space+time by first interpolating space, then time - Generic for 2D/3D
#
#  sub classes must provide constructor and update_cache, which sets interpolation points (gdata0, gdata1) and weights (w0,w1)
#  corresponding to time when_last
        
class GridData_withTime:
    ## --------------------------------------------------------------------
    ## Space+time interpolation front end
    #  This methods parses time argument and performs time interpolation part after
    #  space interpolation which is delegated to grid_data methods
    #  Cache last time frames loade, to enable interpolation speedup
    #  @param self  The object pointer
    #  @param where single position / array of any shap of positions
    #  @param when  datetime instance corresponding to when data should be interpolated
    #  @return single/array of interpolated values
    #
    def interpolate(self, where, when):
        # 1) check cache to see, if we can avoid reloading data
        #    cache assumed valid if when correspond to last data fetch
        if self.when_last != when: # load new frame set
            self.update_cache(when)
        res0 = self.gdata0.interpolate(where) # scalar/array
        res1 = self.gdata1.interpolate(where) # scalar/array
        return self.w0*res0 + self.w1*res1    # inherit shape from res0,res1

    ## redirect call hook in class scope   
    __call__ = interpolate 
    #


#  ====================================================================================
## Super class for interpolation in space+time by first interpolating space, then time - specific for 3D
#
#  compositional relation to GridData_3D
class GridData_3DwithTime(GridData_withTime):
    ## --------------------------------------------------------------------------------------------
    ## Create a GridData_3D instance data snapshot corresponding to time when by time interpolation
    #  sea level corresponding to time when is also obtained by interpolation
    #
    def snapshot_at(self, when):
        if self.when_last != when: # load new frame set
            self.update_cache(when)
        mid_data = self.w0*self.gdata0.data + self.w1*self.gdata1.data  # assume array like 
        mid_grid = copy.copy(self.gdata0.grid)
        z0 = self.gdata0.grid.get_reference_level()
        z1 = self.gdata1.grid.get_reference_level()
        mid_grid.set_reference_level(self.w0*z0 + self.w1*z1)
        return GridData_3D(mid_grid, mid_data)


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
    #  @return GridData_3D instyance with loaded data
    #
    def load_frame(self, fname_data, fname_z):
        if _verbose: print "load_frame: %s" % fname_data
        newgrid = copy.copy(self.dmg.grids[self.gtag]) # shallow copy of ancestor grid 
        newgrid.update_sea_level(fname_z)
        return GridData_3D(newgrid, fname_data)
        
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


        
# ==============================================================================
## Super class for managing a grid data set corresponding to a specific time
#  grid is assumed corresponding to fname. Class is not aware about properties
#  of grid data nor properties of the grid nor the dimensionality of the grid.
#  The class is just a composition of a grid and data associated with the grid.
#  Data/grid/dimension specific issues are delegated to sub classes
#
# <b> Roles </b> 
#    \li  hold loaded data
#    \li  delegate interpolation request, depending on single/array argument
#     \li call hook is redirectd to interpolate, so instances are callable as obj(where)
#
class GridData:
    #  -----------------------------------------------------
    ## Dual constructor of single property ocean data snapshot
    #  Sea level is assumed updated in the passed grid (if it is a 3D grid type),
    #  grid data inputarg is either a suitable array (with shape consistent with grid)
    #  or a file name (if argument is a string type), which is passed to the grid data loader (load_data_frame)
    #  It is a user responsibility to ensure consistence between passed grid and provided data / data file
    #  @param self      The object pointer.
    #  @param grid      grid instance corresponding to data in fname
    #  @param inputarg  either a file name (if string type) of data, or a suitable data array (if not string type), corresponding to grid
    #  
    def __init__(self, grid, inputarg):
        ## grid instance corresponding to data
        self.grid      = grid  # do not copy
        #
        #     split behavior depending on inputarg
        #     test wheter inputarg argument is a string type - if so, assume
        #     it is a file name and pass it to the grid data loader
        #     otherwise assume it is a suitable data array ()
        #     Only set file_name, if provided
        #
        if isinstance(inputarg, basestring): # assume inputarg is a file name    
            ## data buffer
            self.data      = self.grid.load_data_frame(inputarg)
            ## file name corresponding to data (if provided, other wise not set)
            self.file_name = inputarg
        else:  # do not set file_name attribute
            ## data buffer
            self.data      = inputarg # assume it is a suitable data instance

    #  -----------------------------------------------------
    ## interpolate single space point pos
    #  @param  self   The object pointer.
    #  @param  pos  any sequence (lon, lat, at_depth)
    #  @return interpolated value
    #
    def _interpolate_single(self, pos):
        return self.grid.interpolate_data(self.data, pos)

    #  -----------------------------------------------------
    ## interpolate flat array of posistions pos[n,3]
    #  @param  self The object pointer.
    #  @param  pos  any nested sequence pos[n,3] of n positions, where a position correspond to (lon, lat, at_depth)
    #  @return array of interpolated values
    #
    def _interpolate_array(self, allpos):
        # --- pos must respond to numpy protocol
        buffer  = zeros(len(allpos), float) 
        for i,xyz in enumerate(allpos):
            buffer[i] = self._interpolate_single(xyz)
        return buffer

    #  -----------------------------------------------------
    ## spatial interpolation frontend, that delegates
    #  to vector/array interpoaltion
    #  @param  self The object pointer.
    #  @param  pos  any nested sequence (pos[n,3] of n positions (lon, lat, at_depth), or a single position as (lon, lat, at_depth)
    #  @return single/array of interpolated values
    #
    def interpolate(self, pos):
        if  has_rank_2(pos):
            return self._interpolate_array(pos)
        else:
            return self._interpolate_single(pos)

    ## redirect call hook in class scope
    __call__ = interpolate
    
    #  ============== generic writers ==============
    
    ## Write layer front end
    #  Delegate write request to grid, depending on filename extension 
    #  @param self   The object pointer.
    #  @param fname  File name for writing data. Resolve writing format from fname extension
    #  @param kwargs other optional arguments (parsed on format basis)
    def write_data(self, fname, **kwargs):
        (root, ext) = os.path.splitext(fname)
        if ext.lower() == ".nc":
            self.grid.write_data_as_netCDF(fname, self.data, **kwargs) # do not parse kwargs
        else:
            raise exception.ValueError("Filename extension in %s does not map to a writing method" % fname)


# ==============================================================================
## Sub class of GridData for 2-dimensional situations (currently void)
class GridData_2D(GridData):
    pass
    

# ==============================================================================
## Sub class of GridData for 3-dimensional situations
class GridData_3D(GridData):
    
    #  ============== Generic projectors ==============
    #  Project 3D data onto GridData_2D instances
    #  what surface means is defined by the resolution and type of the grid attribute
    #  to vector/array interpoaltion
    #  @param  self The object pointer.
    #
    #  -------------------------------------------------------
    ## Project surface layer of data at native grid resolution
    #  @param  self The object pointer.
    #  @return GridData_2D instance of projected surface layer
    #
    def get_surface_layer(self):
        grid, data = self.grid.get_surface_layer(self.data)
        return GridData_2D(grid, data)
    
    #  -------------------------------------------------------
    ## Project bottom layer of data at native grid resolution
    #  @param  self The object pointer.
    #  @return GridData_2D instance of projected bottom layer
    #
    def get_bottom_layer(self):
        grid, data = self.grid.get_bottom_layer(self.data)
        return GridData_2D(grid, data)
    
    #  -------------------------------------------------------
    ## Generate vertical average of data at native grid resolution
    #  @param  self The object pointer.
    #  @return GridData_2D instance of vertical average layer
    #
    def get_vertical_average(self):
        grid, data = self.grid.get_vertical_average(self.data)
        return GridData_2D(grid, data)

    
##########################################################################################
#    Specific property classes
##########################################################################################        

   

## GridData_3D sub class for temperature data in both space and time
class TemperatureData_3DwithTime(GridData_3DwithTime_DataManager):
    ## data property type 
    prop = "t"   # set class attribute, passed to data manager when soliciting data
    


#################### self test #########################
#
#  mid North Sea test point: ix,iy = 100,150  <->  lon,lat = 4.2083, 56.025
#
if __name__ == "__main__":
    g3D = DMIGrid_3D(os.path.join("DMI_data", "ns_grid.nc"))
    pick = os.path.join("DMI_data", "%s_ns_2015_05_20_00_15_00.nc")
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
    ##
    ## dmg = data_manager.DataManager("DMI_data", DMIGrid_3D)
    ## tt = TemperatureData_3DwithTime(dmg)
    ## pos = [4.2083, 56.025, 1.2]
    ## tt.interpolate(pos, datetime(2015, 05, 20, 00, 16, 00))
    ## print tt
    ## tt.interpolate(pos, datetime(2015, 05, 20, 00, 26, 00))
    ## print tt
    ## tt.interpolate(pos, datetime(2015, 05, 20, 00, 36, 00))
    ## print tt
    
