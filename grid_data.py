#!/usr/bin/env python
# -*- coding: utf-8 -*-  
##########################################################################################
#        @package grid_data
#        Light weight generic data grid classes, which are independed on grid specifics
#
#        2D means horizontal
#        3D means horizontal+vertical
##########################################################################################

from grids import * 

_verbose          = False  # module printing level



## informal, but accepted, test for array 
def has_rank_2(obj):
    return hasattr(obj, "__len__") and hasattr(obj[0], "__len__") and not hasattr(obj[0][0], "__len__")


# =======================================================================================
## Super class for managing generic issues in relation to grid data / vector / ... frames corresponding to a specific time (or time invariant)
#  
# <b> Roles </b> 
#    \li manage generic issues in object frame on a grid
#
class GridObject:
     def pass_time_to_friend(self, other):
         ## pass on time attribute for self to other, if present
         if hasattr(self, "time"):
             other.time = self.time # reference, not copy
        
# ==============================================================================
## Super class for managing a grid (scalar) data set corresponding to a specific time
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
class GridData(GridObject):
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

    
    
    #  ============== whole array operations ==============
    ## gradient of field - delegated to grid
    def gradient(self):
        grad = self.grid.gradient(self.data)
        obj = GridVector(self.grid, grad)
        self.pass_time_to_friend(obj)
        return obj
    
    #  ============== generic writers ==============
    
    ## Write layer front end
    #  Delegate write request to grid, depending on filename extension (currently only netcdf hook)
    #  @param self   The object pointer.
    #  @param fname  File name for writing data. Resolve writing format from fname extension
    #  @param kwargs other optional arguments (parsed on format basis)
    def write_data(self, fname, **kwargs):
        (root, ext) = os.path.splitext(fname)
        if ext.lower() == ".nc":
            kwargs = copy.copy(kwargs) # do not mangle original arg
            if kwargs.has_key("netcfd_attr"):
                netcfd_attr = kwargs["netcfd_attr"]
                del kwargs["netcfd_attr"]
            else:
                netcfd_attr = {}
            self.grid.write_data_as_COARDS(fname, self.data, netcfd_attr, **kwargs) # do not parse kwargs
        else:
            raise exception.ValueError("Filename extension in %s does not map to a writing method" % fname)

    ## Specific writer for netcdf (AnchorLab variant) - delegate write request to grid
    #  @param self   The object pointer.
    #  @param fname  File name / file pointer
    #  @param kwargs other optional arguments (parsed on format basis)    
    def write_data_as_netCDF(self, fname,  **kwargs):
        self.grid.write_data_as_netCDF(fname, self.data, **kwargs) # do not parse kwargs


        
class GridVector(GridObject):
    #  -----------------------------------------------------
    ## Dual constructor of vector property ocean data snapshot
    #  All vector components pertain to same grid (i.e. no staggering for this class)
    #  Sea level is assumed updated in the passed grid (if it is a 3D grid type),
    #  grid data inputarg is either a suitable array (with shape consistent with grid, vector dimension is leading dimension)
    #  or a file name (if argument is a string type), which is passed to the grid data loader (load_vector_frame) 
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
            self.data      = self.grid.load_vector_frame(inputarg)
            ## file name corresponding to data (if provided, other wise not set)
            self.file_name = inputarg
        else:  # do not set file_name attribute
            ## data buffer
            self.data      = inputarg # assume it is a suitable data instance

    #  -----------------------------------------------------
    ## interpolate single space point pos
    #  @param  self   The object pointer.
    #  @param  pos  any sequence (lon, lat, at_depth)
    #  @return interpolated values
    #
    def _interpolate_single(self, pos):
        buffer = []
        for vdat in self.data:
            buffer.append(self.grid.interpolate_data(vdat, pos))
        return array(buffer, float)

    #  -----------------------------------------------------
    ## interpolate flat array of posistions pos[n,3]
    #  @param  self The object pointer.
    #  @param  pos  any nested sequence pos[n,3] of n positions, where a position correspond to (lon, lat, at_depth)
    #  @return array of interpolated values; shape = (n, vectordim)
    #
    def _interpolate_array(self, allpos):
        # --- pos must respond to numpy protocol
        buffer  = []
        for i,xyz in enumerate(allpos):
            buffer.append(self._interpolate_single(xyz))
        return array(buffer, float) # shape = (n, vectordim)

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
    #
    #  ============== whole array operations ==============
    
    
    
    #  ============== generic writers ==============
    #  currently no generic writers defined


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
    def get_surface_layer(self, **kwargs):
        grid, data = self.grid.get_surface_layer(self.data, **kwargs)
        obj = GridData_2D(grid, data)
        self.pass_time_to_friend(obj)
        return obj
    
    #  -------------------------------------------------------
    ## Project bottom layer of data at native grid resolution
    #  @param  self The object pointer.
    #  @return GridData_2D instance of projected bottom layer
    #
    def get_bottom_layer(self, **kwargs):
        grid, data = self.grid.get_bottom_layer(self.data, **kwargs)
        obj = GridData_2D(grid, data)
        self.pass_time_to_friend(obj)
        return obj
    
    #  -------------------------------------------------------
    ## Generate vertical average of data at native grid resolution
    #  @param  self The object pointer.
    #  @return GridData_2D instance of vertical average layer
    #
    def get_vertical_average(self):
        grid, data = self.grid.get_vertical_average(self.data)
        obj = GridData_2D(grid, data)
        self.pass_time_to_friend(obj)
        return obj

    #  -------------------------------------------------------
    ## Generate vertical max of data at native grid resolution
    #  @param  self The object pointer.
    #  @return GridData_2D instance of vertical max layer
    #
    def get_vertical_max(self):
        grid, data = self.grid.get_vertical_max(self.data)
        obj = GridData_2D(grid, data)
        self.pass_time_to_friend(obj)
        return obj

# ==============================================================================
## Sub class of GridData for 2-dimensional situations (currently void)
class GridVector_2D(GridVector):
    #  ============== whole array operations ==============
    ## horizontal curl of vector field - delegated to grid
    def curl(self):
        curl = self.grid.curl(self.data)
        obj = GridData_2D(self.grid, curl)  # only z component of the curl
        self.pass_time_to_friend(obj)
        return obj
    

# ==============================================================================
## Sub class of GridData for 3-dimensional situations
class GridVector_3D(GridVector):
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
    def get_surface_layer(self, **kwargs):
        data = []
        for component in self.data:
            grid, surfdata = self.grid.get_surface_layer(component, **kwargs)
            data.append(surfdata)
        obj = GridVector_2D(grid, data)
        self.pass_time_to_friend(obj)
        return obj

    ## full curl of vector field - delegated to grid
    def curl(self):
        curl = self.grid.curl(self.data)   # 3D grid -> 3D vector
        obj = GridVector_3D(self.grid, curl)
        self.pass_time_to_friend(obj)
        return obj
    
#  ====================================================================================
## Super class for offline interpolation in space+time by first interpolating space, then time - Generic for 2D/3D
#
#  sub classes must provide constructor and update_cache, which sets interpolation points (gdata0, gdata1) and weights (w0,w1)
#  corresponding to time when_last
class GridData_withTime:
    ## --------------------------------------------------------------------
    ## Space+time interpolation front end for scalar+vector data
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

