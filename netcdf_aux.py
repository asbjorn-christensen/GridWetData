#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################################################################
##       @package netcdf_aux
#        Auxillary functions for convenient netcdf I/O
#        
#######################################################################################################################
import os
import exceptions
from numpy import *
import numbers   # python 2.6+

_thisdir         = os.path.dirname(__file__)  # allow remote import
_verbose         = False # True  # log info

## Return the index of the first occurence of item in vec by == operator
#  where vec does not have an index method (like a numpy array)
#  @param item        item to search for
#  @param vec         listlike to search in (must have a length)
#  @return            index of item (if found) or None (if not found)
def find_first(item, vec):
    for i in range(len(vec)):
        if item == vec[i]:
            return i  
    return None # item not found



## Inject/Append data + index along unlimited dimension after last stored frame
#  currently index-name is hardcoded to "index"
#  @param ncfile        NetCDFFile instance
#  @param data          2D data to be written
#  @param index         index value for this data frame (assume index is a unique frame identifier)
#  @param dataParam     Optional keyword argument:  Custom name for data variable. Default: data
#  Assumes netcdf variable and index variable is already created
#  Resolve frame number to write to by checking whether index value exist (frame injection), otherwise append along unlimited dimension
def add_data_to_unlimted_var(ncfile, data, **kwargs):
    if kwargs.has_key("dataParam"):  # Default: data
        dname = kwargs["dataParam"]  # apply custom data variable name
    else:
        dname = "data"               # default data variable name
    index_value                     = kwargs["index"] # mandatory keyword argument
    # --- resolve frame writing position ---
    writepos = find_first(index_value, ncfile.variables["index"])
    if writepos is None:                                        # index_value unique ...
        writepos = ncfile.variables["index"].shape[0]           # ... therefore append frame
        ncfile.variables["index"][writepos]  = index_value      # and store new index_value
    #
    ncfile.variables[dname][writepos,:,:] = data                # inject/append data
    #


## Create NetCDF variable in open NetCDFFile
#  @param ncfile        NetCDFFile instance
#  @param data          2D data to be written
#  @param verbose       True: overwrite (if defined); False: just return existing variable (if defined)
#  @param specification Optional keyword argument:  specification title attribute associated with data. Default: None
#  @param dataParam     Optional keyword argument:  Custom name for data variable. Default: data
#  @param storage_type  Optional keyword argument:  netcdf API supports these types:          Default: select type corresponding
#
#  @return NetCDF variable 
#
#  Data is stored in netCDF variable "data(nx,ny)" in single frame mode and "data(nframes,nx,ny)" indexed by "index(fnframes)" in append mode
#  If dataParam is provided, this variable name is applied instead of default "data"
#  Assumes dimensions nx,ny (and possibly unlimited dimension nframes) are created
#  if verbose == False and variable is already defined, consistency with provided kwargs is not enforced, current definition just returned
#
def create_netCDF_xyvariable(ncfile, data, verbose, **kwargs):
    #
    if kwargs.has_key("dataParam"):  # Default: data
        dname = kwargs["dataParam"]  # apply custom data variable name
    else:
        dname = "data"               # default data variable name
    #
    if verbose == False and ncfile.variables.has_key(dname):  # check if already defined - if so, return existing variable
        return ncfile.variables[dname]
    #
    # -------- define/overwrite variable --------
    #
    if kwargs.has_key("specification"):  # Default: None
        data_specification = kwargs["specification"]
    else:
        data_specification = None
    #
    if kwargs.has_key("storage_type"):
        data_type = kwargs["storage_type"] # do not assess validity
    else:
        data_type = data.dtype # Default: apply array type, assume data has dtype attribute
    #
    #  --------- create netcdf variables ---------
    #
    index = kwargs["index"] # mandatory keyword argument
    if index is not None: # append mode
        var_data  = ncfile.createVariable(dname,  data_type, ('nframes','nx','ny')) 
    else:   #  single frame mode
        var_data = ncfile.createVariable(dname, data_type, ('nx','ny'))
    #
    if data_specification:
        var_data.specification = data_specification
    # 
    return var_data



## Write a NetCDF variable in open NetCDFFile according to COARDS format
#  reference :  http://ferret.pmel.noaa.gov/Ferret/documentation/coards-netcdf-conventions
#  @param ncfile        NetCDFFile instance
#  @param data          2D data to be written (numpy like). If data is None, just declare variable
#  @metadata            metadata for variable (mandatory: lon, lat, variable_name, units, long_name)
#  @writeopt            (optional) write options specified by keyword - data_layout must be specified
def write_lonlatdata_in_COARDS_format(ncfile, data, metadata, **writeopt):   
    # dump data(nlon,nlat) corresponding to data on a positively oriented,
    # uniformly spaced lon-lat grid to open netcdf file ncfile
    # data is transposed on writing, corresponding to COARDS axis requirements (order TZYX)
    # metadata is a dictionary that must contain at least these entries: 
    #   lon:          longitude grid  ([degrees_E])
    #   lat:          latitude  grid  ([degrees_N])
    #   [depth:       depth grid ]    ([m], optional)
    #   [time:        time  grid ]    ([sec since offset], optional)
    #       [time_offset: time_offset ]   (datetime like, optional)
    #       [time_unit:   unit ]          (as string, optional)
    #   variable_name: variablename to be used for data 
    #   long_name:    long title for variable_name (attribute associated with variable)
	#   units:        units for data (standard conformance not assessed; attribute associated with variable)
    # Global attributes must be in entry global_attributes - all other attributes are assumed to be variable attributes
    # writeopt must contain "data_layout" (permutation of 'tzyx', corresponding to axis in data)
    # and possibly "time_frame_number" to write specific time frames
    # variable_name is declared according to appearence of (lon,lat,depth,time) in metadata in COARDS order
    # (lon,lat,depth,time) also becomes netcdf dimensions with same name
    #
    # COARDS notes: dimensions should appear in the relative order T, then Z, then Y, then X in the CDL
    # -----------------------------------------------------------------------------------------
    # --- mandatory dimension lon ---
    if ncfile.dimensions.has_key('lon'): # lon dimension exist, check size
        assert ncfile.dimensions['lon'].size == size(metadata["lon"])
    else:
        ncfile.createDimension('lon', size(metadata["lon"]))
    # ---  mandatory dimension lat --
    if ncfile.dimensions.has_key('lat'): # lat dimension exist, check size
        assert ncfile.dimensions['lat'].size == size(metadata["lat"])
    else:
        ncfile.createDimension('lat', size(metadata["lat"]))
    # ---- optional dimensions time, depth --
    #      for time, signal unlimited dimension by None, if requested
    
    if ncfile.dimensions.has_key('depth'): # lat dimension exist, check size
        assert ncfile.dimensions['depth'].size == size(metadata["depth"])
    elif metadata.has_key("depth"):
        ncfile.createDimension('depth', size(metadata["depth"]))
        
    if metadata.has_key("time") and not ncfile.dimensions.has_key('time'):
        if metadata["time"] is None:
            ncfile.createDimension('time', None)                   # unlimited
        else:
            ncfile.createDimension('time', size(metadata["time"])) # specific size
            
    # -------- create axis variables as needed  --------
    
    if not ncfile.variables.has_key('lon'): # if exist, assume OK and do not redefine
        lon = ncfile.createVariable('lon', 'd', ('lon',))
        lon.long_name      = "Uniformly spaced longitudes"
        lon.cartesian_axis = "X" 
        lon.units          = "degrees_E" 
        lon.ipositive      = 1 
        lon[:] = metadata['lon']        
    # ----     
    if not ncfile.variables.has_key('lat'): # if exist, assume OK and do not redefine
        lat = ncfile.createVariable('lat', 'd', ('lat',))
        lat.long_name      = "Uniformly spaced latitudes"
        lat.cartesian_axis = "Y" 
        lat.units          = "degrees_N" 
        lat.ipositive      = 1 
        lat[:] = metadata['lat']
    # ----
    if metadata.has_key("depth") and not ncfile.variables.has_key('depth'):
        depth = ncfile.createVariable('depth', 'd', ('depth',))
        depth.long_name     = "depth below sea surface"
        depth.units          = "m" 
        depth[:] = metadata['depth']
    # ----
    if metadata.has_key("time") and not ncfile.variables.has_key('time'):
        time = ncfile.createVariable('time', 'i', ('time',))
        if metadata.has_key("time_offset"):
            time_offset = metadata["time_offset"].strftime('%Y-%m-%d %H:%M:%S') # datetime has strftime method
        else:
            time_offset = "undeclared offset"
        if metadata.has_key("time_unit"):
            tunit       = metadata["time_unit"]
        else:
            tunit       = "undeclared unit" # default
        time.long_name      = "time (since %s) corresponding to this time frame" % time_offset
        time.units          = "%s since %s" % (tunit, time_offset)
        # do not set, if time dimension is unlimited
        if metadata["time"] is not None:
            time[:] = metadata['time']
            
    # -------- create data variable, if needed  --------
    
    vardims = () # observe COARDS order
    for dim in ('time','depth','lat','lon'):
        if ncfile.dimensions.has_key(dim):
            vardims = vardims + (dim,)
              
    if not ncfile.variables.has_key(metadata['variable_name']): # if exist, assume OK and do not redefine
        var = ncfile.createVariable(metadata['variable_name'], data.dtype, vardims)
    else:
        var = ncfile.variables[metadata['variable_name']] # variable already exist
        
    # ---- dump data as full array or time frame ----
    #      do not check data rank, to allow a declaration cycle and e.g. setting data to zero
    reorder = reorder_axes_to_COARDS(writeopt["data_layout"])
    if writeopt.has_key("time_frame_number"):
        it = writeopt["time_frame_number"]
        var[it,:] = transpose(data, axes=reorder)
        # store time for this time frame, if it is provided as a single number
        if metadata.has_key("time") and isinstance(metadata['time'], numbers.Number):
            ncfile.variables['time'][it] = metadata['time']
    else:
        var[:] = transpose(data, axes=reorder) 
       
           
    # dump meta data 
    # check mandatory metadata - remaining mandatory keys generate other exceptions, if not present
    if not metadata.has_key('units'):
        raise exceptions.KeyError("metadata has no unit entry")
    if not metadata.has_key('long_name'):    
        raise exceptions.KeyError("metadata has no long_name entry")
    # ------ global attributes ------
    ncfile.Conventions = "COARDS"
    if metadata.has_key('global_attributes'):
        gattr = metadata['global_attributes']
        for key in gattr.keys():
            setattr(ncfile, key, gattr[key])
    # ------ variable attributes ------
    for key in metadata.keys():
        if key in ('lon', 'lat', 'depth','time','time_unit','time_offset','variable_name', 'global_attributes'): # these are handled above
            continue
        if key == 'units':
            var.units = metadata['units']
            continue
        if key == 'long_name':
            var.long_name = metadata['long_name']
            continue
        # all other metadata stored as global attributes
        setattr(var, key, metadata[key])

## Compute axis reordering of an array to conform with COARDS: T(time)-Z(depth)-Y(longitude)-X(latitude)
def reorder_axes_to_COARDS(layout):
    # @layout: axes of an array declared as any combination/selection of 'tzyx'
    # @return: permutation of layout to conform with COARDS ordering 'tzyx'(c-order)
    # --- check dublets  ---
    for a in layout:
        if (layout.count(a)>1):
            raise exceptions.ValueError("reorder_axes_to_COARDS: layout = %s contain dublets" % layout)
    # --- check invalid axis designations  ---
    for a in layout:
        if ('tzyx'.find(a)<0):
            raise exceptions.ValueError("reorder_axes_to_COARDS: layout = %s contain invalid axis designation" % layout)
    # --- somewhat clumsy
    reorder = []
    for a in 'tzyx':
        i = layout.find(a) # -1 if a does not appear in layout
        if (i>=0):
            reorder.append(i)
    return tuple(reorder)
    
        
    
    
    
      
# ========================================================================
#                     self test section
# ========================================================================
if __name__ == "__main__":
    import netCDF4    as netcdf
    from numpy import *
    ncf = netcdf.Dataset("test.nc", "w")
    lon = arange(1, 3, 0.22)
    lat = arange(10, 13, 0.52)
    nx,ny = (size(lon),size(lat))
    var = random.random((nx,ny))
    meta = {'lon':lon,
            'lat':lat,
            'variable_name':'georg',
            'long_name': "GyroGearloose",
            'units':"Gy",
            'global_attributes': {'history': 'xxx'}}
    write_lonlatdata_in_COARDS_format(ncf, var, meta)
    # --- dump second variable in same set
    meta['variable_name'] = 'mask'
    meta['long_name']     = 'maskmask'
    meta['global_attributes'] = {'morehistory': 'xxx'} 
    write_lonlatdata_in_COARDS_format(ncf, ones((nx,ny), int), meta)
    ncf.close()
    #
