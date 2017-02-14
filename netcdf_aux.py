#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################################################################
##       @package netcdf_aux
#        Auxillary functions for convenient netcdf I/O
#        
#######################################################################################################################
import os

_thisdir         = os.path.dirname(__file__)  # allow remote import
_verbose         = False # True  # log info

## Append data + index along unlimited dimension after last stored frame
#  currently index-name is hardcoded to "index"
#  @param ncfile        NetCDFFile instance
#  @param data          2D data to be written
#  @param index         index value for this data frame
#  @param dataParam     Optional keyword argument:  Custom name for data variable. Default: data
#  Assumes netcdf variable and index variable is already created
#
def add_data_to_unlimted_var(ncfile, data, **kwargs):
    if kwargs.has_key("dataParam"):  # Default: data
        dname = kwargs["dataParam"]  # apply custom data variable name
    else:
        dname = "data"               # default data variable name
    index                                = kwargs["index"] # mandatory keyword argument
    nextpos                              = ncfile.variables[dname].shape[0]
    ncfile.variables[dname][nextpos,:,:] = data
    ncfile.variables["index"][nextpos]   = index


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
    
