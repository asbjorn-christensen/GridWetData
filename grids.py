#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################################################################
##       @package grid
#        Light weight grid classes, nucleated from the IBMlib module mesh_grid representing most
#        common grid types.
#        Offer interpolations of data associated with grid. 3D grids may represent sea bed topography
#        Subclasses of this module focuses on the HBM data grids, super classes are generic.
#        Grids are not aware of data that may be associated with them (apart from topography)
#
#        <b> general considerations </b>.
#           grids are assumed clonable by a shallow copy. Dynamic attributes (e.g. sea level, affecting water depth) 
#           of a grid must be updated by assigning new attributes (not inplace modifications), so ancestor grid is not
#           affected by update of child grid (or vice versa).
#           In doc strings, (x,y,z) refer to grid coordinates, (lambda,phi,depth) Cartesian coordinates, verically referenced
#           to dynamic sea surface
#
#        <b> General interpolation principles </b>.
#           \li 1) apply trilinear interpolation whenever possible (interpolation point surrounded by vertices)
#           \li 2) do not extrapolate vertically above first cell center / below bottom / beyond horizontal grid rim cell centers,
#                  but collapse vertices in these boundary cases
#           \li 3) apply only wet point vertices in interpolation
#
#######################################################################################################################
import sys
import os
import exceptions
#from   scipy.io import netcdf   # there is a bug in scipy.io.netcdf < 0.12 in writing part, 0.17.0 has a writing issue 
import netCDF4 as netcdf         # use this API until scipy.io.netcdf is fixed                            

from   numpy import *
import copy                     # for grid copy, notice numpy also provides a copy which is shadowed by this import
from   datetime import *
import re


_thisdir         = os.path.dirname(__file__)  # allow remote import
_verbose         = False # True  # log info

# ===============================================================================
#             Error structure of grid world
# ===============================================================================

## Super class for grid-related error conditions
class GridRangeViolation(exceptions.ValueError): pass

## Raise when horizontal range of grid is exceeded
class HorizontalRangeViolation(GridRangeViolation): pass

## Raise when horizontal range of grid is exceeded
class TopographicalViolation(GridRangeViolation): pass

## Raise when a wet point issue is assessed in a land point (no water at all)
class IsLandPoint(TopographicalViolation) : pass

## Raise when a wet point issue is assessed below the sea bed (but there is water above)
class BelowSeaBed(TopographicalViolation) : pass

## Raise when a wet point issue is assessed above the sea surface (but there is water below)
class AboveSeaSurface(TopographicalViolation): pass


# ===============================================================================
## Horizontal grid regulary spaced in lon-lat super class
#
# Not aware of wet/dry conditions
# grid coordinates are vertex coordinates with vertices at integers >= 0
# (x,y) = (0,0) corresponds to (lon0, lat0)
# Coordinate transformations:
#    (lam,phi) = lon0 + x*dlon,    lat0 + y*dlat
#    ( x , y ) = (lam-lon0)/dlon, (phi-lat0)/dlat
# Notice cell indices != vertex indices - cells are centered around (x,y) = integers
# Currently metric is not implemented, but when this is added, units are degrees E/N
# 
# TODO: Add optional wet/dry awareness, defaulting to all-wet
# ===============================================================================
class LonLatGrid:
    ## -----------------------------------------------------
    ## constructor
    #  @param self   The object pointer.
    #  @param nx     Zonal grid points
    #  @param ny     Meridional grid points
    #  @param lon0   Zonal zero, corresponding to x=0
    #  @param lat0   Meridional zero, corresponding to y=0
    #  @param dlon   Zonal grid spacing
    #  @param dlat   Meridional grid spacing
    #  
    def __init__(self, nx, ny, lon0, lat0, dlon, dlat):
        ## zonal (x) grid dimension
        self.nx   = nx
        ## meridional (y) grid dimension
        self.ny   = ny
        assert nx > 1
        assert ny > 1
        ## Zonal zero, corresponding to x=0
        self.lon0 = lon0
        ## Meridional zero, corresponding to y=0
        self.lat0 = lat0
        ## Zonal grid spacing
        self.dlon = dlon
        ## Meridional grid spacing
        self.dlat = dlat
        
    #  ---------------------------------------------------------------------------------
    ## unrestricted horizontal coordinate transformation from grid coordinates to lon,lat
    #  @param self   The object pointer
    #  @param xy     any sequence (x,y) of grid coordinates
    #  @return shpere cooordinate (lon,lat)
    #
    def get_lonlat(self, xy):
        return self.lon0 + xy[0]*self.dlon,   self.lat0 + xy[1]*self.dlat

    #  ---------------------------------------------------------------------------------
    ## unrestricted horizontal coordinate transformation from lon,lat to grid coordinates
    #  @param self   The object pointer
    #  @param lonlat any sequence (lon,lat) 
    #  @return grid cooordinate (x,y)
    #
    def get_xy(self, lonlat):
        return (lonlat[0]-self.lon0)/self.dlon, (lonlat[1]-self.lat0)/self.dlat

    #  ---------------------------------------------------------------------------------
    ## Resolve vertex coordinates: SW vertex (ix,iy) + NE vertex (ixp1,iyp1) + intra vertex displacement (sx,sy)
    #  get_vertex_indices also define the interpolation choice at the grid rim: confine
    #       0 <= (ix,ixp1) <= nx-1
    #       0 <= (iy,iyp1) <= ny-1
    #  so that (ix,iy, ixp1,iyp1) are referencable vertices.
    #  Vertices collapse at grid rim:
    #      (-0.5 < x < 0 or nx-1 < x < nx-0.5 )
    #      (-0.5 < y < 0 or ny-1 < y < ny-0.5 )
    #  but HorizontalRangeViolation is raised beyond grid rim
    #
    #  @param self The object pointer
    #  @param pos  any sequence (lon,lat <, ...>)
    #  @return vertices + intra vertex displacement (ix,iy, ixp1,iyp1, sx,sy)
    #
    def get_vertex_indices(self, pos):
        x,y = self.get_xy(pos)
        if x < -0.5 or x > self.nx-0.5:  # grid rim
            raise HorizontalRangeViolation("get_vertex_indices: x=%f" % x)
        if y < -0.5 or y > self.ny-0.5:  # grid rim
            raise HorizontalRangeViolation("get_vertex_indices: y=%f" % y)
        ix    = int(x)  # SW corner, int(-1 < z < 1) = 0
        iy    = int(y)  # SW corner, int(-1 < z < 1) = 0
        ixp1  = max(0, min(self.nx-1, ix+1))  # NE corner, possibly collapsed
        iyp1  = max(0, min(self.ny-1, iy+1))  # NE corner, possibly collapsed
        sx    = x-ix   # also for collapsed vertices
        sy    = y-iy   # also for collapsed vertices
        return (ix,iy, ixp1,iyp1, sx,sy)

    #  ---------------------------------------------------------------------------------
    ## Resolve which cell pos belongs to; cells are centered aorund vertices (x,y) = (ix,iy)
    #  confine 0 <= (ix,iy) <= (nx-1,ny-1) so that (ix,iy) is a referencable vertex (project in pos)
    #  Notice cell indices != vertex indices. 
    #  @param self The object pointer
    #  @param pos  any sequence (lon,lat <, ...>)
    #  @return horizontal cell indices (ixc,iyc) of cell surrounding pos
    # 
    def get_cell_indices(self, pos):
        x,y  = self.get_xy(pos)
        #
        ixc  = int(x+0.5)
        if ixc < 0 or ixc > self.nx-1:
            raise HorizontalRangeViolation("get_cell_indices: x=%f" % x)
        #
        iyc  = int(y+0.5)
        if iyc < 0 or iyc > self.ny-1:
            raise HorizontalRangeViolation("get_cell_indices: y=%f" % y)
        #
        return (ixc,iyc)

    #  ---------------------------------------------------------------------------------
    ## Flag whether pos is inside overall horizontal grid (not considering whether point is wet/dry)
    #  @param self The object pointer
    #  @param pos  any sequence (lon,lat <, ...>)
    def is_inside_grid(self, pos):
        try:
            ixc,iyc = self.get_cell_indices(pos)
            return True  # valid cell indices could be resolved
        except HorizontalRangeViolation:
            return False 
        
    
    # -----------------------------------------------------------------------
    ## Trilinear unconditional interpolation between corners of field arr2D sampling f(lam,phi) at pos = (lam,phi)
    # deriv=0 -> interpolation of f at pos = (lam,phi)
    #      deriv=1 -> spatial gradioent of f at pos = (lam,phi)
    #                  df = [df_dl, df_dl]
    #      where dl,dp are lamdba, phi derivatives, respectively
    # numerical deriv test OK (when pos not on vertices/vertex lines, where arr2D is not differentiable)
    # Interpolation is unconditional since LonLatGrid is not aware of wet/dry conditions
    # @param self The object pointer
    # @param arr2D[nx,ny] 2D array, where pos is interpolated
    # @param deriv whether or not to take derivative wrt. (lam,phi) og result
    #         deriv=0 -> interpolation of f at pos = (lam,phi)
    #         deriv=1 -> spatial gradioent of f at pos = (lam,phi) df = [df_dl, df_dl]
    # @return data interpolated at pos
    # 
    def interpolate_data(self, arr2D, pos, deriv=0):
        ix,iy,ixp1,iyp1,sx,sy = self.get_vertex_indices(pos)
        res = self.interpolate_data_from_grid_coordinates(arr2D, ix,iy,ixp1,iyp1,sx,sy, deriv)
        if deriv == 1:  # transform from (sx,sy) derivative to (lamdba,phi) derivative
            res /= array([self.dlon, self.dlat])
        return res


    # -----------------------------------------------------------------------
    ## Trilinear unconditional interpolation between corners of field arr2D sampling f(lam,phi)
    #  Interpolation is unconditional since LonLatGrid is not aware of wet/dry conditions
    #  It is assumed that (ix,iy,ixp1,iyp1) corresponds to referencable vertices in arr2D (not checked)
    #  At the interior of the grid, ixp1 = ix+1 and iyp1 = iy+1
    # @param self The object pointer
    # @param arr2D[nx,ny] 2D array, where pos is interpolated
    # @param ix,iy,sx,sy vertex coordinates for interpolation
    # @param deriv whether or not to take derivative wrt. (lam,phi) og result
    #      deriv=0 -> interpolation of f at grid coordinates (ix,iy,sx,sy)
    #      deriv=1 -> spatial gradient wrt grid coordinates of f at (ix,iy,ixp1,iyp1,sx,sy) df = [df_dsx, df_dsy]
    #      where dsx,dsy are grid coordinate derivatives, respectively
    # @return data interpolated at (ix,iy,ixp1,iyp1,sx,sy)
    # 
    def interpolate_data_from_grid_coordinates(self, arr2D, ix,iy, ixp1,iyp1, sx,sy, deriv=0):
        if deriv == 0:
            f =  arr2D[ix,  iy]   * (1-sx) * (1-sy)   
            f += arr2D[ixp1,iy]   *   sx   * (1-sy)
            f += arr2D[ix,  iyp1] * (1-sx) *   sy
            f += arr2D[ixp1,iyp1] *   sx   *   sy
            return f
        elif deriv == 1: 
            # --- sx derivative 
            df_dsx  = arr2D[ix,  iy]   * (-(1-sy))  
            df_dsx += arr2D[ixp1,iy]   * (1-sy)
            df_dsx += arr2D[ix,  iyp1] * (-sy)
            df_dsx += arr2D[ixp1,iyp1] *  sy 
            # --- sy derivative 
            df_dsy  = arr2D[ix,  iy]   * (-(1-sx))   
            df_dsy += arr2D[ixp1,iy]   * (-sx) 
            df_dsy += arr2D[ix,  iyp1] * (1-sx) 
            df_dsy += arr2D[ixp1,iyp1] *  sx   
            return array([df_dsx, df_dsy], float)
        else:
            print "interpolate_2Ddata: unknown deriv argument value", deriv
            raise ValueError

    ## provide x,y mesh for the grid:
    # matplotlib / matlab compartability function
    # @param self The object pointer
    # @return (x[nx,ny], y[nx,ny]) of grid mesh
    #
    def mgrid(self):
        x0 = self.lon0
        x1 = self.lon0 + (self.nx-1)*self.dlon + 1.0e-12 # include this point
        y0 = self.lat0
        y1 = self.lat0 + (self.ny-1)*self.dlat + 1.0e-12 # include this point
        return mgrid[x0:x1:self.dlon, y0:y1:self.dlat]
    
    # ==================================================== writers ===================================================
    #
    ## Write/append data in netCDF format to file fname
    #  @param self          The object pointer.
    #  @param file          File name / NetCDFFile instance for writing data
    #  @param data          File name provided: write plain 2D array 
    #                       NetCDFFile instance provided: append data to variable data along unlimited dimension
    #  @param undef         Optional keyword argument:  value of data corresponding to undefined grid points (just passed to netcdf file). Default: None
    #  @param bitmask       Optional keyword argument:  1 for valid points, 0 for undefined grid points. Default: None
    #  @param specification Optional keyword argument:  specification title attribute associated with data. Default: None
    #  @param metadata      Optional keyword argument:  a dictionary passed to global attributes. Default: autotitle
    #  @param storage_type  Optional keyword argument:  netcdf API supports these types:          Default: select type corresponding to data <br> 
    #                                                   Supported values: <br> 
    #                                                   'b' : NC_BYTE:   (1 byte)  <br> 
    #                                                   'c' : NC_CHAR:   (1 byte)  <br> 
    #                                                   'h' : NC_SHORT:  (2 bytes) <br> 
    #                                                   'i' : NC_INT:    (4 bytes) <br> 
    #                                                   'f' : NC_FLOAT:  (4 bytes) <br> 
    #                                                   'd' : NC_DOUBLE: (8 bytes) <br>
    #  @param index         Optional/mandatory keyword argument: value for index (mandatory for append, ignored for write)
    #
    #  undef/bitmask are overlapping ways of flagging undefined grid points. If both are defined both are written, do not assess consistency
    #  In append mode, optional meta data is not updated, but pertains to all frames in data
    #  Currently it is hardcoded that the data is stored in netCDF variable "data(nx,ny)" in single frame mode
    #  and "data(nframes,nx,ny)" indexed by "index(fnframes)" in append mode
    # 
    def write_data_as_netCDF(self, file, data, **kwargs):
        if isinstance(file, basestring):   # write single frame
            ncfile  = netcdf.Dataset(file, "w")
            kwargs["index"] = None  # suppress append mode
            self._setup_netCDF_data_set(ncfile, data, **kwargs)
            ncfile.close()
        elif isinstance(file, netcdf.Dataset):
            if file.dimensions == {}:  # empty set, test also OK for python-netcdf
                self._setup_netCDF_data_set(file, data, **kwargs) # index mandatory
            else: # assume properly configured
                self._add_netCDF_data_(file, data, kwargs["index"])  # index mandatory
        else:
            
            raise exceptions.ValueError("argument file inappropriate: file = %s" % str(file))

    # internal method: append data + index along unlimited dimension after last stored frame
    def _add_netCDF_data_(self, ncfile, data, index):
        inext                                = ncfile.variables["data"].shape[0]
        ncfile.variables["data"][inext,:,:]  = data
        ncfile.variables["index"][inext]     = index

          
    # internal method: configure netcdf data set and write data
    # If index is None (default): single frame write
    # If index is number:         store data in append mode, first index=index
    def _setup_netCDF_data_set(self, ncfile, data, **kwargs):
        assert data.shape == (self.nx, self.ny)
        #
        #  --------- parse kwargs ---------
        #
        if kwargs.has_key("undef"):  # Default: None
            ncfile.undef         = kwargs["undef"]
            ncfile.undef_meaning = "value for undefined data points"
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
        if kwargs.has_key("bitmask"):
            bitmask = kwargs["bitmask"]
        else:
            bitmask = None # Default: no bitmask saved
        # metadata saved as global attributes in netCDF set
        if kwargs.has_key("metadata"):
            metadata = kwargs["metadata"]
        else:
            metadata = {"autotitle": "data on regular lon-lat grid exported from LonLatGrid.write_data_as_netCDF"}
        for key in metadata:
            setattr(ncfile, key, metadata[key])
        #
        #  --------- create dimensions ---------
        #
        index = kwargs["index"] # mandatory keyword argument
        if index is not None:   
            ncfile.createDimension('nframes', None)   # None stands for the unlimited dimension
            ncfile.nframes = "number of frames"
        ncfile.createDimension('nx', self.nx)
        ncfile.nx = "Zonal (west-east) grid dimension"
        ncfile.createDimension('ny', self.ny)
        ncfile.ny = "Meridional (south-north) grid dimension"
        #
        #  --------- create netcdf variables ---------
        #
        if index is not None: # append mode
            var_index = ncfile.createVariable("index", 'd',       ('nframes',))
            var_index.specification = "frame index"
            var_data  = ncfile.createVariable("data",  data_type, ('nframes','nx','ny')) 
        else:   #  single frame mode
            var_data = ncfile.createVariable("data", data_type, ('nx','ny'))
        #
        if data_specification:
                var_data.specification = data_specification
        #
        var_lon0 = ncfile.createVariable('lon0', 'd', ())
        var_lon0.specification = "Western grid edge in degrees"
        #
        var_lat0 = ncfile.createVariable('lat0', 'd', ())
        var_lat0.specification = "Southern grid edge in degrees"
        #
        var_dlon = ncfile.createVariable('dlon', 'd', ())
        var_dlon.specification = "Zonal grid increment in degrees, positive west-to-east"
        #
        var_dlat = ncfile.createVariable('dlat', 'd', ())
        var_dlat.specification = "Meridional grid increment in degrees, positive south-to-north"
        #
        if bitmask:
            var_bitmask = ncfile.createVariable('bitmask', 'b', ('nx','ny')) # save as byte
            var_bitmask.specification = "Bitmask for data: 1 = valid, 0 = invalid)"
        #
        #  --------- assign netcdf variables ---------
        #
        var_lon0.assignValue(self.lon0)
        var_lat0.assignValue(self.lat0)
        var_dlon.assignValue(self.dlon)
        var_dlat.assignValue(self.dlat)
        if bitmask:
            var_bitmask[:,:] = bitmask # don't perform shape/value check
        #
        if index is not None: # write data as first frame
            var_index[0]    = index
            var_data[0,:,:] = data
        else:                 # write data as single frame 
            var_data[:,:] = data
        # 
        
        
    
    
# -------------------------------------------------------------------------------------------
## Horizontal lon-lat grid super class with flexible vertical structure, where cells have reference + dynamic thickness
# currently dynamic sea level fluctuations are attributes to surface layer
# cellw0[ix,iy,iz] at cell center line (x,y) = (ix,iy) along with a
# dynamic sea surface elevation z[ix,iy] associated with cellw0[ix,iy,0]
# Based on IBMlib module mesh_grid
#
# <b> Topographic conventions: </b>
#     Grid coordinates are integer values >= 0 at vertices. Vertical cell separations are strictly vertical, so
#     same horizontal grid applies to each horizontal layer.
#     wetmask[ix,iy,iz] tells whether the volume (x,y,z) with grid coordinates [ix+/-0.5, iy +/-0.5, iz +/-0.5] is wet
#     actual water depth at xy is obtained by unconditional trilinear interpolation of layer_sep[:,:, bottom_layer[ix,iy]]
#     where (ix,iy) are the indices of the cell where to xy belongs
#     
#     This implies coast line kinks ("cliffs") in general. Cliffs can be minimized is cellw are manipulated for land points (not done)
#     Sea surface has grid coordinate -0.5, sea bed by z = layer_sep(x,y,bottum_layer(x,y)) 
#
#  TODO: check consistency between grid coordinates and wdepth interpolation (with/without kinks at sea bed
#        due to bottom_layer), and wet/dry flagging
class LonLatZGrid(LonLatGrid):
    #  -----------------------------------------------------
    ## constructor
    #  @param self   The object pointer.
    #  @param nx     Zonal grid points
    #  @param ny     Meridional grid points
    #  @param nz     Number of vertical points (some may be dry)
    #  @param lon0   Zonal zero, corresponding to x=0
    #  @param lat0   Meridional zero, corresponding to y=0
    #  @param dlon   Zonal grid spacing
    #  @param dlat   Meridional grid spacing
    #  @param cellw0[nx,ny,nz]  reference thickness of cells
    #
    #  grid is initialized corresponding to sea level elevation = 0
    #
    def __init__(self, nx, ny, lon0, lat0, dlon, dlat, cellw0):
        LonLatGrid.__init__(self, nx, ny, lon0, lat0, dlon, dlat)
        ## number of vertical grid points (some may be dry)
        self.nz    = cellw0.shape[2]
        self.configure_grid(cellw0)
        
    def export_horizontal_grid(self):
        return LonLatGrid(self.nx, self.ny, self.lon0, self.lat0, self.dlon, self.dlat)
    
    # -------------------------------------------------
    ## Set static auxillaries attributes from reference cell thickness cellw0_in
    #  Dry cells are signalled by cellw0_in <= 0
    #  @param self                 The object pointer.
    #  @param cellw0_in[nx,ny,nz]  reference thickness of cells
    #
    def configure_grid(self, cellw0_in):
        nx,ny,nz = self.nx, self.ny, self.nz
        assert cellw0_in.shape == (nx,ny,nz)
        ## wetmask[ix,iy,iz] = 1/0 if cell is wet/dry
        self.wetmask = where(cellw0_in>0, 1, 0)
        ## cellw0[ix,iy,iz] provided reference cell thickness at cell center line (x,y) = (ix,iy) for cell iz
        self.cellw0  = where(cellw0_in>0, cellw0_in, 0.0) # set cellw0 = 0 for dry cells here after      
        #
        layer_sep0 = 1.0*self.cellw0 # force copy
        for iz in range(1,nz):
            layer_sep0[:,:,iz] = layer_sep0[:,:,iz-1] + self.cellw0[:,:,iz]
        ## layer_sep[ix,iy,iz] is the dynamic depth of the boundary if layers (iz, iz+1) at cell center (ix,iy)
        self.layer_sep0 = layer_sep0
        #
        bottom_layer = sum( self.wetmask, axis=2) - 1 # assume wet layers contiguous from surface and down
        ## bottom_layer[ix,iy] is the lowest wet layer, when bottom_layer[ix,iy] >= 0
        #  bottom_layer[ix,iy] < 0 for dry cells. Layer 0 is surface layer.
        self.bottom_layer = bottom_layer
        # --- derive cell center depth ---
        ccdepth0      = zeros((nx,ny,nz), float)
        for ix in range(self.nx):
            for iy in range(self.ny):
                ibot = bottom_layer[ix,iy]
                if ibot >= 0: # layer 0 is surface layer
                    ccdepth0[ix,iy,0] = 0.5*self.cellw0[ix,iy,0]  # bottom_layer[ix,iy] >= 1
                    for iz in range(1, 1+ibot):
                        thickness          = 0.5*(self.cellw0[ix,iy,iz-1] + self.cellw0[ix,iy,iz])
                        ccdepth0[ix,iy,iz] = ccdepth0[ix,iy,iz-1] + thickness
                    if ibot < nz-1:
                        ccdepth0[ix,iy,(1+ibot):] = ccdepth0[ix,iy,ibot] # terminate column corresponding to cellw0=0
        ## ccdepth0[ix,iy,iz] is reference cell center depth at the cell center line (x,y) = (ix,iy)
        self.ccdepth0     = ccdepth0
        # 
        self.set_reference_level( zeros((nx,ny), float))


    # -----------------------------------------------------------
    ## Update dynamic grid attributes (cellw, ccdepth, layer_sep) corresponding
    # corresponding to sea level elevation z_in[nx,ny] (positive up)
    # Method allows grid cloning, since attributes (cellw, ccdepth, layer_sep) are not
    # updated in-place
    # @param self        The object pointer.
    # @param z_in[nx,ny] sea level elevation with respect to reference level (positive up)
    #
    def set_reference_level(self, z_in):
        nx,ny,nz = self.nx, self.ny, self.nz
        ## z[ix,iy] is the sea surface elevation from reference level (positive up) of surface cell
        self.z             = where(self.wetmask[:,:,0]>0, z_in, 0.0) # make sure z is zero in dry points
        ## cellw[ix,iy,iz] is dynamic cell thickness at cell center line (x,y) = (ix,iy) for cell iz
        self.cellw         = 1.0*self.cellw0 # force copy - does not overwrite existing array self.ccdepth
        self.cellw[:,:,0] += self.z
        #
        zcontrib       = ones(nz, float)
        fluct          = reshape(outer(self.z.flat, zcontrib), (nx,ny,nz))
        ## layer_sep[ix,iy,iz] is the dynamic depth of the boundary if layers (iz, iz+1) at cell center (ix,iy)
        self.layer_sep = self.layer_sep0 + fluct # does not overwrite existing array self.layer_sep
        #
        fluct[:,:,0] *= 0.5 # z only counts half in first layer for ccdepth
        ## ccdepth[ix,iy,iz] is dynamic cell depth at cell center line (x,y) = (ix,iy) for cell iz
        self.ccdepth  = self.ccdepth0 + fluct # does not overwrite existing array self.ccdepth

    # -----------------------------------------------------------
    ## export sea level elevation z (no copy)    
    def get_reference_level(self):
        return self.z

    # ---------------------------------------------------------------------
    ## Resolve grid indices and weights needed for interpolating grid data
    #  get_vertex_indices also define the interpolation choice at the grid rim, surface + bottom
    # Only resolve vertical components for wet cells, because
    # cellw = 0 for dry cells. Return surface position indices for dry cells
    # When vertex indices are resolved (i.e. no exception raised), the
    # returned vertx (ix,iy,iz) is guarantied to be a wet point so an
    # interpolation using returned vertex indices is non-singular. 
    # Special cases: set sz = 0, when we are above center of first layer or below center
    # of bottom layer to avoid unwarrented vertical extrapolaitons
    # @param self    The object pointer.
    # @param pos     any sequence (lon,lat <, ...>)
    # @return vertex + intra vertex displacement (ix,iy,iz, ixp1,iyp1,izp1, sx,sy,sz)
    #
    def get_vertex_indices(self, pos):
        # --- resolve horizontal components ---
        ix,iy,ixp1,iyp1,sx,sy = LonLatGrid.get_vertex_indices(self, pos) # horizontal grid coordinates
        ixc,iyc               = LonLatGrid.get_cell_indices(self, pos)  # cell indices
        if self.wetmask[ixc,iyc,0] == 0:
            raise IsLandPoint(str(pos)) 
        # --- resolve vertical components ---
        # 1) interpolate cell centers horizontally to pos -> ccd, moving down the water column
        depth = pos[2]
        if depth<0:
            raise AboveSeaSurface(str(pos))
        #
        # --- loop for iz,sz ---
        #     scan down layer wise to see, which vertices apply to pos vertically
        #     to identify (iz,izp1,sz). 
        #
        iz = izp1 = 0
        ccd_last  = -1e-3 # above sea surface
        while iz <= self.bottom_layer[ixc,iyc]:
            # ccd is the layer center depth of layer iz, interpolated to horizontal position pos[:2]
            ccd = LonLatGrid.interpolate_data_from_grid_coordinates(self, self.ccdepth[:,:,iz], ix,iy,ixp1,iyp1,sx,sy)
            if ccd > depth: 
                iz = max(0, iz-1)
                sz = (depth - ccd_last)/(ccd - ccd_last) # ccd_last becomes defined at iz == 0
                break  # while loop
            else: # prepare next iteration on iz
                ccd_last = ccd
                iz      += 1
                izp1    += 1
        else: # loop exhausted at bottom_layer, project to valid range
            iz = izp1  = self.bottom_layer[ixc,iyc]
            virtual_seabed = self.layer_sep[:, :, iz] 
            cwd = LonLatGrid.interpolate_data_from_grid_coordinates(self, virtual_seabed, ix,iy,ixp1,iyp1,sx,sy)
            if depth>cwd:
                raise BelowSeaBed(str(pos))
            else:
                sz = (depth-ccd)/cwd # assign a well-def value
        return ix,iy,iz,ixp1, iyp1,izp1, sx,sy,sz
           
    # ---------------------------------------------------------------------------------
    ## Resolve which cell pos belongs to; cells are centered around vertices (x,y,z) = (ix,iy,iz)
    #  confine 0 <= (ix,iy,iz) <= (nx-1,ny-1,nz) so that (ix,iy,iz) is a referencable vertex (project in pos)
    #  Raise IsLandPoint, if point is dry, AboveSeaSurface/BelowSeaBed if pos is not in the water column
    #  @param self The object pointer
    #  @param pos  any sequence (lon,lat, depth)
    #  @return cell indices (ixc,iyc, izc) of cell surrounding pos
    #
    def get_cell_indices(self, pos):
        # --- resolve horizontal components 
        ix,iy,ixp1,iyp1,sx,sy = LonLatGrid.get_vertex_indices(self, pos) # horizontal grid coordinates
        ixc,iyc               = LonLatGrid.get_cell_indices(self, pos)  # cell indices
        if self.wetmask[ixc,iyc,0] == 0:
            raise IsLandPoint(str(pos))
        # --- resolve izc ---
        # 1) interpolate cell layer separations horizontally to pos -> dsep, moving down the water column
        depth = pos[2]
        if depth<0:
            raise AboveSeaSurface(str(pos))
        # --- loop for iz,sz ---
        izc    = 0
        while izc < self.nz:
            dsep = LonLatGrid.interpolate_data_from_grid_coordinates(self, self.layer_sep[:,:,iz], ix,iy,ixp1,iyp1,sx,sy)
            if dsep > depth: #  > 0
                break 
            izc      += 1     
        else: # izc = nz, project to valid range
            raise BelowSeaBed(str(pos))
        #
        return ixc,iyc,izc

    # ------------------------------------------------------------------------------
    ## The water depth in any wet position is defined by unconditional trilinear interpolation of layer_sep[:,:, bottom_layer[ix,iy]]
    # Here (ix,iy) are the indices of the cell where to pos belongs
    # Seabed kinks are defined by boundaries of centered cells, if number of
    # wet layers are varying
    # 
    # @param self The object pointer
    # @param pos  any sequence (lon,lat <, ...>)
    # @return water depth at pos, None if dry point
    #
    def interpolate_wdepth(self, pos):
        ix,iy,ixp1,iyp1,sx,sy = LonLatGrid.get_vertex_indices(self, pos)  # horizontal grid coordinates
        ixc,iyc               = LonLatGrid.get_cell_indices(self, pos)    # cell indices
        if self.wetmask[ixc,iyc,0] == 0:
            raise IsLandPoint(str(pos))
        else:
            # virtual_seabed is the continuation of the sea bed in the cell corresponding to pos
            virtual_seabed = self.layer_sep[:, :, self.bottom_layer[ixc,iyc]] 
            return LonLatGrid.interpolate_data_from_grid_coordinates(self, virtual_seabed, ix,iy, ixp1,iyp1, sx,sy)
        
    # ------------------------------------------------------------------------------
    ## Test whether pos is a wet point
    #  HorizontalRangeViolation exceptions are passed to calling contexts
    #  @param self The object pointer
    #  @param pos  any sequence (lon,lat,depth)
    #  @return if pos is inside a we cell
    def is_wet(self, pos):
        try:
            cwdpt = self.interpolate_wdepth(pos) # may raise IsLandPoint
            return 0 <= pos[2] <= cwdpt
        except IsLandPoint:
            return False
        

    # ------------------------------------------------------------------------------
    ## Interpolate grid data arr3D at position pos
    #  @param self            The object pointer
    #  @param arr3D[nx,ny,nz] Grid data to be used for interpolation
    #  @param pos             any sequence (lon,lat,depth)
    #  @return arr3D interpolated at pos, None is pos is in a dry cell
    #
    def interpolate_data(self, arr3D, pos): # currently drop deriv option
        if not self.is_wet(pos):
            return None
        else:
            ix,iy,iz,ixp1,iyp1,izp1,sx,sy,sz = self.get_vertex_indices(pos)  # horizontal grid coordinates
            # remember to divide with metric tensor when adding deriv
            return self.interpolate_data_from_grid_coordinates(arr3D, ix,iy,iz, ixp1,iyp1,izp1, sx,sy,sz)


    # ------------------------------------------------------------------------------
    ## As interpolate_data, but with grid vertex coordinates as input
    #  Perform trilinear interpolation of arr3D at for indices (ix,iy,iz, ixp1,iyp1,izp1, sx,sy,sz)
    #  where weighting is masked by wetmask so only wet points contribute to
    #  interpolation
    #  @param self            The object pointer
    #  @param arr3D[nx,ny,nz] Grid data to be used for interpolation
    #  @param ix,iy,iz        SWU corner of interpolation cube. ix,iy,iz must be within grid range (referencable), but not necessary wet
    #  @param sx,sy,sz        intra cube displacement coordinates. Interior ponts satisfy 0 < sx,sy,sz < 1
    #  @param deriv           0: value interpolation 1: derivative with respect to (sx,sy,sz)
    #  @return                arr3D interpolated at pos
    #
    def interpolate_data_from_grid_coordinates(self, arr3D, ix,iy,iz, ixp1,iyp1,izp1, sx,sy,sz, deriv=0):
        assert arr3D.shape == (self.nx, self.ny, self.nz)
        if deriv == 0:
            wsum = 0
            f    = 0
            for (jx,wx) in ((ix, 1-sx), (ixp1, sx)):
                for (jy,wy) in ((iy, 1-sy), (iyp1, sy)):
                   for (jz,wz) in ((iz, 1-sz), (izp1, sz)):
                       w     = self.wetmask[jx,jy,jz] * wx*wy*wz
                       wsum += w
                       f    += arr3D[jx,jy,jz]*w
            if wsum < 1e-20:
                bad_point = "ix,iy,iz = %d %d %d  sx,sy,sz = %f %f %f" % (ix,iy,iz, sx,sy,sz)
                raise exceptions.ValueError("rank deficit interpolation at " + bad_point)
            else:
                return f/wsum
        elif deriv == 1: 
            # --- sx derivative 
            # df_dsx  = ..
            # df_dsy  = ..
            # df_dsz  = ..
            raise exceptions.ValueError("deriv == 1 not implemented")
        else:
            raise exceptions.ValueError("deriv == " + str(deriv) + " not implemented")

    #  ============== Generic projectors ==============
    #  -------------------------------------------------------
    ## Project surface layer of data at native grid resolution (both wet/dry points)
    #  @param  self The object pointer.
    #  @param  arr3D is the 3D array to be projected
    #  @param  padval is used for dry points
    #  @return LonLatGrid, array(nx,ny)
    def get_surface_layer(self, arr3D, padval = 0.0):
        data2D = where(self.wetmask[:,:,0]>0, arr3D[:,:,0], padval)
        grid2D = self.export_horizontal_grid()
        return grid2D, data2D
    
    #  -------------------------------------------------------
    ## Project bottom layer of data at native grid resolution (both wet/dry points)
    #  projection may be speeded up avoiding explicit loops later (e.g. a projection mask + axis sum)
    #  @param  self The object pointer.
    #  @param  arr3D is the 3D array to be projected
    #  @param  padval is used for dry points
    #  @return LonLatGrid, array(nx,ny)
    #
    def get_bottom_layer(self, arr3D, padval = 0.0):
        data2D = padval*ones((self.nx, self.ny), float)
        # transfer data at wet points
        for ix in range(self.nx):
            for iy in range(self.ny):
                if self.wetmask[ix,iy,0]>0: data2D[ix,iy] = arr3D[ix, iy, self.bottom_layer[ix,iy]]
        grid2D = self.export_horizontal_grid()
        return grid2D, data2D
    
    #  -------------------------------------------------------
    ## Generate vertical average of data at native grid resolution (both wet/dry points)
    #  projection may be speeded up avoiding explicit loops later (e.g. a dynamic projection mask + axis sum)
    #  Apply dynamic cell thickness cellw as weighting factor for avarage
    #  @param  self The object pointer.
    #  @param  arr3D is the 3D array to be projected
    #  @param  padval is used for dry points
    #  @return LonLatGrid, array(nx,ny)
    #
    def get_vertical_average(self, arr3D, padval = 0.0):
        data2D = padval*ones((self.nx, self.ny), float)
        # transfer data at wet points
        for ix in range(self.nx):
            for iy in range(self.ny):
                if self.wetmask[ix,iy,0]>0:
                    w  = self.cellw[ix, iy, :] # weight for this column
                    data2D[ix,iy] = sum(arr3D[ix, iy, :] * w)/sum(w)
        grid2D = self.export_horizontal_grid()
        return grid2D, data2D

    
    def write_data_as_netCDF(self, file, data, **kwargs):
        raise exceptions.NotImplementedError("to be implemented for class LonLatZGrid")



#################################################################################################
#    Specific grids which are data aware
#
#    The specific grids must provide a these methods:
#        load_data_frame:  loads data from a provided file name corresponding to the grid
#
#    Additionally sea level aware grids must provide a these methods:
#        update_sea_level: loads sea level from provided file name and invoke implied updates
#                          in super class (e.g. set_reference_level for LonLatZGrid)
#################################################################################################

# =================================================================
## Minimalistic construction to allow adding custom attributes (e.g. units) to an ndarray sub class
#  see http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
#  Define hooks __new__, __array_finalize__
#
class InfoArray(ndarray):
    ## hook required by numpy
    def __new__(cls, input_array, info=None):
        obj = asarray(input_array).view(cls)
        obj.info = info
        return obj
    ## hook required by numpy
    def __array_finalize__(self, obj):
        if obj is None: return
        self.info = getattr(obj, 'info', None)

