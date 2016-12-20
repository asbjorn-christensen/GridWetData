#!/usr/bin/env python

from GridWetData.grids import *

# =================================================================
## Mix-in class for shared 2D/3D operations
class HBM_wetpoint_data:
    # --------------------------------------------------------------------
    ## Generate a informative string representation of object
    #  @param self The object pointer
    #
    def __str__(self):
        return "%s instance generated from file %s" % (self.__class__.__name__, self.fname)

    #  --------------------------------------------------------------
    ## Load wet point frame from a netcdf file and inflate to regular array
    #  set attribute unit expected in data set
    #  @param self  The object pointer
    #  @param fname file name of netCDF data
    #
    def load_data_frame(self, fname):
        ncfile      = netcdf.Dataset(fname)
        assert self.nwet == ncfile.dimensions["nwet"].size
        wetptdata   = ncfile.variables["data"][:]   # in Scientific.IO.NetCDF getValue() should be applied to retrieve array
        data        = InfoArray(self.inflate_array(wetptdata)) # 2D/3D sub class defines inflate_array
        ## set unit from file meta data
        data.unit     = ncfile.unit
        ## set data pad value 
        data.padvalue    = wetptdata[0]        # specific convention for this format
        data.minwetvalue = wetptdata[1:].min() # specific convention for this format
        data.maxwetvalue = wetptdata[1:].max() # specific convention for this format
        ncfile.close()
        return data

# =============================================================
## 2D surface wet point HBM grid
class HBMGrid_2D(LonLatGrid, HBM_wetpoint_data):
    #  -----------------------------------------------------
    ## constructor from file data
    #  @param self        The object pointer
    #  @param grid_desc   file name of grid descriptor in netCDF format
    #
    def __init__(self, grid_desc):
        ncfile     = netcdf.Dataset(grid_desc)
        ## file name corresponding to loaded data set
        self.fname = grid_desc
        nx        = ncfile.dimensions["nx"].size # integer 
        ny        = ncfile.dimensions["ny"].size # integer 
        ## number of wet grid point
        self.nwet = ncfile.dimensions["nwet2d"].size # integer 
        lon0      = ncfile.variables["lon0"].getValue() # Retrieve a scalar value from a scipy.io.netcdf.netcdf_variable
        lat0      = ncfile.variables["lat0"].getValue() # Retrieve a scalar value from a scipy.io.netcdf.netcdf_variable
        dlon      = ncfile.variables["dlon"].getValue() # Retrieve a scalar value from a scipy.io.netcdf.netcdf_variable
        dlat      = ncfile.variables["dlat"].getValue() # Retrieve a scalar value from a scipy.io.netcdf.netcdf_variable
        ## map buffer to 2 array: data[ix,iy] = buffer[map[ix,iy]]
        self.map  = ncfile.variables["mapsurf_to_1d"][:]     # in Scientific.IO.NetCDF getValue() should be applied to retrieve array
        ncfile.close()
        #
        LonLatGrid.__init__(self,  nx, ny, lon0, lat0, dlon, dlat)

    #  --------------------------------------------------
    ## raw inflate of 2D wetpoit grid to generic 2D float ndarray
    #  pad value for dry points is contained in first element of wetptdata
    #  @param self        The object pointer
    #  @param wetptdata   wet point vector
    #  @param astype      optional cast type of ndarray (default float)
    #  @return            inflated 2D array
    #
    def inflate_array(self, wetptdata, astype=float):
        data = wetptdata[0]*ones((self.nx, self.ny), astype) # set pad value
        for ix in range(self.nx):
            for iy in range(self.ny):
                    data[ix,iy] = wetptdata[self.map[ix,iy]]
        return data


#  =====================================================================
## 3D wet point HBM grid
#  compositional relation to the 2D surface wet point cmod grid
#         
class HBMGrid_3D(LonLatZGrid, HBM_wetpoint_data):
    #  -----------------------------------------------------
    ## constructor from file data
    #  @param self        The object pointer
    #  @param grid_desc   file name of grid descriptor in netCDF format
    #
    def __init__(self, grid_desc):
        ncfile     = netcdf.Dataset(grid_desc)
        ## file name corresponding to loaded data set
        self.fname = grid_desc
        nx        = ncfile.dimensions["nx"].size   # integer 
        ny        = ncfile.dimensions["ny"].size   # integer 
        nz        = ncfile.dimensions["nz"].size   # integer
        ## number of wet grid point
        self.nwet = ncfile.dimensions["nwet3d"].size   
        lon0      = ncfile.variables["lon0"].getValue()  # Retrieve a scalar value from a scipy.io.netcdf.netcdf_variable
        lat0      = ncfile.variables["lat0"].getValue()  # Retrieve a scalar value from a scipy.io.netcdf.netcdf_variable
        dlon      = ncfile.variables["dlon"].getValue()  # Retrieve a scalar value from a scipy.io.netcdf.netcdf_variable
        dlat      = ncfile.variables["dlat"].getValue()  # Retrieve a scalar value from a scipy.io.netcdf.netcdf_variable
        ## map buffer to 3D array: data[ix,iy,iz] = buffer[map[ix,iy,iz]]
        self.map  = ncfile.variables["map3d_to_1d"][:]     # in Scientific.IO.NetCDF getValue() should be applied to retrieve array
        cellw0    = ncfile.variables["cell_thickness"][:]  # in Scientific.IO.NetCDF getValue() should be applied to retrieve array
        ncfile.close()
        #
        LonLatZGrid.__init__(self,  nx, ny, lon0, lat0, dlon, dlat, cellw0)
        self.surf_grid = HBMGrid_2D(grid_desc) # compositional relation, no heritance

    #  --------------------------------------------------
    ## raw inflate of 3D wetpoit grid to generic 3D float ndarray
    #  pad value for dry points is contained in first element of wetptdata
    #  @param self        The object pointer
    #  @param wetptdata   wet point vector
    #  @param astype      optional cast type of ndarray (default float)
    #  @return            inflated 3D array
    #
    def inflate_array(self, wetptdata, astype=float):
        data = wetptdata[0]*ones((self.nx, self.ny, self.nz), astype) # set pad value
        for ix in range(self.nx):
            for iy in range(self.ny):
                for iz in range(self.nz):
                    data[ix,iy,iz] = wetptdata[self.map[ix,iy,iz]]
        return data

    # -----------------------------------------------------------
    ## Update dynamic grid attributes corresponding to sea level in file zname
    #  @param self    The object pointer
    #  @param zname   file name of wetpoint vector for sea level elevation in netCDF format
    #
    def update_sea_level(self, zname):
        z = self.surf_grid.load_data_frame(zname)
        self.set_reference_level(z) # super class handles dynamic attributes
        

        
#################### self test #########################
#
#  mid North Sea test point: ix,iy = 100,150
#
if __name__ == "__main__":
    _verbose = True
    ## g2D = HBMGrid_2D(os.path.join("../DMI_data", "ns_grid.nc"))
    ## print "dir(g2D)=", dir(g2D)                 
    g3D = HBMGrid_3D(os.path.join("../DMI_data", "ns_grid.nc"))
    ## g3D.update_sea_level(os.path.join("../DMI_data", "z_ns_2015_05_20_00_15_00.nc"))
    ## for ix in range(g3D.nx):
    ##      for iy in range(g3D.ny):
    ##          if g3D.wetmask[ix,iy,0] == 1: print g3D.cellw[ix,iy,0]
    ## print "cellw0      =",g3D.cellw[100,150,:]  
    ## print "ccdepth0    =",g3D.ccdepth[100,150,:] 
    ## print "layer_sep0  =",g3D.layer_sep[100,150,:] 
    ## for ix in range(g3D.nx):
    ##     for iy in range(g3D.ny):
    ##         if  g3D.wetmask[ix,iy,0] == 0: print ix,iy
    ##
    ## for x in arange(-3,9,0.1):
    ##     for y in arange(51,59,0.1):
    ##         if not g3D.is_wet((x,y,0.1)):  print x,y
    ##
    ##g3D.update_sea_level(os.path.join("../DMI_data", "z_ns_2015_05_20_00_15_00.nc"))
    ## for x in arange(-3,9,0.001):
    ##     pos = (x,54)
    ##     wd = g3D.interpolate_wdepth(pos)
    ##     ixc,iyc = g3D.surf_grid.get_cell_indices(pos)
    ##     if wd != None:
    ##         print x,wd, g3D.bottom_layer[ixc,iyc]
                
                    
         
