#!/usr/bin/env python

#####################################################################################
#  Create grid descriptor of HBM/ERGOM wetpoint formatted data
#
#  usage:     make_grid_desc.py   <grid_summary_file>   <grid_description>   <netCDF_descriptor>  
#
#  grid_summary_file contains the ASCII header:
#
#       nx   ny   nz   # grid dimensions (longitude,latitude,vertical)
#       lon0  lat0     # South-west corner cell center (longitude,latitude) in degrees E,N
#       dlon  dlat     # grid spacings along  (longitude,latitude) in decimal degrees
#
#  grid_description is the DMI ASCII formatted grid description file containing 
#  wet point map and static cell thicknesses (dry cells are given thickness = -1)
#  maintain DMIs original wet point format, where first element is a pad value for dry points
#
#  Based on the netCDF4 API
#
#  make_grid_desc.py   ../Grids/ns_grid.sum    /home/data/GUDP-VIND_test/VINDdata/setup/ns_grid_depth_v5   ns_grid.nc
#
#####################################################################################
from   numpy    import *
import netCDF4  as netcdf 

import sys
import exceptions
import datetime

usage_str = "usage: %s  <grid_summary_file>   <grid_description>  <netCDF_descriptor> " % sys.argv[0]

if len(sys.argv) != 4:
    raise exceptions.IOError(("%s expects 3 arguments\n\n" % sys.argv[0]) + usage_str)
#
# ------------ read header file ------------
#
fin          = open(sys.argv[1])
nx, ny, nz   = map(int,   fin.readline().split()[:3])
dlon, dlat   = map(float, fin.readline().split()[:2])
lon0, lat0   = map(float, fin.readline().split()[:2])
fin.close()

#
# ------------ read grid_description ------------
#

fgr    = open(sys.argv[2])
data   = fgr.read().split()

# --- parse 3D to 1D map ---

n3d = nx * ny * nz
map3d_to_1d = array(map(int, data[:n3d]), int)     # m[nz,nx,ny],  ix,iy,iz offset 0
map3d_to_1d = reshape(map3d_to_1d, (nz, nx, ny))   # fastest index opposite fortran, fortran declaration order (ny,nx,nz)
#
map3d_to_1d = swapaxes(map3d_to_1d, 0,1)   # m[nz,nx,ny] -> m[nx,nz,ny]
map3d_to_1d = swapaxes(map3d_to_1d, 1,2)   # m[nx,nz,ny] -> m[nx,ny,nz]
map3d_to_1d = map3d_to_1d[:, ::-1, :]      # flip y-axis (native grid scans N->S)

nwet3d = int(1 + ma.max(map3d_to_1d))         # incl element 0 is pad value
nwet2d = int(1 + ma.max(map3d_to_1d[:,:,0]))  # incl element 0 is pad value; data scanned layer wise, starting at surface layer

# --- read cell thickness table ---

cwd = array(map(float, data[n3d:]), float)       # include pad value = data[n3d]
cell_thickness = -ones((nx,ny,nz), float)        # pad with -1 to signal dry cells

# print "nwet3d = %d should len(cwd) = %d" % (nwet3d, len(cwd)) 

for ix in range(nx):
    for iy in range(ny):
        for iz in range(nz):
            i = map3d_to_1d[ix,iy,iz]
            if i > 0: # i=0 is dry points
                cell_thickness[ix,iy,iz] = cwd[i]
                
##                 if iz==0:
##                     print ix,iy,cwd[i]

## #
## # ------ print wet points (ix,iy) ------
## #               
## for ix in range(nx):
##     for iy in range(ny):
##         if map3d_to_1d[ix,iy,0] >= 0:
##             print ix,iy

## for ix in range(nx):
##      for iy in range(ny):
##         if map3d_to_1d[ix,iy,0] > 0:
##             print cell_thickness[ix,iy,0]       

#
# ------------ dump to netCDF ------------
#

ncfile    = netcdf.Dataset(sys.argv[3], 'w')               # netCDF4 constructor

# --- define mode ---

ncfile.createDimension("nx", nx)
ncfile.createDimension("ny", ny)
ncfile.createDimension("nz", nz)
ncfile.createDimension("nwet3d", nwet3d)
ncfile.createDimension("nwet2d", nwet2d)

ncfile.createVariable("lon0", 'd', ())
ncfile.createVariable("lat0", 'd', ())
ncfile.createVariable("dlon", 'd', ())
ncfile.createVariable("dlat", 'd', ())

ncfile.createVariable("map3d_to_1d",   'i',  ("nx", "ny", "nz"))
ncfile.createVariable("mapsurf_to_1d", 'i',  ("nx", "ny"))
ncfile.createVariable("cell_thickness", 'd', ("nx", "ny", "nz"))

bcmd  = "\n creation  : "    + " ".join(sys.argv)
btime = "\n created at: "  + datetime.datetime.now().strftime("%Y - %m - %d")
ncfile.description = "grid descriptor file for HBM/ERGOM wetpoint grid from underlying regular lon-lat grid" + bcmd + btime
ncfile.variables["dlon"].unit            = "decimal degrees"
ncfile.variables["dlat"].unit            = "decimal degrees"
ncfile.variables["lon0"].unit            = "decimal degrees E"
ncfile.variables["lat0"].unit            = "decimal degrees N"
ncfile.variables["cell_thickness"].unit  = "meters"

ncfile.variables["dlon"].description           = "longitude step"
ncfile.variables["dlat"].description           = "latitude step"
ncfile.variables["lon0"].description           = "SW corner cell center, longitude" 
ncfile.variables["lat0"].description           = "SW corner cell center, latitude" 
ncfile.variables["cell_thickness"].description = "vertical cell thickness (without dynamic sea level elevation). Dry cells has negative thickness"
ncfile.variables["map3d_to_1d"].description    = "bulk map: field_3d(i,j,k) <- field_1d(map3d_to_1d(i,j,k)) - field_1d(0) is pad value for dry points"
ncfile.variables["mapsurf_to_1d"].description  = "surface layer map: field_2d(i,j) <- field_1d(map2d_to_1d(i,j)) - field_1d(0) is pad value for dry points"

ncfile.variables["cell_thickness"].padvalue = -1

# --- assign mode ---

ncfile.variables["lon0"].assignValue(lon0)
ncfile.variables["lat0"].assignValue(lat0)
ncfile.variables["dlon"].assignValue(dlon)
ncfile.variables["dlat"].assignValue(dlat)
ncfile.variables["map3d_to_1d"][:]    = map3d_to_1d            # netCDF4 array assignment
ncfile.variables["mapsurf_to_1d"][:]  = map3d_to_1d[:,:,0]     # netCDF4 array assignment
ncfile.variables["cell_thickness"][:] = cell_thickness         # netCDF4 array assignment
ncfile.close()
 
