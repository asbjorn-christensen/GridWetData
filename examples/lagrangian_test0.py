#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################################################################
## Usage example of GridWetData.Lagrangian
## Remember to set your python path, e.g. (bash)
##   export PYTHONPATH="/home/asbjorn/tmp10/":${PYTHONPATH}   
#######################################################################################################################

from GridWetData.Lagrangian import *
from GridWetData.Lagrangian.test_environments import StaticFlow3D
from GridWetData.Lagrangian.lagrangian_base import *   # numpy etc

# --- setting up a particle ensemble 
myens = ensembles.Ensemble()                               # ensemble instance
myens.env = StaticFlow3D(array([1.,1.,0.], float), 20.0)   # artificial test hydrography
myens.now = datetime.datetime(2005,1,3)                    # set ensemble clock

# --- creating 100 particles in the ensemble
edge0 = array([3.0, 55, 2], float)                                 # corner of box, where particles are released           
ensembles.uniform_distribution_of_wetpts(myens, 100, edge0, edge0) # release particles in corner (default: passive particles)

# --- running a simulation with 100 time steps forward
for i in range(100):
    myens.propagate(datetime.timedelta(seconds=300))    # use time step of 300 sec

# --- display particle positions    
print myens.get_positions()  # numpy array

# --- cast particle distribution on a grid 
from GridWetData.grids import LonLatGrid
agrid = LonLatGrid(100, 100, 3, 55, 1./6, 0.1)          # plotting grid
gd2d = myens.evaluate_2D_distribution(agrid)            # GridData_2D instance
gd2d.grid.write_data_as_netCDF("jj.nc", gd2d.data)      # example of saving particle distribution to netcdf

