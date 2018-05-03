#!/usr/bin/env python
# -*- coding: utf-8 -*-
########################################################################
##   @package derived_layers
#    External layer-generation utilities
#######################################################################
from numpy import *
#  ---------------------------------------------------------------------------------
## Generate a stratification front index of any griddata object
#  Relatively crude algorithm algortihm, looking at absolute surface-bottom
#  difference, and taking an isotropic, non-cerntered absolute gradient
#  TODO: Apply a natural scale to the index; currently, the frontal index is not scaled
#        Gradient evaluation should be done with a generic method, e.g. a projection opeartor provided by the grid
#  @param  obj A GridData_3D like object
#  @param  padval (default = 0) applies to polits, where front index is not defined (e.g. the grid rim)
#  @return A layer dstrat[nx,ny] with frontal index
#  
def StratificationFront(obj, padval = 0):
    tsurf = obj.get_surface_layer() # GridData_2D
    tbott = obj.get_bottom_layer()  # GridData_2D
    strat = abs(tsurf.data-tbott.data) 
    dstrat = padval*ones(strat.shape, float) # padval -> rim 
    # assume strat mesh is sufficiently isotropic
    dstrat[1:-1, 1:-1] = abs(strat[0:-2, 1:-1]-strat[1:-1, 1:-1]) + \
                         abs(strat[2:  , 1:-1]-strat[1:-1, 1:-1]) + \
                         abs(strat[1:-1, 0:-2]-strat[1:-1, 1:-1]) + \
                         abs(strat[1:-1, 2:  ]-strat[1:-1, 1:-1])
    return tsurf.__class__(tsurf.grid, dstrat) # inherit grid

