#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################################################################
##       @package boundary_conditions
#        Boundary conditions for light weight Lagrangian classes for 2D/3D individual-based modelling
#
#       
#
#        TODO: ReflectiveBoundaryConditions : canonical horizontal (and vertical, if 3D)
#                                             reflective boundary conditions
#                                             
#######################################################################################################################

from lagrangian_base import *

_verbose         = False # True  # log info

# ===============================================================================
#             Error sub structure of boundary conditions 
# ===============================================================================


# ===============================================================================
#             Protocols for boundary condition objects
#
# __init__      : must accept the environmental context object
# validate_step : first  arg == current 2D/3D position
#                 second arg == spatial additive 2D/3D step from current position to be explored
#                 return: corrected step
# ===============================================================================

# -------------------------------------------------------
## Solve the vertical bouncing problem for
#  0 < z < wd where dz is proposed sted
#  @param  z     vertical current position
#  @param  wd    water depth
#  @param  dz    proposed step (oriented as z)
#  @return dzcorr  corrected step, corresponding to z+dz bounced bewteen surface and bottom (if needed)
# -------------------------------------------------------
def vertical_bounce(z, wd, dz):
    assert 0 <= z <= wd
    # znew is the position z+dz, bounced if needed
    nb, znew = divmod(abs(z+dz), wd)
    nb   = int(nb+0.5)
    # resolve whether last bounce is surface or bottom and set znew accordingly
    if nb%2 == 1:
        znew = wd-znew # bottom bounce
    dzcorr = znew-z  # convert to corrected step 
    return  dzcorr

# -------------------------------------------------------
## generate a random point in a circle with radius r
# -------------------------------------------------------
def random_point_in_circle(r):
    u   = random.rand(2)
    rad = r*sqrt(u[0]) # correct radial distribution
    theta = 2*pi*u[1]
    x   = cos(theta)
    y   = sin(theta)
    return rad*array([x,y])

def random_point_in_sphere(r):
    u       = random.rand(3)
    rad     = r * u[0]**(1./3.)  # correct radial distribution
    theta   = 2*pi*u[1]
    phi	    = arccos(2*u[2] - 1) # correct meridional distribution
    x       = sin(phi)*cos(theta)
    y       = sin(phi)*sin(theta)
    z       = cos(phi)
    return rad*array([x,y,z])


# -------------------------------------------------------
## dimension independent wet point generator
#
# search in successively smaller spheres corresponding to
# proposed step dx direction
# This algorithm may fail, if x becomes dried at time t
# (no wet convergence point)
# -------------------------------------------------------
def locate_wetpt(x, dx0, t, wet_test, randomptgen, nout=100, ninner=3):
    dx = 1.0*dx0 # force copy
    for i in range(nout):
        for j in range(ninner):
            d  = 0.5*dx
            rd = sqrt(sum(d*d))
            dxtest = d + randomptgen(rd)
            if wet_test(x+dxtest, t):
                return dxtest
        else:
            dx *= 0.5 # reduce radius of search sphere
    else:
        dim = len(x)
        raise WetPointUnresolvable("locate_wetpt(%dD) bailiing out" % dim)


# ===============================================================================
## Pesudo-reflective boundary conditions
#
#  Poor mans version of ReflectiveBoundaryConditions - to be developed and optimized
# ===============================================================================
class NonSticky:
    def __init__(self, env):
        self.env = env
        
    #  ---------------------------------------------------------------------------------
    ## Validate a proposed 2D/3D step dpos from current position currpos
    #
    #  This crude placeholder algorithm is not guaratied to produce
    #  appropriate spatial distribution near boundaries
    #  No path checking (step may cross land)
    #  Basic strategy:
    #       1) correct trial step so it stays wet by stochastic sampling
    #       2) horizontal/vertical problem split if 3D (rely on small-slope bottum)
    #  @param  self    The object pointer
    #  @param  currpos    current position
    #  @param  dpos       proposed step
    #  @param  attime     time corresponding to after proposed step (when supposedly pos = currpos+dpos)
    #  @param  env        environmental query object 
    #  @return corr_dpos  corrected step
    def validate_step(self, currpos, dpos, attime, env):
        assert len(currpos) == len(dpos)
        assert env.is_wet(currpos, attime) # test for dried base point
        is3D = len(dpos)==3
        # apply vertical BC first, if 3D 
        if is3D:
            dpos[2] = vertical_bounce(currpos[2], env.wdepth(currpos, attime), dpos[2])
        # apply horizontal BC (3D: on vertically prefiltered step)
        if env.is_wet(currpos+dpos, attime):
            return dpos # current step is OK
        else:
            if is3D:
                return locate_wetpt(currpos, dpos, attime, self.env.is_wet, random_point_in_sphere)
            else:
                return locate_wetpt(currpos, dpos, attime, self.env.is_wet, random_point_in_circle)
     
            
            
if __name__ == "__main__":
    ## -------------------------------
    #wd = 2.0
    #for dz in arange(-6, 10, 0.01):
    #    print dz, vertical_bounce(0.5, wd, dz)
    for i in range(1000000):
        xyz = random_point_in_sphere(2)
        if abs(0.5 < xyz[2] < 0.51): print "%12.7f %12.7f" % tuple(xyz[:2])
