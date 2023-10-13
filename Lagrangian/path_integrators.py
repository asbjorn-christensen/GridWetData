#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################################################################
##       @package path_integrators
#        Light weight Lagrangian classes for 2D/3D individual-based modelling
#
#        Emphasis on flexibility and transparency rather than optimal speed (but
#        having performance in mind)
#
#        TODO: integrators with turbulent step 
#######################################################################################################################

from lagrangian_base import *
from GridWetData.constants import EarthMeanRadius, deg2rad # GridWetData directory must be visible in path

_verbose         = False # True  # log info

# ===============================================================================
#             Error sub structure of path integrators
# ===============================================================================


# ===============================================================================
# Integrator protocol: all path integrators should apply this calling sequence
#
# newx =  integrator(x, env, now, dt)
#
# x   : is the current position of the particle (2D/3D)
# env : is the environmental query object
# now : is the current time (datetime like)
# dt  : is the time step (timedelta like) - may be positive or negative (backtracking)
#
# integrator should generate the new position newx corresponding to
# time now + dt, without considering boundary conditions (currently we
# split propagation and boundary condition application - this paradigm may
# change in the future).
# Advection field is accessed like env.currents(at_pos, at_time) and
# expected output is a Cartesian 2D/3D vector oriented along (longitudeE, latitudeN [, depth_down]) in units m/2
# Dimensionality of Advection field and x must match (otherwise a numerical exception is raised)
# ===============================================================================

## add sec seconds to a timedelta - sec is a float/int
def addsecs(timeobj, sec):
    return timeobj + datetime.timedelta(seconds=sec)
#  ------------------------------------------------------------------------------------------------
## Transform Cartesian current vector u at position x in units m/s to units (degE/s, degN/s [,m/s])
#  
#  works for horizontal current vectors (2D) and full 3D current vectors
#  NB: inplace transformation of u is performed
#  ------------------------------------------------------------------------------------------------
def to_deglonlat_per_sec(u,x):
    u[0] /= EarthMeanRadius*cos(x[1]*deg2rad)*deg2rad  # bugfix 13 Oct 2023
    u[1] /= EarthMeanRadius*deg2rad # x = (lambda[degE], phi[degN], z[m_down])
    return u
    
def euler(x, env, now, dt):
    h  = dt.total_seconds()
    u1 = to_deglonlat_per_sec(env.currents(x, now), x)
    return x + u1*h

def rk2(x, env, now, dt):
    h  = dt.total_seconds()
    u1 = to_deglonlat_per_sec(env.currents(x, now), x)
    t2 = addsecs(now, 0.5*h)
    x2 = x+0.5*u1*h
    u2 = to_deglonlat_per_sec(env.currents(x2, t2), x)
    return x + u2*h

def rk4(x, env, now, dt):
    h  = dt.total_seconds()
    t2 = addsecs(now, 0.5*h)
    t4 = addsecs(now,     h)
    u1 = to_deglonlat_per_sec(env.currents(x           , now), x)
    u2 = to_deglonlat_per_sec(env.currents(x + 0.5*h*u1, t2 ), x)
    u3 = to_deglonlat_per_sec(env.currents(x + 0.5*h*u2, t2 ), x)
    u4 = to_deglonlat_per_sec(env.currents(x +     h*u3, t4 ), x)
    return x + h*(u1 + 2*(u2 + u3) + u4)/6
