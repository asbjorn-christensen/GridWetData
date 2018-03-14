#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################################################################
##       @package particles
#        Light weight Lagrangian classes for 2D/3D individual-based modelling
# 
#        Emphasis on flexibility and transparency rather than optimal speed (but
#        having performance in mind)
#
#        
#######################################################################################################################

from boundary_conditions import NonSticky # to be amended
from path_integrators import euler, rk2, rk4
from lagrangian_base import *

_verbose         = False # True  # log info

# ===============================================================================
#             Error sub structure of particles
# ===============================================================================


# ===============================================================================
## Super class for particles which inherits from a passive particle
#
# Main attributes
#    pos         : 2D/3D position (numpy vector)
#    env         : particle environment providing environmental query protocol
#    motion_mask : 2D/3D vector with 0/1 for projecting steps (numpy vector)
#    bc          : boundary condition object
# ===============================================================================
class PassiveParticle:
    #  -----------------------------------------------------
    ## constructor
    #  @param self   The object pointer.
    #  @param pos    current position (2D/3D)
    #  @param env    Particle environment object providing environmental query protocol
    def __init__(self, pos, env, **kwargs):
        self.pos = 1.0*asarray(pos, float) # force copy
        self.env = env
        #
        if kwargs.has_key("motion_mask"):
            self.motion_mask = asarray(kwargs["motion_mask"], float)
        else:
            self.motion_mask = ones(len(self.pos), float) # default: free particle
        #
        if kwargs.has_key("boundary_conditions"):
            self.bc = kwargs["boundary_conditions"]
        else:
            self.bc = NonSticky(self.env)   # default to be amended
        #
        if kwargs.has_key("integrator"):
            self.integrator = kwargs["integrator"]
        else:
            self.integrator = rk2
        
    ## -----------------------------------------------------
    ## Conform a new proposed particle position to boundary conditions and topography
    #
    #  @param self      The object pointer.
    #  @param newpos    New proposed position (2D/3D)
    #  @return corrpos  New position conforming to boundary conditions and topography
    def correct_step(self, newpos, attime):
        dpos   = newpos - self.pos # assume linear space
        dpos  *= self.motion_mask
        dpos   = self.bc.validate_step(self.pos, dpos, attime, self.env)
        return dpos
    
    ## -----------------------------------------------------
    ## Unified 2D/3D spatial update of passive particle
    #
    #  @param self   The object pointer.
    #  @param now    current time (datetime like object)
    #  @param dt     time step (timedelta like object)
    def update(self, now, dt):
        newpos = self.integrator(self.pos, self.env, now, dt)
        dpos   = self.correct_step(newpos, now+dt)
        self.pos += dpos

        
#######################################################################
#
#   Particles with properties beyond that of passive particles
#   are represented as sub classes of PassiveParticle
#
#######################################################################

# ===============================================================================
## Passive particle with age
#
# Example on biological decoration of PassiveParticle by subclassing
# ===============================================================================
class AgedParticle(PassiveParticle):
    ## -----------------------------------------------------
    ## constructor
    #  @param self   The object pointer.
    #  @param pos    current position (2D/3D)
    #  @param env    Particle environment object providing environmental query protocol
    def __init__(self, pos, env, age=0, **kwargs):
        self.age = age
        PassiveParticle.__init__(self, pos, env, **kwargs)

    ## -----------------------------------------------------
    ## Update of particle state
    #
    #  @param self   The object pointer.
    #  @param now    current time (datetime like object)
    #  @param dt     time step (timedelta like object)
    def update(self, now, dt):
        self.age += dt
        PassiveParticle.update(self, now, dt)


# ===============================================================================
## Particle with swimming behavior
#
# Example on biological decoration of PassiveParticle by subclassing
# ===============================================================================
class SwimmingParticle(PassiveParticle):
    ## -----------------------------------------------------
    ## constructor
    #  @param self   The object pointer.
    #  @param pos    current position (2D/3D)
    #  @param env    Particle environment object providing environmental query protocol
    #  @param swim   function describing spatial swimming/active behavior
    def __init__(self, pos, env, swim, **kwargs):
        self.swim = swim
        PassiveParticle.__init__(self, pos, env, **kwargs)

    ## -----------------------------------------------------
    ##  Update of particle state
    #
    #  Apply active + passive motion components sequentially
    #  @param self   The object pointer.
    #  @param now    current time (datetime like object)
    #  @param dt     time step (timedelta like object)
    def update(self, now, dt):
        newpos = self.swim(self, now, dt)
        dpos   = self.correct_step(newpos)
        self.pos += dpos
        PassiveParticle.update(self, now, dt)
