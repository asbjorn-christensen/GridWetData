#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################################################################
##       @package ensembles
#        Light weight Lagrangian classes for 2D/3D individual-based modelling
#
#        Emphasis on flexibility and transparency rather than optimal speed (but
#        having performance in mind)
#
#        TODO:  Ensemble.evaluate_2D_distribution(self, agrid):
#######################################################################################################################


import particles
from lagrangian_base import *
from GridWetData.grid_data import GridData_2D # for Ensemble.evaluate_2D_distribution

_verbose         = False # True  # log info

# ===============================================================================
#             Error sub structure of ensembles
# ===============================================================================



# ===============================================================================
#             Protocols for environmental objects
#
# currents(x,t) : (2D/3D contexts) interpolate the currents at (x,t)
# wdepth(x,t)   : (in 3D contexts) interpolate the water depth at (x,t)
# is_wet(x,t)    : (2D/3D contexts) assess wheter (x,t) is a wet point

# coast line crossing: to be implemented
# ===============================================================================


# -------------------------------------------------------
## Represent an ensemble of particles (may be instances of different classes)
#
#  Important attributes of an ensemble instance:
#    env :        environmental query object
#    now :        current time, in relation to particle propagation (datetime like)
#    generators:  (optional) list of particle generators that should be invoked during particle propagation 
# -------------------------------------------------------

class Ensemble(list):
    #  ---------------------------------------------------------------------------------
    ## Advance particles in ensemble by time step
    #
    #  @param self    The object pointer
    #  @param dt      time step (timedelta like) to propagate ensemble from current time
    #  
    def propagate(self, dt):
        # --- update each particle ---
        for part in self:
            part.update(self.now, dt)
        # --- activate generators ---
        if hasattr(self, "generators"):
            for gen in self.generators:
                gen(self) # may inject particles
        self.now += dt # advance clock
   
    #  ---------------------------------------------------------------------------------
    ## Probe positions of particles in the ensemble
    #  
    #  @param self     The object pointer
    #  @return parpos  Postions of particles (2D or 3D, depending of instantiation)
    #                  Mixture of 2D/3D particles will currently throw an exception
    def get_positions(self):
        parpos = []
        for par in self:
            parpos.append(par.pos)
        # cast will throw an exception if 2D/3D particles are mixed
        return asarray(parpos)
    
    #def evaluate_3D_distribution(self, agrid):
    
    #  ---------------------------------------------------------------------------------
    ## Create 2D distribution of particles on provided grid
    #
    #  If ensemble is 3D, vertical integration is understood
    #  @param self      The object pointer
    #  @param agrid     The grid, where 2D distribution should be evaluated
    #  @return          a griddata instance, corresponding to agrid
    #
    #  agrid must feature: dimension attributes nx,ny
    #                      method get_cell_indices(lonlatpos)
    #                      method get_cell_areas
    def evaluate_2D_distribution(self, agrid):
        distrib = zeros((agrid.nx, agrid.ny), float)
        for par in self:
            ix,iy = agrid.get_cell_indices(par.pos)[:2] # waste vertical component
            distrib[ix,iy] += 1
        distrib /= agrid.get_cell_areas()
        agrid_2Dcast = agrid.export_horizontal_grid() # may be agrid, if 2D
        return GridData_2D(agrid_2Dcast, distrib)

    
###################################################################################
#            Ensemble generators / factories
###################################################################################

#  ---------------------------------------------------------------------------------
## Handy particle factory to inject a number of particles into ensemble
#
#  Works, if particle is instantiated like a passive particle
#  @param ensb     ensemble, where environmental object has been set. injection inplace
#  @param npart    number of particles to inject
#  @param edge0    lower 2D/3D edge of coordinates of starting box
#  @param edge1    upper 2D/3D edge of coordinates of starting box
#  @param pargen   particle class; must be instantiated like pargen(position, env)
#  @param maxfail  maximum number of attempts to pick particle position in starting box per particle
#  
def uniform_distribution_of_wetpts(ensb, npart, edge0, edge1,
                                   pargen=particles.PassiveParticle, maxfail=1000):
    diff = edge1-edge0
    for ipar in range(npart):
        for itry in range(maxfail):
            pos = edge0 + random.rand(len(edge0))*diff
            if ensb.env.is_wet(pos, ensb.now):         
                ensb.append(pargen(pos, ensb.env))
                break
        else: # we exhausted allowed number of trails
            raise WetPointUnresolvable("range cap exceeded")
    
