#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################################################################
##       @package test_environments
#        Internal testing environments for Lagrangian setup (not included in __init__.py)
#######################################################################################################################


from lagrangian_base import *

_verbose         = False # True  # log info

class StaticFlow3D:
    def __init__(self, usurf, wd):
        self.usurf     = zeros(3, float)
        self.usurf[:2] = usurf[:2] # surface current, ignore potential 3rd component
        self.wd        = wd        # water depth
        
    def currents(self,x,t):
        if self.is_wet(x,t):
            return self.usurf*(1.0 - x[2]/self.wd) # zero at bottom
        else:
            return zeros(3, float)
        
    def wdepth(self,x,t):
        return self.wd

    def is_wet(self,x,t):
        return 0 <= x[2] <= self.wd
        
