#!/usr/bin/env python
# -*- coding: utf-8 -*-
##################################################################################################
##       @package constants
#  Misc Geophysical constants
##################################################################################################
from numpy import pi
EarthMeanRadius = 6371.0088e3  # m, International Union of Geodesy and Geophysics (IUGG) 
deg2rad         = pi/180.0     # degree to radian conversion factor
g               = 9.80665      # m/s2, nominal gravitational acceleration of an object in a vacuum near the surface of the Earth

seconds_per_sidereal_day = 86400 - 4*60 + 4.1             # 23 hr 56 m 4.1 s
omega                    = 2*pi/seconds_per_sidereal_day  # rotation rate of the Earth = 2pi per sidereal day
