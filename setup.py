#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################################################################
#    Setup script to install package (non-python components)- currently no options
#
#    TODO: make DTU10 build optional
#          DTU10 optional compiler for f2py
#          optional tidal constituents
#          if verbose, create install log
#######################################################################################################################
import sys
import os
import exceptions
import platform
from subprocess import Popen, PIPE

# ------------------------------------------------------
#               Install DTU10 API
# 
# Linux: fortran compiler currently relies on f2py defaults
#        expect gunzip available
# ------------------------------------------------------

if platform.system() == 'Linux':
    # --- build dtu10API.so 
    h = Popen("f2py -c -m dtu10API fortran_sources/perth3.f", stdout=PIPE, stderr=PIPE, stdin=PIPE, shell=True)
    stdout, stderr = h.communicate()
    # --- install tidal constituents
    h = Popen("cat data/NorthSea_tidal_constituents.dat.gz | gunzip > tidal_constituents.dat", stdout=PIPE, stderr=PIPE, stdin=PIPE, shell=True)
    stdout, stderr = h.communicate()
else:
    raise exceptions.OSError("non-linux install not implemented")


print 63*"-"
print 10*">", " Installation of GridWetData successful  ", 10*"<"
print 63*"-"
