#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################################################################
#    Setup script to install package - currently no options
#
#    TODO: make DTU10 build optional
#          DTU10 optional compiler for f2py
#          optional tidal constituents
#######################################################################################################################
import sys
import os
import exceptions
import platform

# ------------------------------------------------------
#               Install DTU10 API
# 
# Linux: fortran compiler currently relies on f2py defaults
# ------------------------------------------------------
if platform.system() == 'Linux':
    os.system("f2py -c -m dtu10API fortran_sources/perth3.f")                                    # build dtu10API.so
    os.system("cat data/NorthSea_tidal_constituents.dat.gz | gunzip > tidal_constituents.dat")   # tidal constituents
else:
    raise exceptions.OSError("non-linux install cnot implemented")
