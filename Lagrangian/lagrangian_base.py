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

import exceptions
import os
import datetime
from numpy import *

_thisdir         = os.path.dirname(__file__)  # allow remote import
_verbose         = False # True  # log info

# ===============================================================================
#             Error structure 
# ===============================================================================

## Super class for errors originating from Lagrangian context

class LagrangianError(exceptions.Exception):
    pass

# many algorithms in the Lagrangian context relies on
# resolving a wet point by stochastic / deterministic processes, which may fail
class WetPointUnresolvable(LagrangianError):
    pass


