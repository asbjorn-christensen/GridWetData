## \mainpage GridWetData
# GridWetData is a package to access data in a water body, which is described by an underlying grid
# embedding the spatial area of the water body
#
# \image latex tfront.png "Output example: frontal structure in the North Sea" width=\textwidth
#
# <br>
# Copyright 2015-2016 Asbjorn Christensen
# <br>
# ---------------------------------------------------------------------------------------------------------------
# <br>
# LICENSE:
#    GridWetData is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation. A copy of GNU Lesser General Public License pertaining to GridWetData is provided in
#    LICENSE/lgpl.txt, referring to the GNU General Public License provided in LICENSE/gpl.txt
#
#    GridWetData is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    Further information about GNU licenses are given at http://www.gnu.org/licenses
#
#    perth3_src.f is written by R. Ray, which have all rights to perth3_src.f, the source code is was distributed by NASA
#    perth3_src.f is part of the DTU10 model, distributed by DTU SPACE. 
#    perth3_src.f is not part of GridWetData, but only co-distributed with GridWetData for user convenience
#    tidal_constituents.dat is a sub set of the glabal DTU10 data set of tidal constituents
# <br>
# ---------------------------------------------------------------------------------------------------------------
# <br>
#
from grid_data import *   # triggers import of grids and support libraries
import derived_layers
import grids           # offer encapsulated environment for import GridWetData.grids
import grid_data       # offer encapsulated environment for import GridWetData.grid_data
#
