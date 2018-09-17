## \mainpage GridWetData
# GridWetData.CMEMS.BAL is a package to access data output from the BAL circulation model(s) output
# grid definitions:    GridWetData/CMEMS/BAL/BAL_grids
# grid data:           GridWetData/CMEMS/BAL/BAL_grid_data
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
# <br>
# ---------------------------------------------------------------------------------------------------------------
# <br>
#
from BAL_grid_data import *   # triggers import of data_manager, grids and support libraries
import BAL_grids              # offer encapsulated environment for import GridWetData.CMEMS.BAL.BAL_grids
import BAL_grid_data          # offer encapsulated environment for import GridWetData.CMEMS.BAL.BAL_grid_data
