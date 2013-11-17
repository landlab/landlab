######### radiation_driver.py ###########
##
## This is an example for implementation
## of radiation.py
##
##
## Written by Sai Nudurupati & Erkan Istanbulluoglu, 2013
#############################################

# Import Required Classes

import os
from landlab.components.radiation.radiation import *
from landlab.io import read_esri_ascii

# Input DEM
_DEFAULT_INPUT_FILE_2 = os.path.join(os.path.dirname(__file__),
                                 'HugoSite.asc')

# Create & Initialize Radiation Object
RD = Radiation()
RD.initialize()

# Importing Grid and Elevations from DEM
(RMG,data) = read_esri_ascii(_DEFAULT_INPUT_FILE_2)

# Transfer Elevation Data
RD._Z = RMG.zeros(centering='node')
RD._Z = data  

# Temporary Adjustment of 'No Data' values
Max = max(RD._Z)
RD._Z[RD._Z == -9999] = Max     

# Update Radiation
RD.update( RMG )

# Plot Radiation Data
RD.plott( RMG )
show()
