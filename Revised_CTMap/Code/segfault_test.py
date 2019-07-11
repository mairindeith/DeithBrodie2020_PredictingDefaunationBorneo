#------------------------
# Code to create a resistance surface from standardized files
#     for the Malaysian states of Borneo
#
# This script uses maps created with 'MSYBorneo_ResistanceMapLayers_Preparation.py'
# Created by Mairin Deith on Nov12 2017
# Last edited for revised mapping effort on Jul4 2019
#------------------------
# Hierarcy of resistance values:
#   Walking speed - modified by land cover type (incl. plantations) and slope
#   Roadways - logging roads have a max speed of 20kph, major roads 60 kph
#   Rivers - passable rivers (Strahler order <3), navigable rivers (Strahler order > 5 is navigable at 20kph)

#!/usr/bin/env python

# Import libraries
import numpy as np

array_a = np.array([
    [0,0,0],
    [0,1,1],
    [0,2,2]
    ])
    
array_b = np.array([
    [0,1,2],
    [0,1,0.5],
    [0,1,1.25]
    ])
    
array_c = np.empty(array_b.shape)
    
array_c[:] = np.nan
    
array_c[np.where(np.isnan(array_c))] = array_a[np.where(np.isnan(array_c))]*array_b[np.where(np.isnan(array_c))]
