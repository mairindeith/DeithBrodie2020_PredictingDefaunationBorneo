#! /usr/bin/env python

# -----------------------
# Supplementary Material for Deith and Brodie 2020; “Predicting defaunation – accurately mapping bushmeat hunting pressure over large areas”
#    doi: 10.1098/rspb.2019-2677
#------------------------
# Code to create a resistance surface from standardized files
#     for the Malaysian states of Borneo
#
# This script uses maps created with 'MSYBorneo_ResistanceMapLayers_Preparation.py'
# Input files should be standard resolution, extent, and projection.
# Created by Mairin Deith on Nov12 2017
# Last edited for revised mapping effort on Jul26 2019
#------------------------
# Hierarcy of resistance values:
#   Walking speed - modified by land cover type (incl. plantations) and slope
#       - presumed maximum at ideal slope is 6kph
#   Roadways - logging roads have a max speed of 20kph, major roads 60 kph
#   Rivers - passable rivers (Strahler order <3), navigable rivers (Strahler order > 5 is navigable at 20kph)

#!/usr/bin/env python

# Import libraries
import numpy as np
import gdal # for opening files as np arrays
import gc
import os, sys
import math
from pathlib import Path
import codecs
import rasterio as rio
from datetime import datetime

# Global parameters

output_path = Path('/home/mairin/Documents/GradSchool/Research/CircuitTheory_Borneo/PRSB_Revision2/Revised_CTMap/ResistanceMap')
input_path = output_path / 'ModifiedSpatialData'
date = datetime.now().strftime('%d.%b.%Y')

# Input file locations

fc_loc = input_path / 'REGIONBorneo_FCDefDeg_1973to2010_CIFOR_clipped_rs.tif'
plantation_loc = input_path / 'GaveauPlantationClip_raster.tif'
slope_loc = input_path / 'RDSlope_pc.tif'
strahler_loc = input_path / 'StrahlerOrder_rs.tif'
road_loc = input_path / 'rasterizedRoads.tif'

### LANDUSE

# Use the first dataset, in this case forest cover, as template data for raster creation
#   tmp_shape indicates the ideal shape for the output rasters (based on fc_arr)
#   tmp_meta contains information needed to write the geospatial .tif file

with rio.open(str(fc_loc), 'r') as ds:
    tmp_meta = ds.profile
    fc_arr = ds.read()  # read all raster values

tmp_shape = fc_arr.shape

# Modify meta data

tmp_meta['count'] = 1
tmp_meta['dtype'] = "float64"

# Create empty array to populate with speed values
speed_arr = np.empty(tmp_shape)
speed_arr[:] = np.nan

fc_arr[fc_arr == -9999] = np.nan

with rio.open(str(slope_loc), 'r') as ds:
    slope_arr = ds.read()  # read all raster values

if np.shape(slope_arr) != tmp_shape:
    quit("Error! Could not create resistance map.\n Reason: slope_arr and fc_arr do not have matching dimensions.")

slope_arr[slope_arr == -9999] = np.nan
# Convert percent to dh/dx
slope_arr = slope_arr/100

# Calculate walking speed from Tobler's function based on slope (in dh/dx)

tobler_modifier = np.where(np.isnan(slope_arr), 1.0, np.exp(-3.5*(slope_arr + 0.05)))

# New speed on land sources from wonderful MAP project
# https://www-nature-com.ezproxy.library.ubc.ca/articles/nature25181
# A global map of travel time to cities to assess inequalities in accessibility in 2015, Weiss et al. 2018

# Broadleafed evergreen

speed_arr[fc_arr == 1] = 1.62 * tobler_modifier[fc_arr == 1]

# Logged forest

speed_arr[fc_arr == 2] = 1.0 * tobler_modifier[fc_arr == 2]

# Cropland (will be overwritten by plantations in much of the area)

speed_arr[fc_arr == 3] = 2.5 * tobler_modifier[fc_arr == 3]

# Urban areas?

speed_arr[fc_arr == 4] = 5.0 * tobler_modifier[fc_arr == 4]

# Clear up memory

del fc_arr
del slope_arr
del tobler_modifier
gc.collect()

### PLANTATIONS

with rio.open(str(plantation_loc), 'r') as ds:
    pl_arr = ds.read()  # read all raster values

if pl_arr.shape != tmp_shape:
    quit("Error! Could not create resistance map.\n Reason: pl_arr and does not match template dimensions.")

# Try 10 kph
speed_arr[pl_arr == 1] = 10.0

del pl_arr
gc.collect()

# RIVERS

with rio.open(str(strahler_loc), 'r') as ds:
    st_arr = ds.read()

if st_arr.shape != tmp_shape:
    quit("Error! Could not create resistance map.\n Reason: slope_arr and does not match template dimensions.")

# Strahler order: 1 and 2 - walking speed
#                 3 and 4 - impassible - maybe this is too much...
#                 5+      - 20 kph
# Modification: intermediate rivers
speed_arr[st_arr == 4] = 10.0
speed_arr[st_arr == 5] = 20.0

del st_arr
gc.collect()

# ROADS

with rio.open(str(road_loc), 'r') as ds:
    rd_arr = ds.read()

if rd_arr.shape != tmp_shape:
    quit("Error! Could not create resistance map.\n Reason: rd_arr does not match template dimensions.")

# Roads: 1 - main roads
#        2 - logging
speed_arr[rd_arr == 1] = 60.0
speed_arr[rd_arr == 2] = 20.0

del rd_arr
gc.collect()

# GENERATE RESISTANCE SURFACE

# Save speed file

with rio.open(str(output_path / str('Speed_kph_' + date + '.tif')), 'w', **tmp_meta) as ds:
    ds.write(np.squeeze(speed_arr), 1) # squeeze removes the empty 1st dimension for writing

# Save resistance (in hpk)

with rio.open(str(output_path / str('Resistance_hpk_' + date + '.tif')), 'w', **tmp_meta) as ds:
    ds.write(np.squeeze(np.reciprocal(speed_arr)), 1) # squeeze removes the empty 1st dimension for writing
