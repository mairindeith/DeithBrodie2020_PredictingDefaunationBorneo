# Quick script to add together maps after they have been processed from the cluster

#!/usr/bin/env python

# Import libraries
import os
import gdal
import re
import rasterio as rio
import numpy as np
import glob
from datetime import datetime

### GLOBAL INFO
# File paths/gdalcalc paths

basedir = os.path.abspath('/home/mairin/Documents/GradSchool/Research/CircuitTheory_Borneo/PRSB_Revision2/Revised_CTMap/')
out_dir_j = os.path.join(basedir, 'ClusterFlexOutputs/TmpOutputs_Jul30_2019')
out_dir_a = os.path.join(basedir, 'ClusterFlexOutputs/TmpOutputs_Aug05_2019')

# Output cumulative map location
date = datetime.now().strftime('%d.%b.%Y') 
cum_outmap_j = os.path.join(out_dir_j, date + '_final.asc')
cum_outmap_a = os.path.join(out_dir_a, date + '_final.asc')

file_count = 0

# Template info
with rio.open(os.path.join(out_dir_j,"output_temp_17_230.193424225.asc"), 'r') as ds:
    tmp_meta = ds.profile
    tmp_shape = ds.shape # read in height and width

for cum_outmap in [cum_outmap_j, cum_outmap_a]:
    if os.path.exists(cum_outmap): # Don't bother if the file doesn't exist
        exists = True
        delete = str(input("\n\nShould the existing map (%s) be deleted?\n'Y' for yes, 'N' for no > " %(os.path.basename(cum_outmap))))
        delete = delete.upper()
        while delete!='N' and delete!='Y':
            delete = str(input("Sorry, I was expecting either 'Y' or 'N'. \nShould the existing map %s be deleted?\n'Y' for yes, 'N' for no >" %(os.path.basename(cum_outmap))))
            delete = delete.upper()
        if delete=='Y':
            os.remove(cum_outmap)
            exists = False
        if delete=='N':
            exists = True
    else:
        exists = False
    if exists == False:
        intmap = np.empty(tmp_shape, dtype=np.float32)
        print("...Creating blank output map: %s" %(os.path.basename(cum_outmap)))
        with rio.open(cum_outmap, 'w', **tmp_meta) as ds:
            ds.write_band(1, intmap)
    elif exists == True:
        print("...Reading existing cumulative map %s" %(os.path.basename(cum_outmap)))    
        with rio.open(cum_outmap, 'r') as ds:
            intmap = ds.read()
    ### BEGIN CALCULATIONS
    files = glob.glob([out_dir_j, out_dir_a][file_count]+"/*Aug*")
    for f in files:
        print("...Reading %s, adding to %s" %(os.path.basename(f), os.path.basename(cum_outmap)))
        with rio.open(f, 'r') as ds:
            tmp_map = ds.read()
        intmap = intmap + tmp_map
    print("Saving...")
    with rio.open(cum_outmap, 'w', **tmp_meta) as ds:
        ds.write(intmap)
