# What this code needs to do:

# 1. Iterate through GFlow output files
#   2. Use the name of the file to find the associated file in:
#       ../../SourceSinks/17.Jul.2019_Nodes/
#   3. Move the .asc file to the < Processed > folder 

# .asc file structure: 
#       output_temp_{iter}_{population size}
#   If there are multiple source files, I don't know which was _0/_1/_2, etc.
#       Solution: add up all .asc files, source counts, iteration counts and
#           Use this as the divisor. 

#!/usr/bin/env python

# Import libraries
import os
import gdal
import re
import glob
import rasterio as rio
import numpy as np
from datetime import datetime
from collections import defaultdict

### HELPER FUNCTIONS
# From https://stackoverflow.com/questions/5419204/index-of-duplicates-items-in-a-python-list
# This function passes through a list and keeps a list of locations seen for each item, 
#   and returns items seen more than once

def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() 
                            if len(locs)>1)

def list_singles(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() 
                            if len(locs)==1)

### GLOBAL INFO
# File paths/gdalcalc paths

basedir = os.path.abspath('/home/mairin/Documents/GradSchool/Research/CircuitTheory_Borneo/PRSB_Revision2/Revised_CTMap/')
out_dir = os.path.join(basedir, 'ClusterFlexOutputs/TmpOutputs_Jul30_2019')
asc_dir = os.path.join(out_dir, 'ASCs')
src_dir = os.path.join(basedir, 'SourceSinks/17.Jul.2019_Nodes/NodesTSV')

# Output cumulative map location
date = datetime.now().strftime('%d.%b.%Y') 
cum_outmap = os.path.join(out_dir, date + '_cumulative.asc')

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

### WORK ON FILES TO GENERATE CALCULATION INFORMATION

asc_list = glob.glob(asc_dir+"/output_temp_*")

# Modify file names to isolate just numbers
iterpop = [a.replace('output_temp_', '').replace('.asc','') for a in asc_list]

# Separate the number of iterations from population size (for lookups)

niter = []
p = []
p_mult = []
p_mult_tmp = []

for a in iterpop:
    ntmp, ptmp = re.split("_", os.path.basename(a))
    niter.append(int(ntmp))
    p.append(float(ptmp))
    multtmp = np.log10(float(ptmp)/78.54)*-0.054357+0.179789
    if multtmp < 0:
        multtmp = 0.0001
    p_mult_tmp.append(multtmp)
    p_mult.append(multtmp*float(ptmp))
    # Calculates the expected number of hunters based on mini lit review
    #   Each population density was calculated within a 5km buffer; so divide
    #   by 78.54 (approximate area of a 5km-radius circle) for density per km2

# Find singles/duplicates:
single_src = sorted(list_singles(p))
dup_src = sorted(list_duplicates(p))

### CREATE EMPTY MAP OR READ EXISTING
# Use same dimensions/transformation as asc_list[0]
if exists == False:
    with rio.open(os.path.join(asc_dir,asc_list[0]), 'r') as ds:
        tmp_meta = ds.profile
        tmp_shape = ds.shape  # read in height and width
    intmap = np.empty(tmp_shape, dtype=np.float32)
    print("Creating blank output map: %s" %(os.path.basename(cum_outmap)))
    with rio.open(cum_outmap, 'w', **tmp_meta) as ds:
        ds.write_band(1, intmap)

elif exists == True:
    print("Reading existing cumulative map %s" %(os.path.basename(cum_outmap)))
    with rio.open(cum_outmap, 'r') as ds:
        intmap = ds.read()
    with rio.open(os.path.join(asc_dir,asc_list[0]), 'r') as ds:
        tmp_meta = ds.profile

### BEGIN CALCULATIONS 

print("Calculating single-source file maps...")

save_counter = 0 # save every 5 steps 
singles_counter = 0
total = len(single_src)

for s in single_src:
    singles_counter += 1
    print("\n...Processing map %s of %s...\n" %(singles_counter, total))
    ntmp = niter[int(s[1][0])]
    pmtmp = p_mult[int(s[1][0])]
    # Find corresponding source TSV file:
    src_f = glob.glob(src_dir+"/*%s*" %(str(ptmp)+"_"))[0]
    with open(src_f) as file:
        head = [next(file) for x in range(1)]
    # Number of sources is just the first # of the second column minus 1
    nsource = int(re.split("\t", head[0])[1].replace("\n", ""))-1
    asc_tmp = os.path.join(asc_dir, asc_list[int(s[1][0])])
    print("......Opening file: %s" %(os.path.basename(asc_tmp)))
    with rio.open(asc_tmp, 'r') as ds:
        tmp = ds.read()
    print(np.unique(tmp))
    tmp[tmp == -9999] = 0
    intmap = intmap + ((tmp * nsource * pmtmp) / ntmp)
    os.rename(asc_tmp, os.path.join(out_dir, "Processed", os.path.basename(asc_tmp)))
    save_counter += 1
    if save_counter == 10:
        with rio.open(cum_outmap, 'w', **tmp_meta) as ds:
            ds.write(intmap)
        save_counter = 0

if total!=0:
    with rio.open(cum_outmap, 'w', **tmp_meta) as ds:
        ds.write(intmap)

total = len(dup_src)
dup_counter = 0

for d in dup_src:
    dup_counter += 1
    print("\n...Processing map %s of %s...\n" %(dup_counter, total))
    # Population at 0th index, location at 1st
    ptmp = d[0]
    pmtmp = p_mult[int(d[1][0])]
    for f in range(len(d[1])):
        ntmp = niter[int(d[1][f])]
        src_f = glob.glob(src_dir+"/*%s*" %(str(ptmp)+"_"))[f]
        with open(src_f) as file:
            head = [next(file) for x in range(1)]
        nsource = int(re.split("\t", head[0])[1].replace("\n", ""))-1    
        asc_tmp = os.path.join(asc_dir, asc_list[int(d[1][f])])
        print("......Opening file: %s" %(os.path.basename(asc_tmp)))
        with rio.open(asc_tmp, 'r') as ds:
            tmp = ds.read()
        tmp[tmp == -9999] = 0
        intmap = intmap + ((tmp * nsource * pmtmp) / ntmp)
        os.rename(asc_tmp, os.path.join(out_dir, "Processed", os.path.basename(asc_tmp)))
    save_counter += 1
    if save_counter == 10:
        with rio.open(cum_outmap, 'w', **tmp_meta) as ds:
            ds.write(intmap)
        save_counter = 0

with rio.open(cum_outmap, 'w', **tmp_meta) as ds:
    ds.write(intmap)
