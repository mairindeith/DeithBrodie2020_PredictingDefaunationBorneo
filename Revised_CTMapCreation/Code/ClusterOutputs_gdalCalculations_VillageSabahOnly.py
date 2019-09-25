#!/usr/bin/env python

# What this code needs to do:
# 1. Iterate through the INPUT files (i.e. the .tsv and .txt files)
# 2. Use the name of the file to find the associated OUTPUT file in:
#       ../../SourceSinks/17.Jul.2019_Nodes/
#   3. Move the .asc file to the < Processed > folder 
# .asc file structure: 
#       output_temp_{iter}_{population size}
#   If there are multiple source files, I don't know which was _0/_1/_2, etc.
#       Solution: add up all .asc files, source counts, iteration counts and
#           Use this as the divisor.

# Import libraries

import os
import shutil
import gdal
import re
import glob
import rasterio as rio
import numpy as np
from pathlib import Path
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

### CHANGE THESE PARAMETERS 
# File paths/gdalcalc paths

basedir = os.path.abspath('/home/mairin/Documents/GradSchool/Research/CircuitTheory_Borneo/PRSB_Revision2/Revised_CTMapCreation/')
basefname = 'ClusterOutput_Summation'
map_dir = os.path.join(basedir, 'GFlowOutputs', 'ClusterFlexOutputs')

### Calculated maps (ie. those not sent to the cluster)
# Passable rivers (August outputs)
maps_pass_dir = os.path.join(basedir, 'GFlowOutputs/ClusterFlexOutputs/AugCluster_GazetteerPassableRivers_AllNodes/Preprocessing')
# Impassable rivers (July outputs)
maps_impass_dir = os.path.join(basedir, 'GFlowOutputs/ClusterFlexOutputs/JulCluster_GazetteerImpassableRivers_AllNodes/Preprocessing')

# Output cumulative map location
date = datetime.now().strftime('%d.%b.%Y')
out_dir = os.path.join(basedir, 'GFlowOutputs/ClusterFlexOutputs/')

complex_calcs = ['VillageSabah', 'VillageSabahOther']
rivers = ['passable','impassable']
for cc in complex_calcs:
    print("CC = %s" %(cc))
    # Set paths:
    for r in rivers:
        print("R = %s" %(r))
        # Directories where processed .ascs live:
        if r=='passable':
            old_maps_dir = maps_pass_dir
        elif r=='impassable':
            old_maps_dir = maps_impass_dir
        new_maps_dir = os.path.join(map_dir, 'SepTmpOutputs_Gazetteer'+cc, r + 
            'Rivers')
        # List of already calculated .ascs
        asc_list = glob.glob(old_maps_dir+"/output_temp_*")
        
        # Old and new source directories
        #   calc_dir = calculated before September 2019 (i.e. including all sources)
        #   uncalc_dir = calculated after Sept 2019
        source_dir = os.path.join(basedir, 'SourceSinks_'+cc) 
        # Calculated nodes:
        calc_dir = os.path.join(source_dir, 'NodesTSV','Calculated_' + r)
        calc_n = glob.glob(calc_dir + '/*.tsv')
        # Setup outfile and check its existance 
        output_fname = basefname + '_' + cc + '_' + r + '.asc'
        cum_outmap = os.path.join(map_dir, 'SepTmpOutputs_Gazetteer'+cc, r + 
            'Rivers', output_fname)
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
        ### CREATE EMPTY MAP OR READ EXISTING
        if exists == False:
            with rio.open(os.path.join(old_maps_dir,asc_list[0]), 'r') as ds:
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
            with rio.open(os.path.join(old_maps_dir,asc_list[0]), 'r') as ds:
                tmp_meta = ds.profile
        
        ### Iterate through all of the already calculated sources
        print("...iterating through calculated maps")
        calc_sources_pops = []
        for n in calc_n:     # for each source/sink file:
            # Isolate just the population density
            ptmp = float(re.split("_", os.path.basename(n))[1])
            calc_sources_pops.append(ptmp)
        calc_sources_singles = sorted(list_singles(calc_sources_pops))
        calc_sources_dups = sorted(list_duplicates(calc_sources_pops))
        
        save_counter = 0 # save every 5 steps 
        singles_counter = 0
        total = len(calc_sources_singles)
        
        for s in calc_sources_singles:
            ptmp = s[0]
            multtmp = np.log10(float(ptmp)/78.54)*-0.054357+0.179789
            popmult = multtmp*float(ptmp)
            maptmp = glob.glob(old_maps_dir+"/*%s*.asc" %(ptmp))
            if len(maptmp)!=1:
                print("Warning! Mismatching maps for pop density %s \n...skipping..." %(ptmp))
                for f in range(len(maptmp)):
                    Path(os.path.join(old_maps_dir, cc+"Problems", os.path.basename(maptmp[f]))).touch()
                singles_counter += 1
                continue
            singles_counter += 1
            print("\n...Processing map %s of %s...\n" %(singles_counter, total))
            with open(calc_n[s[1][0]]) as file:
                head = [next(file) for x in range(1)]
            nsources = int(re.split("\t", head[0])[1].replace("\n", ""))-1
            niter = int(re.split("_", os.path.basename(maptmp[0]))[2])  # number of iterations
            print("......Opening file: %s" %(os.path.basename(maptmp[0])))
            with rio.open(maptmp[0], 'r') as ds:
                tmp = ds.read()
            # print(np.unique(tmp))
            tmp[tmp == -9999] = 0
            intmap = intmap + ((tmp * nsources * popmult) / niter)
            Path(os.path.join(old_maps_dir, cc+"Processed", os.path.basename(maptmp[0]))).touch()
            save_counter += 1
            if save_counter == 10:
                with rio.open(cum_outmap, 'w', **tmp_meta) as ds:
                    ds.write(intmap.astype(np.float32))
                save_counter = 0
        if total!=0:
            with rio.open(cum_outmap, 'w', **tmp_meta) as ds:
                ds.write(intmap.astype(np.float32))
                
        # Duplicates for calculated maps
        save_counter = 0
        total = len(calc_sources_dups)
        dup_counter = 0
        for s in calc_sources_dups:
            dup_counter += 1
            print("\n...Processing map %s of %s...\n" %(dup_counter, total))
            ptmp = s[0]
            multtmp = np.log10(float(ptmp)/78.54)*-0.054357+0.179789
            popmult = multtmp*float(ptmp)
            maptmp = glob.glob(old_maps_dir+"/*%s*.asc" %(ptmp))
            if len(maptmp)!= len(s[1]):
                print("Warning! Mismatching maps for pop density %s \n...skipping..." %(ptmp))
                for f in range(len(maptmp)):
                    Path(os.path.join(old_maps_dir, cc+"Problems", os.path.basename(maptmp[f]))).touch()
                dup_counter += 1
                continue
            dup_counter += 1
            tmp_array = np.zeros(shape = tmp.shape)
            niter = 0
            nsources = 0
            for f in range(len(s[1])):
                with open(calc_n[s[1][f]]) as file:
                    head = [next(file) for x in range(1)]
                nsources += int(re.split("\t", head[0])[1].replace("\n", ""))-1
                niter += int(re.split("_", os.path.basename(maptmp[f]))[2]) 
                print("......Opening file: %s" %(os.path.basename(maptmp[f])))
                with rio.open(maptmp[f], 'r') as ds:
                    tmp = ds.read()
                tmp[tmp == -9999] = 0
                tmp_array = tmp_array + tmp
                Path(os.path.join(old_maps_dir, cc+"Processed", os.path.basename(maptmp[f]))).touch()
            intmap = intmap + ((tmp_array * nsources * popmult) / niter)
            save_counter += len(s[1])
            if save_counter == 10:
                with rio.open(cum_outmap, 'w', **tmp_meta) as ds:
                    ds.write(intmap.astype(np.float32))
                save_counter = 0
        with rio.open(cum_outmap, 'w', **tmp_meta) as ds:
            ds.write(intmap.astype(np.float32))
            
        # Now, for each map type, perform simpler calculations 
        #   over the newly calculated files
        
print("Now calculating newly generated GFlow maps")

complex_calcs = ['VillageSabahOther'] #, 'VillageSabah']
rivers = ['passable'] # 'impassable', 
for cc in complex_calcs:
    print("CC = %s" %(cc))
    # Set paths:
    for r in rivers:
        print("R = %s" %(r))
        # Save to a new file, just to be safe
        output_fname = basefname + '_' + cc + '_' + r + '_newGflowMaps.asc'
        cum_outmap = os.path.join(map_dir, 'SepTmpOutputs_Gazetteer' + cc, r + 
            'Rivers', output_fname)
        new_maps_dir = os.path.join(map_dir, 'SepTmpOutputs_Gazetteer' + cc, r + 
            'Rivers')
        # asc_list = glob.glob(new_maps_dir+cc+"_"+r+"RiverOutput"+"/output_temp_*")    
        asc_list = glob.glob(new_maps_dir + "/output_temp_*")    
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
        ### CREATE EMPTY MAP OR READ EXISTING
        if exists == False:
            with rio.open(os.path.join(new_maps_dir,asc_list[0]), 'r') as ds:
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
            with rio.open(os.path.join(new_maps_dir,asc_list[0]), 'r') as ds:
                tmp_meta = ds.profile
        # Source/sink files
        source_dir = os.path.join(basedir, 'SourceSinks_'+cc) 
        # Uncalculated nodes:
        uncalc_dir = os.path.join(source_dir, 'NodesTSV','NotCalculated_' + r)
        ### uncalc_n = glob.glob(uncalc_dir + '/*.tsv')
        iterpop = [a.replace('output_temp_', '').replace('nodes_', '').replace('.asc','') for a in asc_list]
        niter = []
        p = []
        p_mult = []
        p_mult_tmp = []
        for a in iterpop:
            ntmp, ptmp = re.split("_", os.path.basename(a))[0:2]
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
            src_f = glob.glob(uncalc_dir+"/*%s*.tsv" %(str(ptmp)+"_"))[0]
            with open(src_f) as file:
                head = [next(file) for x in range(1)]
            # Number of sources is just the first # of the second column minus 1
            nsource = int(re.split("\t", head[0])[1].replace("\n", ""))-1
            asc_tmp = os.path.join(new_maps_dir, asc_list[int(s[1][0])])
            print("......Opening file: %s" %(os.path.basename(asc_tmp)))
            with rio.open(asc_tmp, 'r') as ds:
                tmp = ds.read()
            tmp[tmp == -9999] = 0
            intmap = intmap + ((tmp * nsource * pmtmp) / ntmp)
            save_counter += 1
            Path(os.path.join(new_maps_dir, "Processed", os.path.basename(asc_tmp[0]))).touch()
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
                src_f = glob.glob(uncalc_dir+"/*%s*" %(str(ptmp)+"_"))[f]
                with open(src_f) as file:
                    head = [next(file) for x in range(1)]
                nsource = int(re.split("\t", head[0])[1].replace("\n", ""))-1    
                asc_tmp = os.path.join(new_maps_dir, asc_list[int(d[1][f])])
                print("......Opening file: %s" %(os.path.basename(asc_tmp)))
                with rio.open(asc_tmp, 'r') as ds:
                    tmp = ds.read()
                tmp[tmp == -9999] = 0
                intmap = intmap + ((tmp * nsource * pmtmp) / ntmp)
                save_counter += 1
                Path(os.path.join(new_maps_dir, cc+"Processed", os.path.basename(asc_tmp[0]))).touch()
            if save_counter >= 10:
                with rio.open(cum_outmap, 'w', **tmp_meta) as ds:
                    ds.write(intmap)
                save_counter = 0

with rio.open(cum_outmap, 'w', **tmp_meta) as ds:
    ds.write(intmap)
