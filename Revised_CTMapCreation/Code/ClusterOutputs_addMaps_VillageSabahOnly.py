#! /usr/bin/env python3

# Quick script to add together maps after they have been processed from the cluster

# NOTE! This version adds together only those GFlow output maps that were 
#    created from source points from Sabah and those labelled "Village" or "Other"
#    from Sarawak

#!/usr/bin/env python

# Import libraries
import os
import gdal
import re
import rasterio as rio
import numpy as np
import glob
import shutil
from datetime import datetime
from collections import defaultdict

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
basedir = os.path.abspath('/home/mairin/Documents/GradSchool/Research/CircuitTheory_Borneo/PRSB_Revision2/Revised_CTMapCreation/')

# Output maps file names:
#   Calculate these from maps created in August (where rivers are passable)
cum_outmap_vs_pr = os.path.join(out_dir, date + '_vs_nodes_passrivers_cumulative.asc')
cum_outmap_vso_pr = os.path.join(out_dir, date + '_vso_nodes_passrivers_cumulative.asc')

#   Calculate these from July maps (where rivers are impassable)
cum_outmap_vs_pr = os.path.join(out_dir, date + '_vs_nodes_impassrivers_cumulative.asc')
cum_outmap_vso_pr = os.path.join(out_dir, date + '_vso_nodes_impassrivers_cumulative.asc')

# Output cumulative map location
date = datetime.now().strftime('%d.%b.%Y') 
cum_outmap = os.path.join(out_dir, date + '_rivers_cumulative.asc')

### Node files
# Original: Impassable rivers:
out_dir_j = os.path.join(basedir, 'ClusterFlexOutputs/TmpOutputs_Jul30_2019')
tsv_dir_j = os.path.join(basedir, 'SourceSinks/17.Jul.2019_Nodes/NodesTSV/')
txt_dir_j = os.path.join(basedir, 'SourceSinks/17.Jul.2019_Nodes/NodesTXT/')

# Original: Passable rivers:
out_dir_a = os.path.join(basedir, 'ClusterFlexOutputs/TmpOutputs_Aug05_2019/')
tsv_dir_a = os.path.join(basedir, 'SourceSinks/26.Jul.2019_Nodes/NodesTSV/')
txt_dir_a = os.path.join(basedir, 'SourceSinks/26.Jul.2019_Nodes/NodesTXT/')

# Subset: Note! These were created from a resistance map with impassable rivers
#   This is important when running the follow-up CT simulations, not for matching
txt_dir_vs = os.path.join(basedir, 'SourceSinks_VillageSabah/09.Sep.2019_Nodes/NodesTXT/')
tsv_dir_vs = os.path.join(basedir, 'SourceSinks_VillageSabah/09.Sep.2019_Nodes/NodesTSV/')

txt_dir_vso = os.path.join(basedir, 'SourceSinks_VillageSabahOther/10.Sep.2019_Nodes/NodesTXT/')
tsv_dir_vso = os.path.join(basedir, 'SourceSinks_VillageSabahOther/10.Sep.2019_Nodes/NodesTSV/')

txt_files_vs = glob.glob(txt_dir_vs + "*.txt") # village sabah = 205 files
txt_vs = [os.path.basename(x) for x in txt_files_vs]
tsv_files_vs = glob.glob(tsv_dir_vs + "*.tsv") # village sabah = 205 files
tsv_vs = [os.path.basename(x) for x in tsv_files_vs]

txt_files_vso = glob.glob(txt_dir_vso + "*.txt") # village sabah other = 251 files
txt_vso = [os.path.basename(x) for x in txt_files_vso]
tsv_files_vso = glob.glob(tsv_dir_vso + "*.tsv") # village sabah other = 251 files
tsv_vso = [os.path.basename(x) for x in tsv_files_vso]

txt_files_a = glob.glob(txt_dir_a + "*.txt") # village sabah = 205 files
txt_a = [os.path.basename(x) for x in txt_files_a]
tsv_files_a = glob.glob(tsv_dir_a + "*.tsv") # village sabah = 205 files
tsv_a = [os.path.basename(x) for x in tsv_files_a]

txt_files_j = glob.glob(txt_dir_j + "*.txt") # village sabah = 205 files
txt_j = [os.path.basename(x) for x in txt_files_j]
tsv_files_j = glob.glob(tsv_dir_j + "*.tsv") # village sabah = 205 files
tsv_j = [os.path.basename(x) for x in tsv_files_j]

#########################################
# Possible options: 
#   Match up the SOURCE nodes in the original node folders, and identify which of these
#       are also identified in the limited node files
#   If the files are entirely matching, this map can be included
#   Move the map file into a folder?
#   Then, with the remaining nodes in the subfolders, run CS on these

usable_a = [] # array to keep track of which nodes have been properly calculated in GFlow
unusable_a = []

# 1. Figure out how many source nodes are in each of the files

# August results first (passable rivers)
for node in tsv_a:
    # Update here
    tsv_files_vs = glob.glob(tsv_dir_vs + "*.tsv") # village sabah = 205 files
    tsv_vs = [os.path.basename(x) for x in tsv_files_vs]
    tsv_files_vso = glob.glob(tsv_dir_vso + "*.tsv") # village sabah = 205 files
    tsv_vso = [os.path.basename(x) for x in tsv_files_vso]
    # node is the original node file, find corresponding files in the filtered set
    if node in tsv_vs:
        with open(tsv_dir_a + node) as file:
            head = [next(file) for x in range(1)]
            nnode_a = int(re.split("\t", head[0])[1].replace("\n", ""))-1
        with open(tsv_dir_vs + node) as file:
            head = [next(file) for x in range(1)]
            nnode_vs = int(re.split("\t", head[0])[1].replace("\n", ""))-1
        if nnode_a == nnode_vs:
            ### Check that the source nodes are the same
            with open(txt_dir_a + os.path.splitext(node)[0] + '.txt') as file:
                source_nodes_a = [next(file) for x in range(nnode_a)]
            with open(txt_dir_vs + os.path.splitext(node)[0] + '.txt') as file:
                source_nodes_vs = [next(file) for x in range(nnode_vs)]
            if source_nodes_a == source_nodes_vs:            
                usable_a.append(node)
                shutil.copy(tsv_dir_vs + node, os.path.join(tsv_dir_vs, "Calculated_a", node))
                shutil.copy(txt_dir_vs + os.path.splitext(node)[0] + '.txt', os.path.join(txt_dir_vs, "Calculated_a", os.path.splitext(node)[0]))
            else:
                unusable_a.append(node)
                shutil.copy(tsv_dir_vs + node, os.path.join(tsv_dir_vs, "NotCalculated_a", node))
                shutil.copy(txt_dir_vs + os.path.splitext(node)[0] + '.txt', os.path.join(txt_dir_vs, "NotCalculated_a", os.path.splitext(node)[0] + '.txt'))
        else:
            unusable_a.append(node)
            shutil.copy(tsv_dir_vs + node, os.path.join(tsv_dir_vs, "NotCalculated_a", node))
            shutil.copy(txt_dir_vs + os.path.splitext(node)[0] + '.txt', os.path.join(txt_dir_vs, "NotCalculated_a", os.path.splitext(node)[0] + '.txt'))
    else:
        unusable_a.append(node)
    if node in tsv_vso:
        with open(tsv_dir_a + node) as file:
            head = [next(file) for x in range(1)]
            nnode_a = int(re.split("\t", head[0])[1].replace("\n", ""))-1
        with open(tsv_dir_vso + node) as file:
            head = [next(file) for x in range(1)]
            nnode_vso = int(re.split("\t", head[0])[1].replace("\n", ""))-1
        if nnode_a == nnode_vso:
            ### Check that the source nodes are the same
            with open(txt_dir_a + os.path.splitext(node)[0] + '.txt') as file:
                source_nodes_a = [next(file) for x in range(nnode_a)]
            with open(txt_dir_vso + os.path.splitext(node)[0] + '.txt') as file:
                source_nodes_vso = [next(file) for x in range(nnode_vso)]
            if source_nodes_a == source_nodes_vso:            
                usable_a.append(node)
                shutil.copy(tsv_dir_vso + node, os.path.join(tsv_dir_vso, "Calculated_a", node))
                shutil.copy(txt_dir_vso + os.path.splitext(node)[0] + '.txt', os.path.join(txt_dir_vso, "Calculated_a", os.path.splitext(node)[0]))
            else:
                unusable_a.append(node)
                shutil.copy(tsv_dir_vso + node, os.path.join(tsv_dir_vso, "NotCalculated_a", node))
                shutil.copy(txt_dir_vso + os.path.splitext(node)[0] + '.txt', os.path.join(txt_dir_vso, "NotCalculated_a", os.path.splitext(node)[0] + '.txt'))
        else:
            unusable_a.append(node)
            shutil.copy(tsv_dir_vso + node, os.path.join(tsv_dir_vso, "NotCalculated_a", node))
            shutil.copy(txt_dir_vso + os.path.splitext(node)[0] + '.txt', os.path.join(txt_dir_vso, "NotCalculated_a", os.path.splitext(node)[0] + '.txt'))
    else:
        unusable_a.append(node)

# July results
usable_j = [] # array to keep track of which nodes have been properly calculated in GFlow
unusable_j = []

for node in tsv_j:
    # Update here
    tsv_files_vs = glob.glob(tsv_dir_vs + "*.tsv") # village sabah = 205 files
    tsv_vs = [os.path.basename(x) for x in tsv_files_vs]
    # node is the original node file, find corresponding files in the filtered set
    if node in tsv_vs:
        with open(tsv_dir_j + node) as file:
            head = [next(file) for x in range(1)]
            nnode_j = int(re.split("\t", head[0])[1].replace("\n", ""))-1
        with open(tsv_dir_vs + node) as file:
            head = [next(file) for x in range(1)]
            nnode_vs = int(re.split("\t", head[0])[1].replace("\n", ""))-1
        if nnode_j == nnode_vs:
            ### Check that the source nodes are the same
            with open(txt_dir_j + os.path.splitext(node)[0] + '.txt') as file:
                source_nodes_j = [next(file) for x in range(nnode_j)]
            with open(txt_dir_vs + os.path.splitext(node)[0] + '.txt') as file:
                source_nodes_vs = [next(file) for x in range(nnode_vs)]
            if source_nodes_j == source_nodes_vs:            
                usable_j.append(node)
                shutil.copy(tsv_dir_vs + node, os.path.join(tsv_dir_vs, "Calculated_j", node))
                shutil.copy(txt_dir_vs + os.path.splitext(node)[0] + '.txt', os.path.join(txt_dir_vs, "Calculated_j", os.path.splitext(node)[0]))
            else:
                unusable_j.append(node)
                shutil.copy(tsv_dir_vs + node, os.path.join(tsv_dir_vs, "NotCalculated_j", node))
                shutil.copy(txt_dir_vs + os.path.splitext(node)[0] + '.txt', os.path.join(txt_dir_vs, "NotCalculated_j", os.path.splitext(node)[0] + '.txt'))
        else:
            unusable_j.append(node)
            shutil.copy(tsv_dir_vs + node, os.path.join(tsv_dir_vs, "NotCalculated_j", node))
            shutil.copy(txt_dir_vs + os.path.splitext(node)[0] + '.txt', os.path.join(txt_dir_vs, "NotCalculated_j", os.path.splitext(node)[0] + '.txt'))
    else:
        unusable_j.append(node)
    unusable_j_vs = unusable_j
    usable_j_vs = usable_j
    unusable_j = []
    usable_j = []
    if node in tsv_vso:
        with open(tsv_dir_j + node) as file:
            head = [next(file) for x in range(1)]
            nnode_j = int(re.split("\t", head[0])[1].replace("\n", ""))-1
        with open(tsv_dir_vso + node) as file:
            head = [next(file) for x in range(1)]
            nnode_vso = int(re.split("\t", head[0])[1].replace("\n", ""))-1
        if nnode_j == nnode_vso:
            ### Check that the source nodes are the same
            with open(txt_dir_j + os.path.splitext(node)[0] + '.txt') as file:
                source_nodes_j = [next(file) for x in range(nnode_j)]
            with open(txt_dir_vso + os.path.splitext(node)[0] + '.txt') as file:
                source_nodes_vso = [next(file) for x in range(nnode_vso)]
            if source_nodes_j == source_nodes_vso:            
                usable_j.append(node)
                shutil.copy(tsv_dir_vso + node, os.path.join(tsv_dir_vso, "Calculated_j", node))
                shutil.copy(txt_dir_vso + os.path.splitext(node)[0] + '.txt', os.path.join(txt_dir_vso, "Calculated_j", os.path.splitext(node)[0]))
            else:
                unusable_j.append(node)
                shutil.copy(tsv_dir_vso + node, os.path.join(tsv_dir_vso, "NotCalculated_j", node))
                shutil.copy(txt_dir_vso + os.path.splitext(node)[0] + '.txt', os.path.join(txt_dir_vso, "NotCalculated_j", os.path.splitext(node)[0] + '.txt'))
        else:
            unusable_j.append(node)
            shutil.copy(tsv_dir_vso + node, os.path.join(tsv_dir_vso, "NotCalculated_j", node))
            shutil.copy(txt_dir_vso + os.path.splitext(node)[0] + '.txt', os.path.join(txt_dir_vso, "NotCalculated_j", os.path.splitext(node)[0] + '.txt'))
    else:
        unusable_j.append(node)

# Then do calculations I guess
