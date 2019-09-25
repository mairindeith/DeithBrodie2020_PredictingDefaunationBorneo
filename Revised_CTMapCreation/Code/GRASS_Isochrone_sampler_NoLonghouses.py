#!/usr/bin/env python

#------------------------
# Python script to create isochrones and sample sink points from within these isochrones
# This script should be run from within a GRASS session as it uses GRASS libraries
# e.g. After opening a GRASS session in Terminal:
#    $ python isochrone_sampler.py
#
# First created by Mairin Deith in Spring 2018
# Last edited for revised mapping effort on Jul26 2019
#------------------------

from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import vector as v
import grass.script as grass
from osgeo import gdal, ogr
from osgeo import gdal_array
from osgeo import gdalnumeric
import struct
import math
import os, sys
import errno
import numpy as np
import shutil
from datetime import datetime

# Command-line arguments
global overwrite_nodes
global average_res
if "overwrite_nodes" in sys.argv:
    overwrite_nodes = True
else:
    overwrite_nodes = False

if "avg" in sys.argv:
    average_res = True
else:
    average_res = False

# Global parameters
date = datetime.now().strftime('%d.%b.%Y') 
global node_count
global node_total

# Resistance file should be resistance (measured in hours/pixel)
# If needed, the function "RasterConvert" can convert a h/km map into an h/pixel map
infiles = os.path.join('/home/mairin/Documents/GradSchool/Research/CircuitTheory_Borneo/PRSB_Revision2/Revised_CTMapCreation/')

### Jul26 has passable rivers
### resist_infile = os.path.join(infiles,'ResistanceMap','Resistance_hpk_26.Jul.2019.tif')

### Aug3 has impassable rivers
resist_infile = os.path.join(infiles,'ResistanceMap','Resistance_hpk_03.Aug.2019.tif')

# Path to a folder where you want your NodesTSV and NodesTXT folders to be created
# If these folders already exist, this code DELETES THESE and recreates them
#    Change the "SUBFOLDER" to a different value in each instance if you want to
#    run this code in multiple instances

node_path = os.path.join(infiles,'SourceSinks_VillageSabahOther','impassableRivers')
# node_path = os.path.join(infiles,'SourceSinks_SwkNoLongHouses_',str(date+"_Nodes"))

percent_sinks_included = 0.01 # Sample 1% of sinks

# Mean/average one-way travel time in a hunting trip (via Levi et al.)
mean_hunt_travel_h = 3.62
# convert to sigma used to parameterize Rayleigh distribution
sigma = mean_hunt_travel_h/((math.pi/2)**0.5)

# Start and end points (for when the code fails and needs to restart at an intermediate level)
# Set to "False" to ignore this
# start_density = 105.253
start_density = 1660.512
#end_density = np.inf
end_density = False

#------------------------------------------------------------------------
#               Define classes
#------------------------------------------------------------------------

def samplePopDensity(srcfile, shpfile, targetfile, pbuff = 5):
    """
    Sample population density raster (with 5 pixel buffer) cost-distance based isochrones from a resistance surface.
    Inputs: srcfile - filepath (as string) to raster of population density
        shpfile - filepath (as string) to population density shapefile (as multipoints)
        targetfile - filepath (as string) with resolution/transformation equal to the resistance surface
        pbuff - buffer size in # of pixels, defaults 5 pixels (equal to ~5km with this script's srcfile resolution)
        
    Outputs a 5-column numpy array:
        0: Population density summed within buffer
        1: X-position of point according to the targetfile's resolution/transformation
        2: Y-position of above
        Col 3: Lon-coordinate (georeferenced)
        Col 4: Lat-coordinate (georeferenced)
    """
    src_ds=gdal.Open(srcfile) 
    gt=src_ds.GetGeoTransform()
    rb=src_ds.GetRasterBand(1)
    
    tar_ds=gdal.Open(targetfile) 
    tt=tar_ds.GetGeoTransform()
    
    ds=ogr.Open(shpfile)
    lyr=ds.GetLayer()
    
    # Iteration 1: ONLY VILLAGES, COMMENT OUT
        # Sabah settleme_1 ('Sabah'), Swk should include only 'Village'
    # lyr.SetAttributeFilter("settleme_1 IN ('Village', 'Sabah')")
    # Should create 1947 sources 
    
    # Iteration 2: Villages AND Other, no Longhouses
    lyr.SetAttributeFilter("settleme_1 != ('Long House')")
    # Should create 2704 sources
    
    # Initialize arrays for buffered pop density and coordinate location in XY
    buf_popdensities = np.empty([len(lyr),1])
    idx = 0
    
    buf_popcoordinates = np.empty([len(lyr), 2])
    
    # First column = x
    # Second column = y
    buf_popxy = np.empty([len(lyr), 2])
    
    
    for feat in lyr:
        geom = feat.GetGeometryRef()
        mx = geom.GetGeometryRef(0).GetX()
        my = geom.GetGeometryRef(0).GetY()  #coord in map units
        buf_popcoordinates[idx,0], buf_popcoordinates[idx,1] = mx, my
        
        #Convert from map to pixel coordinates.
        #Only works for geotransforms with no rotation.
        tx = int((mx - tt[0]) / tt[1]) #x pixel (target)
        ty = int((my - tt[3]) / tt[5]) #y pixel (target)
        px = int((mx - gt[0]) / gt[1]) #x pixel (sample)
        py = int((my - gt[3]) / gt[5]) #y pixel (sample)
        buf_popxy[idx,0], buf_popxy[idx,1] = tx, ty
        # Sample population density at each settlement point, using a 5 pixel (~5km) buffer
        bufferval = np.array(rb.ReadAsArray(px,py,1,1,buf_xsize=pbuff, buf_ysize=pbuff))
        bufferval[bufferval < 0] = np.nan
        buf_popdensities[idx] = np.nansum(bufferval)
        idx += 1
        print idx
    
    # Remove any 0 population density points
    # Also remove any values not within the target extent
    
    popdens = np.delete(buf_popdensities, np.where((buf_popdensities <= 0) | (buf_popxy < 0)), axis = 0)
    popcoords = np.delete(buf_popcoordinates, np.where((buf_popdensities <= 0) | (buf_popxy < 0)), axis = 0)
    popxy = np.delete(buf_popxy, np.where((buf_popdensities <= 0) | (buf_popxy < 0)), axis = 0)
    
    pop_xyz = np.column_stack((popdens, popxy[:,0], popxy[:,1], popcoords[:,0], popcoords[:,1]))
    return pop_xyz

def isochronesToSinks(coordinate_array, multiple, multiple_number=0, file_path=node_path):
    """
    Calculate cost-distance based isochrones from a resistance surface.
    Inputs: coordinate_array - array of georeferenced source coordinates
        multiple - Boolean indicating whether the source ID is a duplicate (useful when converting an entire population density grid rather than discrete points)
        multiple_number - how many files have already been created for this     population density?
        file_path - where to save output files that contain source-sink pairs
        
    No output, writes files in two new folders, NodesTSV and NodesTXT:
        Nodes .txt files: contains XY positions of each source/sink point to be used by GFlow (format: "Xposition Yposition")
        Nodes .tsv files: contains source-->sink identifiers (format: "sourceID \tab\ sinkID", where ID is the line of the .txt file that contains the source/sink point). Note: this file identifies source and sink locations based on their position in the .txt file: the first point listed in the .txt file is identified as point 1 in the .tsv file
    """
    global node_count
    global node_total
    # Creates output .tsv files into NodesTSV folder
    outpath_tsv=file_path+"/NodesTSV/"
    outpath_txt=file_path+"/NodesTXT/"
    source_id = float(np.unique(coordinate_array[:,0]))
    
    coordinates = []
    for c in coordinate_array:
        coordinates.append(float(c[3]))
        coordinates.append(float(c[4]))

    sources_XY = XYtoPixels(np.column_stack((coordinate_array[:,3], coordinate_array[:,4])))
    print(sources_XY)
    sys.exit
    # Starting node positions for the TSV node-to-node file
    source_position = 0
    sink_position = len(sources_XY)
    
    # Uses grass to create a cost surface for that source based on the resistance map input
    grass.run_command('r.cost', input='resist_input', output='temp_cost', start_coordinates=coordinates, flags='', overwrite=True, max_cost=(mean_hunt_travel_h*2))# , memory=50, null_cost=100)
    cost_surface=grass.read_command('r.out.xyz', input='temp_cost')
    cost_surface_array=(np.array([[float(j) for j in i.split('|')] for i in cost_surface.splitlines()]))
    if cost_surface_array.shape[0]<10:
        print("Empty cost array (<10 sink nodes), skipping")
        node_count+=1
        return
    cost_surface_array[:,2]=np.absolute(cost_surface_array[:,2])
    cost_surface_sort=cost_surface_array[cost_surface_array[:,2].argsort()]
    if multiple==False:
        filename="nodes_"+str(source_id)+"_0"
    elif multiple==True:
        filename="nodes_"+str(source_id)+"_"+str(multiple_number)
    print("Writing source node(s)")
    source_node_counter=1
    try:
        os.remove(outpath_txt+filename+'.txt')
        os.remove(outpath_txt+filename+'.tsv') 
        print("\n Removed exising files:", str(outpath_txt+filename))
    except Exception:
        pass
    with open((outpath_txt+filename+'.txt'), 'w') as f:
        for s in sources_XY:
            print str(source_node_counter)+"...",
            f.write("%d %d\n" %(s[1], s[0]))
            source_node_counter+=1
    print("Writing sink node(s)")
    for s in sources_XY:
        source_position += 1
        # if s != sources_XY[len(sources_XY)-1]:
        #    print(str(source_position)+"...")
        # else:
        #    print(str(source_position))
        rayleigh_values=np.random.rayleigh(scale=sigma, size=int(round(cost_surface_sort.shape[0]*percent_sinks_included)))
        index_list=[]
        for r in rayleigh_values:
            index=findNearest(cost_surface_sort[:,2], r)
            index_list.append(index)
            cost_surface_sort[index,2]=cost_surface_sort[index,2]*-1
        sink_coordinates=cost_surface_sort[index_list]
        sink_XY_o=XYtoPixels(sink_coordinates)
        # Remove any repeated sink locations:
        sink_XY = np.vstack({tuple(row) for row in sink_XY_o})
        # if np.array(sink_XY_o).shape != np.array(sink_XY).shape:
			# print('Original sink dim: %s x %s' %(np.array(sink_XY_o).shape[0],np.array(sink_XY_o).shape[1]))
			# print('New sink dim: %s x %s' %(sink_XY.shape[0], sink_XY.shape[1]))
        for sink in sink_XY:
            with open((outpath_txt+filename+'.txt'), 'a') as f:
                f.write("%d %d\n" %(sink[1], sink[0]))
                sink_position+=1
            with open((outpath_tsv+filename+'.tsv'), 'a') as f:
                f.write('%d\t%d\n' %(source_position, sink_position))
        node_count+=1

def findNearest(array, value):
    """
    Find the nearest value in the cost distance array to a randomly sampled Rayleigh value.
        Input: array - cost-distance array to sample from
               value - randomly sampled value to compare to
        Output: idx  - index of the cost distance value in 'array' closest to 'value'
    """
    idx = np.searchsorted(array, value, side='left')
    if idx > 0 and (idx==len(array) or math.fabs(value - array[idx-1]) < math.fabs(value-array[idx])):
        return idx-1
    else:
        return idx

### Helper function to create subfolders NodesTXT and NodesTSV
#     if they don't already exist
# If these folders do already exist, they can be deleted and repopulated
# If you restart the program, you will have to create new folders and merge
def ensureDir(file_path):
    exist=False
    try:
        os.makedirs(file_path+"/NodesTSV")
        os.makedirs(file_path+"/NodesTXT")
        exist=False
    except OSError as exception:
        print '......Node directories exist. /NodesTXT and /NodesTSV already present in %s' %(file_path)
        exist=True
        if exception.errno != errno.EEXIST:
            raise
    if exist==True and not overwrite_nodes:
        print "\n\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "! Files present in node folders - this can cause this process to fail in GFlow !"
        print "!            Should the files be deleted now?                                  !"
        print "!     (Yes=Y, No=N (appends to existing dirs), print affected path=P)          !"
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n\n"
        response=False

        ### FOR TESTING:
        var=raw_input("Is it safe to delete these files? (Y, N, or P, then Enter) --> ")
        while response==False:
            if var=='Y' or var=='y':
                print 'Deleting folder contents'
                shutil.rmtree(file_path+'/')
                os.makedirs(file_path+"/NodesTSV")
                os.makedirs(file_path+"/NodesTXT")
                response=True
            elif var=='N' or var=='n':
                return
            elif var=='P' or var=='p':
                print '/NodesTSV/ and /NodesTXT/ in %s' %(file_path)
                var=raw_input("Is it safe to delete these files? (Y or N, then Enter) --> ")
            else:
                var=raw_input("Please enter either Y, N, or P --> ")

def XYtoPixels(xy_array, convert=True):
# Change from source_input, originally, to resist_input, the new template file
    mapInfo=mapMatch('resist_input')
    ymax=mapInfo[0]
    xmin=mapInfo[3]
    pixelsize=mapInfo[4]
    pixel_array=[]
    if convert==True:
        for coordinate in xy_array:
            # Modify X, round up fit the grid (coordinate is on the centrepoint of the cell, such that cell 1,1 would be identified as 0.5, 0.5)
            x_grid=math.ceil((coordinate[0]-xmin)/pixelsize+0.0000001)
            # Modify Y
            y_grid=math.ceil((ymax-coordinate[1])/pixelsize+0.0000001)
            pixel_array.append([int(x_grid), int(y_grid)])
    else:
        for coordinate in xy_array:
            pixel_array.append([int(coordinate[0]), int(coordinate[1])])
    return pixel_array

def RasterConvert(raster=None, output_name=None):
    """
    Helper function to convert h/km rasters in decimal degrees into h/pixel
    Inputs: 
        raster - GRASS' name for the map to be converted from h/km into h/pixel
        output_name - GRASS name for the map to be created
    """
    if output_name==None or raster==None:
        print("Please provide input and output map names for raster conversion. Quitting, please try again.")
        return
    region_info=grass.read_command('g.region', rast=raster, flags='m')
    # Resolution in meters
    m_resolution=np.array([[str(j) for j in i.split('=')] for i in region_info.splitlines()])
    ns_res_m=float(m_resolution[6][1])
    ew_res_m=float(m_resolution[7][1])  
    if ns_res_m==ew_res_m:
        m_res=ns_res_m
        next
    elif round(ns_res_m,2)==round(ew_res_m, 2) or average_res:
#        print("Resolution (m) of the NS/EW directions within .01 decimal places. RasterConvert() will average them.")
        m_res=(ns_res_m+ew_res_m)/2
    else:
        print("\n\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!              Resolution (m) of the NS/EW directions NOT within 0.01 decimal places:        !\n!     NS resolution (m): %s   |   EW resolution (m): %s                   !" %(ns_res_m, ew_res_m))
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n")
        print("! Do you want RasterConvert() to average them, or do you want to quit? ([A]verage or [Q]uit) !\n\n")
        answer=False
        avg=raw_input('[A] or [Q], then Enter --> ')
        while answer==False:
            if avg=='A' or avg=='a':
                print("Averaging...")
                m_res=(ns_res_m+ew_res_m)/2
                answer=True
            elif avg=='Q' or avg=='q':
                sys.exit('Quitting.')
                answer=True
            else:
                avg=raw_input('Answer not recognized. Please type either [A] or [Q], and try again --> ')
    if m_res < 1:
      # Convert from decimal degrees into meters
      m_res2 = 111139*(m_res)
      print("Resolution in m is equal to %s, assuming this is incorrect and applying %s instead" %(m_res, m_res2))
      m_res = m_res2
    # Create the conversion value to translate hours/km as our pixel value to h/pixel
    # based on m_res being the meter-wise resolution of the map
    pixel_m_conversion=m_res/1000
    #return pixel_m_conversion, m_res
    exp=str(output_name)+' = '+raster +'*'+str(pixel_m_conversion)
    map_conversion=grass.run_command('r.mapcalc', expression=exp, overwrite=True)

### Function to return information about raster maps (xmin/max, ymin/max, resolution, etc.)
# If two maps are provided with the c=True flag, this function (c)ompares their values; otherwise this function takes one map and returns its statistics.
def mapMatch(map1, map2=None, c=False):
    stats_string_1=grass.read_command('r.info', map=map1, flags='g')
    stats_array_1=np.array(stats_string_1.split('\n'))
    if c==True:
        print("...Comparing maps")
        stats_string_2=grass.read_command('r.info', map=map2, flags='g')
        stats_array_2=np.array(stats_string_2.split('\n'))
        for i in range(0,8):
            if round(float(stats_array_1[i].split('=')[1]),2)!=round(float(stats_array_2[i].split('=')[1]),2):
                print("The two input maps do not match (within 2 decimal points). Mismatched parameters (from GRASS' r.info):")
                print("Map 1: %s | Map 2: %s" %(stats_array_1[i], stats_array_2[i]))
                answer=False
                while answer==False:
                    ignore=raw_input("Would you like to ignore these problems? (Y or N (quit), then Enter) -->")
                    if ignore=='Y' or ignore=='y':
                        answer=True
                        print("...ignoring. Region will be assigned according to resistance map parameters.")
                    elif ignore=='N' or ignore=='n':
                        answer=True
                        sys.exit("\nPlease address these discrepancies and try again.")
    # If the maps match, or if only one map is provided, return map statistics
        print("...maps match. Returning statistics for for map \'%s\'" %(map1))
    y_max=float(stats_array_1[0].split('=')[1])
    y_min=float(stats_array_1[1].split('=')[1])
    x_max=float(stats_array_1[2].split('=')[1])
    x_min=float(stats_array_1[3].split('=')[1])
    ns_res=float(stats_array_1[4].split('=')[1])
    ew_res=float(stats_array_1[5].split('=')[1])
    if ns_res!=ew_res:
        if round(ns_res,5)!=round(ew_res,5):
            ns_res=np.mean([ns_res,ew_res])
    return (y_max, y_min, x_max, x_min, ns_res)

    #------------------------------------------------------------------------
    #               Program start!
    #------------------------------------------------------------------------

def main():
    print(":::PROGRAM START:::")
    ensureDir(str(node_path))
    print("...Reading in resist/source file inputs...")
    print("...Setting GRASS region to resistance map projection...")
    grass.run_command('r.in.gdal', input=resist_infile, output='resist_input_o', flags='e', overwrite=True)
    
    ### REPLACE
    # grass.run_command('r.in.gdal', input=sources_infile, output='source_input', flags='oe', overwrite=True)
    
    # Convert from kph to hours per pixel
    RasterConvert('resist_input_o','resist_input')
    
    # Focus the region of the GRASS session and apply a mask
    
    grass.run_command('g.region', rast='resist_input')
    grass.run_command('r.mask', rast='resist_input', overwrite = True)
    
    shppath = os.path.join('/home/mairin/Documents/GradSchool/Research/CircuitTheory_Borneo/PRSB_Revision2/distanceToPopulation/')
    shpfile = os.path.join(shppath,'Revision2_PopulationSources','sab_swk_merged_gazetteervillages_sep9.shp')
    # For the sake of testing (Sept 9, 2019) - identify nodes that are not Longhouses
    #   THIS CORRESPONDS TO A MODIFICATION ON LINES 110ish
    srcpath = os.path.join('/home/mairin/Documents/GradSchool/Research/CircuitTheory_Borneo/PRSB_Revision2/Revised_CTMapCreation/SourceSinks')
    srcfile = os.path.join(srcpath,'asuds00ag.tif')
    
    sources = samplePopDensity(srcfile, shpfile, targetfile = str(resist_infile))
    
    global node_total
    node_total = len(sources)
    global node_count
    node_count = 1
    o_sources = None
    
    # Unique densities:
    u_densities = np.unique(sources[:,0])
    
    try:
        len(u_densities)
    except(TypeError):
        u_densities = [u_densities]
    
    print("\n" + str(len(u_densities)) + " unique population densities")
    u_count_start = 0
    u_count_current = 0
    
    if start_density:
        if end_density:
            include = np.where((u_densities >= start_density) & (u_densities <= end_density))
        else:
            include = np.where((u_densities >= start_density))
        # print(include[0])
        # print(u_count_start)
        
        # u_count_start, u_count_current = int(include[0][0]), int(include[0])
        u_count_start = len(np.where((sources[:,0] < start_density))[0])
        u_densities = u_densities[include]
        # print(np.where(sources[:,0] not in u_densities))
        # return
        u_count_current = u_count_start
        # print(u_count_start)
        
    # For each unique population density value, create a set of coordinates
    #   for all points that have that density and make files

    for d in u_densities:
        coordinate_set = sources[sources[:,0]==d]
#        print(coordinate_set)
#        print(len(coordinate_set))
        if len(coordinate_set) < 50:
            u_count_start=u_count_current+1
            u_count_current = u_count_start + len(coordinate_set)
            print("......Converting nodes of density %s to isochrones [nodes %s-%s of %s]" %(round(d,3), u_count_start, u_count_current, node_total))
#           print(coordinate_set)
            isochronesToSinks(coordinate_set, multiple=False,file_path=str(node_path))
            # Otherwise, if > 50, create multiple files with suffixes for ensuring that .tsv and .txt files are paired
        else:
                total_length = len(coordinate_set)
                current_length = 0
                m = 0
                while current_length < total_length:
                    if current_length + 50 < total_length:
                        u_count_current = u_count_start + 50
                        print("......Converting nodes of density %s to isochrones [nodes %s-%s of %s]" %(round(d,3), u_count_start, u_count_current, node_total))
                        coordinate_subset = coordinate_set[current_length:current_length+50]
                        isochronesToSinks(coordinate_subset, multiple=True, multiple_number=m,file_path=str(node_path))
                        current_length += 50
                        m += 1
                        u_count_start = u_count_current + 1
                    elif current_length + 50 >= total_length:
                        u_count_current=u_count_start+(total_length-current_length)
                        print("......Converting nodes of density %s to isochrones [nodes %s-%s of %s]" %(round(d,3), u_count_start, u_count_current, node_total))
                        coordinate_subset=coordinate_set[current_length:total_length]
                        isochronesToSinks(coordinate_subset, multiple=True, multiple_number=m,file_path=str(node_path))
                        current_length = total_length
                        u_count_start = u_count_current+1
    print(":::PROGRAM END, RANDOM SINK POINTS GENERATED:::")

if __name__ == "__main__":
    sys.exit(main())
