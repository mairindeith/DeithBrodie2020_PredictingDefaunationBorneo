#!/usr/bin/env python

# Python script to create isochrones and sample sink points from within these isochrones
# This script should be run from within a GRASS session
# e.g. After opening a GRASS session in Terminal:
#    $ python grass_isochrone_sampler.py

# Import pertinent grass libraries (raster and general)
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import raster as r
import grass.script as grass
import math
import sys, os
import errno
import numpy as np
import shutil

### Program-wide parameters
# Resistance file should be resistance (measured in hours/pixel)
# If needed, the function "RasterConvert" can convert a h/km map into an h/pixel map
resist_infile=os.path.join("PATH_TO_INPUT_FILES","RESISTANCE_MAP_HOURS_PER_PIXEL.asc")
# Sources file should be a raster of sources, where cell value=population density in that pixel
sources_infile=os.path.join("PATH_TO_INPUT_FILES","POPULATION_DENSITY_N_PER_PIXEL.asc")
# Path to a folder where you want your NodesTSV and NodesTXT folders to be created
# If these folders already exist, this code DELETES THESE and recreates them
#    Change the "SUBFOLDER" to a different value in each instance if you want to
#    run this code in multiple instances
node_path=os.path.join("PATH_TO_SOURCE_NODES","SUBFOLDER")

#Range of population densities to include in the node list
begin_point=[1] # lowest population density to consider a discrete point
end_point=[2000] # maximum value if applicable, otherwise any number larger
                 # than maximum population density
# Depreciated - this is only appropriate when running in parallel, keep =[1]
parallel_list=[1]

percent_sources_included=1 # what % of source points should be included?
percent_sinks_included=0.01 # what % percent of pixels within the travel time
                            # window should be included?

# Mean/average one-way travel time in a hunting trip
mean_hunt_travel_h=3.62
# convert to sigma for Rayleigh distribution
sigma=mean_hunt_travel_h/((math.pi/2)**0.5)

global node_count
global node_total

#------------------------------------------------------------------------
#               Define classes
#------------------------------------------------------------------------

### Function to calculate isochrones
# Inputs: coordinates (georeferenced) of sources, costDist in hours
# No output, writes files in two new folders, NodesTSV and NodesTXT:
#       Nodes .txt files: contains XY positions of each source/sink point to be used by GFlow (format: "Xposition Yposition")
#       Nodes .tsv files: contains source-->sink identifiers (format: "sourceID \tab\ sinkID", where ID is the line of the .txt file that contains the source/sink point)
#           Note: this file identifies source and sink locations based on their position in the .txt file: the first point listed in the .txt file is identified as point 1 in the .tsv file

def isochronesToSinks(coordinate_array, multiple, multiple_number=0, file_path=node_path):
    global node_count
    global node_total
    source_id=coordinate_array[0][2]
    sources_XY=XYtoPixels(coordinate_array)
    # Starting node positions for the TSV node-to-node file
    source_position=0
    sink_position=len(sources_XY)
    coordinates=[]
    for c in coordinate_array:
        coordinates.append(float(c[0]))
        coordinates.append(float(c[1]))
# Creates output .tsv files into NodesTSV folder
    outpath_tsv=file_path+"/NodesTSV/"
    outpath_txt=file_path+"/NodesTXT/"
    # Uses grass to create a cost surface for that source based on the resistance map input
    grass.run_command('r.cost', input='resist_input', output='temp_cost', start_coordinates=coordinates, flags='', overwrite=True, max_cost=(mean_hunt_travel_h*2))# , memory=50, null_cost=100)
    cost_surface=grass.read_command('r.out.xyz', input='temp_cost')
    cost_surface_array=(np.array([[float(j) for j in i.split('|')] for i in cost_surface.splitlines()]))
    cost_surface_array[:,2]=np.absolute(cost_surface_array[:,2])
    if cost_surface_array.shape[0]<10:
        print "Empty cost array (<10 sink nodes), skipping"
        node_count+=1
        return
    cost_surface_sort=cost_surface_array[cost_surface_array[:,2].argsort()]
#    print cost_surface_sort[0:5,:]
#    print "..."
#    print np.mean(cost_surface_sort[:,2])
#    print cost_surface_sort[cost_surface_sort.shape[0]-5:cost_surface_sort.shape[0],:]
    if multiple==False:
        filename="nodes_"+str(source_id)+"_0"
    elif multiple==True:
        filename="nodes_"+str(source_id)+"_"+str(multiple_number)
    #    sys.exit("Yep, it's true.")
    print "Writing source node",
    source_node_counter=1
    with open((outpath_txt+filename+'.txt'), 'w') as f:
        for s in sources_XY:
            print str(source_node_counter)+"...",
            f.write("%d %d\n" %(s[1], s[0]))
            source_node_counter+=1
    print 'done.'
    print "Writing sink node",
    for s in sources_XY:
        source_position+=1
        if s != sources_XY[len(sources_XY)-1]:
            print str(source_position)+"...",
        else:
            print str(source_position)
        rayleigh_values=np.random.rayleigh(scale=sigma, size=int(round(cost_surface_sort.shape[0]*percent_sinks_included)))
        index_list=[]
#        print rayleigh_values[rayleigh_values.argsort()]
        for r in rayleigh_values:
            index=findNearest(cost_surface_sort[:,2], r)
            index_list.append(index)
            cost_surface_sort[index,2]=cost_surface_sort[index,2]*-1
        sink_coordinates=cost_surface_sort[index_list]
        # np.savetxt('sink_coordinates_%s.txt' %(source_id), sink_coordinates)
        # np.savetxt('source_coordinates_%s.txt' %(source_id), coordinate_array)
        # grass.run_command('v.in.ascii', input='sink_coordinates_%s.txt' %(source_id),
        # separator=' ', overwrite=True, output='sinks_rayleigh', quiet=True)
        # grass.run_command('v.in.ascii', input='source_coordinates_%s.txt' %(source_id),
        # separator=' ', overwrite=True, output='sources', quiet=True)
        sink_XY_o=XYtoPixels(sink_coordinates)
        # Line to remove any repeated sink locations:
        sink_XY = np.vstack({tuple(row) for row in sink_XY_o})
        if np.array(sink_XY_o).shape != np.array(sink_XY).shape:
			print 'Original sink dim: %s x %s' %(np.array(sink_XY_o).shape[0],np.array(sink_XY_o).shape[1])
			print 'New sink dim: %s x %s' %(sink_XY.shape[0], sink_XY.shape[1])
#        sys.exit("Check unique")
        for sink in sink_XY:
            with open((outpath_txt+filename+'.txt'), 'a') as f:
                f.write("%d %d\n" %(sink[1], sink[0]))
#                if sink[1]==201 & sink[0]==512:
#					print "\n\n\n\n\n 201-512 error \n\n\n\n"
                sink_position+=1
            with open((outpath_tsv+filename+'.tsv'), 'a') as f:
                f.write('%d\t%d\n' %(source_position, sink_position))
        node_count+=1
    print "...done."


### Function to find the nearest value in the cost distance array to the
#   randomly sampled Rayleigh value
# Called by isoChronesToSinks
#
# This function is used in sampling points - the array must be 1D and sorted by the relevant value
# Input: array - cost-distance array to sample from
#        value - randomly sampled value to compare to
# Output: idx - index of the cost distance value in 'array' closest to 'value'

def findNearest(array, value):
    idx = np.searchsorted(array, value, side='left')
    if idx > 0 and (idx==len(array) or math.fabs(value - array[idx-1]) < math.fabs(value-array[idx])):
        return idx-1
    else:
        return idx

### Function to convert georeferenced coordinates to X-Y sink points for GFlow
# Called by the isoChronesToSinks function
# Input: xy_array - Georeferenced coordinate array from GRASS's r.out.xyz function
# Output: pixel_array - the same array with non-georeferenced x-y positions for GFlow input

def XYtoPixels(xy_array):
    mapInfo=mapMatch('source_input')
    ymax=mapInfo[0]
    xmin=mapInfo[3]
    pixelsize=mapInfo[4]
    pixel_array=[]
    for coordinate in xy_array:
        # Modify X, round up fit the grid (coordinate is on the centrepoint of the cell, such that cell 1,1 would be identified as 0.5, 0.5)
        x_grid=math.ceil((coordinate[0]-xmin)/pixelsize+0.0000001)
        # Modify Y
        y_grid=math.ceil((ymax-coordinate[1])/pixelsize+0.0000001)
        pixel_array.append([int(x_grid), int(y_grid)])
    return pixel_array

### Helper function to create subfolders NodesTXT and NodesTSV
#     if they don't already exist
# If these folders do already exist, they can be deleted and repopulated
# If you restart the program, you will have to create new folders and merge
def ensureDir(file_path):
    exist=False
    try:
        os.makedirs(file_path+"/NodesTSV")
        os.makedirs(file_path+"/NodesTXT")
        print '......Creating node directories in %s.' %(file_path)
        exist=False
    except OSError as exception:
        print '......Node directories exist. /NodesTXT and /NodesTSV already present in %s' %(file_path)
        exist=True
        if exception.errno != errno.EEXIST:
            raise
    if exist==True:
        print "\n\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "! Files present in node folders - this can cause this process to fail in GFlow !"
        print "!            Can the files be deleted now? Please check first!                 !"
        print "!             (Yes=Y, No=N, print affected path=P)                             !"
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
                sys.exit("Please move these folders, or change the node_path variable in the Python script.")
            elif var=='P' or var=='p':
                print '/NodesTSV/ and /NodesTXT/ in %s' %(file_path)
                var=raw_input("Is it safe to delete these files? (Y or N, then Enter) --> ")
            else:
                var=raw_input("Please enter either Y, N, or P --> ")

### Helper function to convert h/km rasters in decimal degrees into h/pixel
# Inputs: raster - GRASS' name for the map to be converted from h/km into h/pixel
#         output_name - GRASS name for the map to tbe created
def RasterConvert(raster=None, output_name=None):
    if output_name==None or raster==None:
        print "Please provide input and output map names for raster conversion. Quitting, please try again."
        return
    region_info=grass.read_command('g.region', rast=raster, flags='m')
    # Resolution in meters
    m_resolution=np.array([[str(j) for j in i.split('=')] for i in region_info.splitlines()])
    ns_res_m=float(m_resolution[6][1])
    ew_res_m=float(m_resolution[7][1])
    if ns_res_m==ew_res_m:
        m_res=ns_res_m
        next
    elif round(ns_res_m,2)==round(ew_res_m, 2):
        print "Resolution (m) of the NS/EW directions within .01 decimal places. RasterConvert() will average them."
        m_res=(ns_res_m+ew_res_m)/2
    else:
        print "\n\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!              Resolution (m) of the NS/EW directions NOT within 0.01 decimal places:        !\n!     NS resolution (m): %s   |   EW resolution (m): %s                   !" %(ns_res_m, ew_res_m)
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
        print "! Do you want RasterConvert() to average them, or do you want to quit? ([A]verage or [Q]uit) !\n\n"
        answer=False
        avg=raw_input('[A] or [Q], then Enter --> ')
        while answer==False:
            if avg=='A' or avg=='a':
                print "Averaging..."
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
      print "Resolution in m is equal to %s, assuming this is incorrect and applying %s instead" %(m_res, m_res2)
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
        print "...Comparing maps",
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
                        print "...ignoring. Region will be assigned according to resistance map parameters."
                    elif ignore=='N' or ignore=='n':
                        answer=True
                        sys.exit("\nPlease address these discrepancies and try again.")
    # If the maps match, or if only one map is provided, return map statistics
        print "...maps match. Returning statistics for for map \'%s\'" %(map1)
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
    print ":::PROGRAM START:::"
    print "...Reading in resist/source file inputs..."
    print "\n\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "!   This overwrites maps resist_input & source_input in GRASS session  !"
    print "!     (local files outside the session will not be affected)           !"
    print "!        Are you sure you want to continue? (Yes=Y, No=N)              !"
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n\n"
    var=raw_input("Do you want to continue? (Y or N, then Enter) --> ")
    valid_response=False
    while valid_response==False:
        if var=='N' or var=='n':
            sys.exit("...exiting.")
        elif var=='Y' or var=='y':
            print ("...proceeding...")
            valid_response=True
        else:
            var=raw_input("Please enter either Y or N --> ")
            print "...Confirming/creating source/sink node folder structure..."

    ensureDir(str(node_path))
# Uncomment this to see a list of all rasters loaded into the GRASS session
#     g.list(type='raster')
    print "...Setting GRASS computational region to resistance map projection... "
    grass.run_command('r.in.gdal', input=resist_infile, output='resist_input_o', flags='e', overwrite=True)
    grass.run_command('r.in.gdal', input=sources_infile, output='source_input', flags='oe', overwrite=True)

    RasterConvert('resist_input_o','resist_input')
    # Compare the extent and resolution of source_input and resist_input
    maxY, minY, maxX, minX, res = mapMatch('source_input', 'resist_input', c=True)
    # Extending the region of the GRASS session based on the now checked resist_input
    grass.run_command('g.region', rast='resist_input')
    # Create an X-Y-Z file of the X and Y coordinates with population size Z in the input file
    print "...Converting sources to XYZ georeferenced points:"
    source_xyz=grass.read_command('r.out.xyz', input='source_input')
    # Remove sources less than the cutoff
    o_sources=np.array([[float(j) for j in i.split('|')] for i in source_xyz.splitlines()])

    for parallel in parallel_list:
        print "Generating node files for population densities between lower_cutoff and upper_cutoff"
        print "...Lower cutoff: %s\n...Upper cutoff: %s" %(begin_point, end_point)
        delete=[]
        upper_delete_count=0
        lower_delete_count=0

        lower_cutoff=begin_point
        upper_cutoff=end_point
        # Count number of nodes to be deleted/included
        for i in range(0, len(o_sources)):
            if(o_sources[i,2]<lower_cutoff):
                delete.append(i)
                lower_delete_count+=1
            if(o_sources[i,2]>upper_cutoff and upper_cutoff!=False):
                delete.append(i)
                upper_delete_count+=1

        sources=np.delete(o_sources,delete,0)

        print "\n...%s nodes below density cutoff of %s\n...%s nodes above density cutoff of %s\n...(out of %s originally); removed for a total of %s source nodes." %(lower_delete_count, lower_cutoff, upper_delete_count, upper_cutoff, len(o_sources), len(sources))

        global node_total
        node_total=len(sources)
        global node_count
        node_count=1
        o_sources=None
        # Unique population density values:
        u_densities=np.unique(sources[:,2])
        try:
            len(u_densities)
        except(TypeError):
            u_densities=[u_densities]
        # For each unique population density value, create a set of coordinates
        #   for all points that have that density and make files
        u_count_start=1
        u_count_current=0
        for d in u_densities:
            coordinate_set_o=[]
            for s in sources:
                if s[2]==d:
                    # Make a list, coordinate_set_o, populated with the
                    # coordinates of points where density=the current population density
                    u_count_current+=1
                    coordinate_set_o.append(s)
            coordinate_set_np=np.array(coordinate_set_o)
            # print coordinate_set_np.shape[0]
            coordinate_set=coordinate_set_np[np.random.choice(coordinate_set_np.shape[0], size=int(math.ceil(percent_sources_included*coordinate_set_np.shape[0])), replace=False),:]
            # print coordinate_set.shape[0]
            # print "\n"+str(coordinate_set)

            # If there are less than 50 coordinates for a given population density,
            #   make a single .txt and a single .tsv file for that population density
            if len(coordinate_set)<50:
                print "......Converting nodes of density %s to isochrones [nodes %s-%s of %s]" %(round(d,3), u_count_start, u_count_current, node_total)
                # Call isoChronesToSinks with the set of coordinates for the current set of sources
                isochronesToSinks(coordinate_set, multiple=False,file_path=str(node_path))
                u_count_start=u_count_current+1
            # Otherwise, if > 50, create multiple files with suffixes for ensuring that .tsv and .txt files are paired
            else:
                total_length=len(coordinate_set)
                current_length=0
                m=0
                while current_length<total_length:
                    if current_length+50<total_length:
                        u_count_current=u_count_start+50
                        print "......Converting nodes of density %s to isochrones [nodes %s-%s of %s]" %(round(d,3), u_count_start, u_count_current, node_total)
                        coordinate_subset=coordinate_set[current_length:current_length+50]
                        isochronesToSinks(coordinate_subset, multiple=True, multiple_number=m,file_path=str(node_path))
                        current_length+=50
                        m+=1
                        u_count_start=u_count_current+1
                    elif current_length+50>=total_length:
                        u_count_current=u_count_start+(total_length-current_length)
                        print "......Converting nodes of density %s to isochrones [nodes %s-%s of %s]" %(round(d,3), u_count_start, u_count_current, node_total)
                        coordinate_subset=coordinate_set[current_length:total_length]
                        isochronesToSinks(coordinate_subset, multiple=True, multiple_number=m,file_path=str(node_path))
                        current_length=total_length
                        u_count_start=u_count_current+1
    print(":::PROGRAM END, RANDOM SINK POINTS GENERATED:::")

if __name__ == "__main__":
    sys.exit(main())
