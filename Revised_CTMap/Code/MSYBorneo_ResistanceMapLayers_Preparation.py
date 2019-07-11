#!/usr/bin/env python

#------------------------
# Code to standardize input maps used to generate
#    a human-based resistance map for the Malaysian states of Borneo
# This code prepares and saves standardized map files for further
#   processing
#
# First created by Mairin Deith on Nov12 2017
# Last edited for revised mapping effort on Jul4 2019
#------------------------

# Import libraries
import numpy as np
from osgeo import gdal, gdal_array, osr, ogr, gdalconst
import geopandas as gpd # for shapefile clipping
import richdem as rd
import elevation
import gc
from PIL import Image
import os, sys
import math
from pathlib import Path
import codecs
from datetime import datetime

startTime = datetime.now()

# Global parameters

output_path = Path('/home/mairin/Documents/GradSchool/Research/CircuitTheory_Borneo/PRSB_Revision2/Revised_CTMap/ResistanceMap')
input_path = output_path / 'InputSpatialData'
modified_path = output_path / 'ModifiedSpatialData'
global_nodata = -9999

# Helper functions
def rasterize(vector_data_path, outfile_path, template_path=None, cols=None, rows=None, geo_transform=None, projection=None, target_value=1):
	"""Rasterize the given vector (wrapper for gdal.RasterizeLayer)."""
	if(template_path != None):
		template = gdal.Open(template_path)
		cols = template.RasterXSize
		rows = template.RasterYSize
		geo_transform = template.GetGeoTransform()
		projection = template.GetProjection()
		template.GetRasterBand(1).SetNoDataValue(global_nodata)
	else:
		if(projection == None or geo_transform == None or cols == None or rows == None):
			print("ERROR: 'rasterize' needs either a template raster, or rows, columns, projection, and geotransformation must be provided.")
			return()
	data_source = gdal.OpenEx(vector_data_path, gdal.OF_VECTOR)
	layer = data_source.GetLayer(0)
	driver = gdal.GetDriverByName('GTiff')  # In memory dataset
	target_ds = driver.Create(outfile_path, cols, rows, 1, gdal.GDT_UInt16)
	target_ds.SetGeoTransform(geo_transform)
	target_ds.SetProjection(projection)
	gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[target_value])
	target_ds.GetRasterBand(1).SetNoDataValue(0)
	target_ds.GetRasterBand(1).SetNoDataValue(global_nodata) 
	target_ds.FlushCache()
	target_ds = None

def clipraster(source_raster, output_file, polygon):
	"""Clips a raster given an input polygon using a call to the operating system's command line.
	Arguments: source_raster, the raster that is to be clipped;
	output_file, the path + name of the clipped raster to be saved;
	polygon, a path to the polygon used in sclipping.
	Warning: this will overwrite any clipped rasters given the same name/path.
	"""
	os.system('gdalwarp -overwrite -cutline ' + str(polygon) + ' -crop_to_cutline ' + str(source_raster) + ' ' + str(output_file))

def resample(input_file, template_file, output_file, method = None):
	"""Resamples the given raster based on template file using method.
	If not method provided, uses nearest neighbour.
	Other methods: 'mode','mean','bilinear','max', and 'min'"""
	input = gdal.Open(input_file, gdalconst.GA_ReadOnly)
	inputProj = input.GetProjection()
	inputTrans = input.GetGeoTransform()
	reference = gdal.Open(template_file, gdalconst.GA_ReadOnly)
	referenceProj = reference.GetProjection()
	referenceTrans = reference.GetGeoTransform()
	ref_na = reference.GetRasterBand(1).GetNoDataValue()
	x = reference.RasterXSize
	y = reference.RasterYSize
	driver= gdal.GetDriverByName('GTiff')
	output = driver.Create(output_file,x,y,1,reference.GetRasterBand(1).DataType)
	output.SetGeoTransform(referenceTrans)
	output.SetProjection(referenceProj)
	output.GetRasterBand(1).SetNoDataValue(global_nodata)
	# Note: reproject image maps NO_DATA to 0
	if method == None or method == 'neighbour':
		gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_NearestNeighbour)
	elif method == 'mode':
		gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_Mode)
	elif method == 'mean':
		gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_Average)
	elif method == 'max':
		gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_Max)
	elif method == 'min':
		gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_Min)
	elif method == 'bilinear':
		gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_Bilinear)
	else:
		print("Oops! Method must be one of 'neighbour' (default), 'mode', 'max','min', 'mean', or 'bilinear'.")
	del output

# Main code for processing
def main():
    print("Begin")

    ### ADMINISTRATIVE BOUNDARY
    boundaryShape_loc = input_path / 'AdminBoundaries' / 'MSY_adm' / 'MSY_Borneo_Merged.shp'

	### ROADS
	# 	Already clipped and rasterized to target resolution.
	#	1 = highway, 2 = logging road.
    roadraster_loc = input_path / 'BorneoRoadDensity' / 'rasterizedRoads.tif'

	### LAND COVER
	#	Forest cover by Gaveau et al. 2014
    forestcover_loc_o = input_path / 'BorneoForest_Gaveau' / 'Gaveau_ForestCover' / 'REGIONBorneo_FCDefDeg_1973to2010_CIFOR.tif'
    forestcover_loc_clip = modified_path / 'REGIONBorneo_FCDefDeg_1973to2010_CIFOR_clipped.tif'
	# Creates a new file after clipping
    print("...Clipping and resampling forest cover raster...")
    clipraster(source_raster = str(forestcover_loc_o),
	    output_file = str(forestcover_loc_clip),
	    polygon = str(boundaryShape_loc))
    fc_resample_loc = modified_path / "REGIONBorneo_FCDefDeg_1973to2010_CIFOR_clipped_rs.tif"
    resample(input_file = str(forestcover_loc_clip),
		template_file = str(roadraster_loc),
		output_file = str(fc_resample_loc))
    # 	Plantation locations by Gaveau et al. 2014
    print("...Clipping and rasterizing plantations...")
    plantation_loc_o = input_path / 'BorneoForest_Gaveau'/'Gaveau_plantations' / 'REGIONBorneo_IndustrialPlantation_2010_CIFOR.shp'
    # First clip to boundary
    plantation_o = gpd.read_file(str(plantation_loc_o))
    boundary = gpd.read_file(str(boundaryShape_loc))
    if (plantation_o.crs == boundary.crs):
	    print("Both layers are in the same crs!", plantation_o.crs, boundary.crs)
    plantation_clip = gpd.overlay(plantation_o, boundary, how = 'intersection')
    plantation_clip_loc = modified_path / 'GaveauPlantationClip.shp'
    plantation_clip.to_file(str(plantation_clip_loc))
    # Rasterize
    plantation_raster_loc = modified_path / 'GaveauPlantationClip_raster.tif'
    rasterize(vector_data_path = str(plantation_clip_loc),
    outfile_path = str(plantation_raster_loc),
    template_path = str(roadraster_loc))

    ###	DEM elevation file
    print("...Clipping and resampling DEM...")
    # Use new elevation data
    dem_loc = input_path/'DEM_STRM'/'MergedSTRMScenes.tif'
    dem_clip_loc = modified_path / "DEM_clip.tif"
    dem_resample_loc = modified_path / "DEM_rs.tif"

    clipraster(source_raster = str(dem_loc),
		output_file = str(dem_clip_loc),
		polygon = str(boundaryShape_loc))

    resample(input_file = str(dem_clip_loc),
    	template_file = str(roadraster_loc),
    	output_file = str(dem_resample_loc),
    	method = "bilinear")
    
    slope_loc = modified_path / "RDSlope_pc.tif"
    
    # Calculate slope using 111120 as the scale factor
    os.system('gdaldem slope -compute_edges -p -s 111120 ' + str(dem_resample_loc) + ' ' + str(slope_loc))
	
    ### RIVERS
    print("...Clipping and resampling Strahler stream order raster...")
    strahler_loc = input_path/'DEM_StreamData'/'DEM_strahlerorder.tif'
    strahler = gdal.Open(str(strahler_loc), gdal.GA_Update)
    strahler.GetRasterBand(1).SetNoDataValue(255)
    strahler.FlushCache()
    strahler = None
    strahler_resample_loc = modified_path/'StrahlerOrder_rs.tif'
    resample(input_file = str(strahler_loc),
	    template_file = str(roadraster_loc),
	    output_file = str(strahler_resample_loc),
	    method = 'max')
    print("Done!")
    print('Total running time: ' + str(datetime.now() - startTime))


if __name__ == "__main__":
	main()
