#------------------------
# Code to create a resistance surface for the Indonesian island of Buton
#
# Created by MD on Nov12 2017
# Last edited on Feb17 2018
#------------------------

#!/usr/bin/env python

# Import libraries
import numpy as np
from osgeo import gdal, gdal_array, osr, ogr
from PIL import Image
import os, sys
import math
import codecs

# Designate the out_path (path where files will be output) and input_path (the file path where geospatial data can be found)
out_path = '/home/mairin/Documents/UBC/ButonResistanceMapping'
input_path = out_path+'/InputSpatialData'

def rasterize(shapefile_infile, targetName, p_size):
	# This function identifies which ID field is the one being used
	# to generate values for the eventual raster file created from the
	# .shp file

	# Get short names for the infile
	(in_path, in_name) = os.path.split(shapefile_infile)
	short_in_name = in_name[0:len(in_name)-4]
	out_name = str(short_in_name + ".tif")

	print "Rasterizing " + in_name + " to create " + out_name

	driver = ogr.GetDriverByName('ESRI Shapefile')
	dataSource = driver.Open(shapefile_infile, 0)

	layer = dataSource.GetLayer()
	layer_defn = layer.GetLayerDefn()
	id_index=0
	found=False
	while found==False:
		try:
			id_field = layer_defn.GetFieldDefn(id_index).GetName()
			if id_field==targetName:
				found=True
			else:
				id_index+=1
		except AttributeError as e:
			print(e)
			found=True
			sys.exit("Error rasterizing %s, target ID field %s not found in shapefile." %(short_in_name, targetName))

	# id_field = layer_defn.GetFieldDefn(id_index).GetName()
	print "ID field is: " + id_field
	feature = layer.GetNextFeature()

	unique = []

	while feature:
		feature_id = feature.GetField(id_index)
		if feature_id in unique:
			continue
		else:
			unique.append(feature_id)
			#print "Feature ID is: " + str(feature_id)
		feature = layer.GetNextFeature()

	print "Number of unique values: " + str(len(unique))
	for q in range(0,len(unique)):
		print "...Unique value " + str(i) + ": " + str(unique[i])

	# Create a new raster based on the id_field selected and the pixel size given in main:

	print "Generating new raster using gdal.rasterize:"
	os.system('gdal_rasterize -a_nodata -999 -a %s -te 122.5578708 -5.7024932 123.2101382 -4.4178918 -tr %s %s -l %s %s %s/%s' %(id_field, p_size, p_size, short_in_name, shapefile_infile, out_path, out_name))

	# Return the raster path that was created
	return str(out_path + "/" +  out_name)

	layer = None
	feature = None
	unique = None

# Inputs:
	### ROADS (need rasterization)
		# File: OSM_ButonClip.shp, converted into raster format in OSM_ButonClip.tif
		# Downloaded from OSMaps, field "highway" identifies road type
# Speed; classifications (fields - rdsV_ID: 'highway')
	### x 1: 'footway' - For designated footpaths; i.e., mainly/exclusively for pedestrians. This includes walking tracks and gravel paths.
	### x 2: 'path' - A non-specific path. Use highway=footway for paths mainly for walkers
	### x 3: 'pedestrian' - For roads used mainly/exclusively for pedestrians in shopping and some residential areas which may allow access by motorised vehicles only for very limited periods of the day
	### 60; 4: 'primary'  - The next most important roads in a country's system after motorways and highways, often link larger towns
	### 30; 5: 'primary_link' - The link roads (sliproads/ramps) leading to/from a primary road from/to a primary road or lower class highway.
	### 15; 6: 'residential' - Roads which serve as an access to housing, without function of connecting settlements. Often lined with housing.
	### 15; 7: 'road' - A road/way/street/motorway/etc. of unknown type
	### 30; 8: 'secondary' - The next most important roads in a country's system after primary (often link towns)
	### 15; 9: 'service' - alley / driveway / parking_aisle etc.
	### 10: 'steps' - For flights of steps (stairs) on footways
	### 15; 11: 'tertiary' - The next most important roads in a country's system after secondary (Often link smaller towns and villages)
	### 15; 12: 'track' - Roads for mostly agricultural or forestry uses.
	### 15; 13: 'unclassified' - The least most important through roads in a country's system, often link hamlets

# In calculation file: "WalkingSpeed_RoadOverlay":
#	4 = 60kph
#	5,7,8 = 40 kpg
#   6, 9, 11, 12, 13

	### ELEVATION (SRTM downloaded from ??? - received from Adrienne)
		# Elevation is in meters, ranges from -8 to 1050m above sea level
		# Apply Tobler's walking function according to land cover type after converting into slope
	### LAND COVER - from GLOBCover:
		# File: GLOBCOVER_ButonClip.tif
		# Downloaded from ESA's GlobCover database (http://due.esrin.esa.int/page_globcover.php)
		# Classifications:
			### 11: Post-flooding or irrigated croplands (or aquatic)
			### 14: Rainfed croplands
			### 20: Mosaic croplands (50-70% crops, 20-50% grassland/shrubland/forest)
			### 30: Mosaic vegetation (grassland/shrubland/forest, 50-70%)/cropland(20-50%)
			### 40: Closed to open broadleaf evergreen/semi-deciduous forest (>5m)
			### 110: Mosaic forest or shrubland (50-70%, grassland 20-50%)
			### 130: Closed to open (>15%) broadleaved forest, shrubland (<5m)
			### 160: Closed to open (>15%) broadleaved forest, regularly flooded temporarily
			### 170: Closed (>40%) broadleaved forest or shrubland permanently flooded, saline/brackish water
			### 210: Water bodies

def main():
   print('BEGIN')
   # Note: the pixel size is set to one in this case, but will depend on your map resolution
   pixel_size = 1 # Standardize pixel size for all rasters
   			      #      modified to equal 1 for the sake of GFlow (smaller pixel sizes caused GFlow to fail)
   # Locations for the various spatial datasets used to create the resistance surface
	adminShape_loc = input_path+'/boundaries_shapefile.shp'
	dem_loc = input_path+'/srtm_clip.tif'
	slope_loc_o = input_path+'/srtm_slope_clip.tif'
	# print("Calculating slope from %s" % dem_loc)
	# os.system('gdaldem slope -p -s 111120 %s %s' %(dem_loc, slope_loc_o))


	osmRoad_loc = input_path+'/OSM_ButonClip_RoadID.shp'
	osmRoad_rasterLoc = rasterize(osmRoad_loc, 'rdsV_ID', pixel_size)
	sys.exit("rasterize?")
	strahler_loc_o = '/home/mairin/deithvader/Documents/UBC/MScThesis/ASTERData/DEM_merged_attempt.tif'

######## MODIFIED FOR TESTING
#	for q in range(1,5):
	for q in range(1,2):
		print('................SECTION %s OF 4................') %(str(q))
#		landuse_quarters = quarters(landuse_loc_o, q, pixel_size, "raster")
#		landuse_q=landuse_quarters[0]
#		landuse_q_dim=landuse_quarters[1]
#		print(str(landuse_q))
#		print(str(landuse_q_dim))

		print('...Rasterizing shapefiles...')
#		adminShape_quarters = quarters(adminShape_loc, q, pixel_size, "shp")
#		adminShape_q=adminShape_quarters[0]
#		adminShape_q_dim=adminShape_quarters[1]
#		print(str(adminShape_q))
#		print(str(adminShape_q_dim))

		##### TEMPORARY RASTER FOR SABAH ONLY
		adminShape_q = '/home/mairin/deithvader/Documents/UBC/MScThesis/ResistanceData/MSY_Borneo.shp'
		sabah_ds = ogr.Open(adminShape_q,0)
		sabah_layer = sabah_ds.GetLayer()
		sabah_extent = sabah_layer.GetExtent()

### FOR SABAH ONLY
		x_ll = sabah_extent[0]
		y_ll = sabah_extent[2]
		y_ul = sabah_extent[3]
### WHEN RUNNING FULL CODE:
	# The minimum x, min y, max y
#		x_ll = adminShape_q_dim[0]
#		y_ll = adminShape_q_dim[1]
#		x_ul = adminShape_q_dim[2]
#		y_ul = adminShape_q_dim[3]

		# Procedure to create a geospatially explicit Tiff output (from http://gis.stackexchange.com/questions/37238/writing-numpy-array-to-raster-file)
		geotransform=(x_ll,pixel_size,0,y_ul,0, -pixel_size)

		# That's (top left x, w-e pixel resolution, rotation (0 if North is up),
		#         top left y, rotation (0 if North is up), n-s pixel resolution)

#		adminShape_q = str(out_path + '/MSY_Borneo_%s.shp' %(i))
		### Field for minor roads is "Country"
#		minRoad_loc = rasterize(quarters(minRoadPoly_loc_o, q, pixel_size, "shp"), 5, pixel_size)
		### Field for major roads is "TYPE"
#		majRoad_loc = rasterize(quarters(majRoadPoly_loc_o, q, pixel_size, "shp"), 0, pixel_size)
		mergedRoad_loc = cropRaster(mergedRoad_loc_o, adminShape_q)
#		landuse_loc = cropRaster(landuse_loc_o, adminShape_q)
		landuse_loc = cropRaster(landuse_loc_o, adminShape_q)
		strahler_loc = cropRaster(strahler_loc_o, adminShape_q)
#		strahler_loc = str(out_path + '/DEM_patch_t100_strahler_%s.tif' %(i))
		slope_loc = cropRaster(slope_loc_o, adminShape_q)
		print('...Finished raster modification, converting to numpy...')

		print 'Road location:%s\nLanduse location:%s\nStrahler location:%s\nSlope location: %s\n' %(mergedRoad_loc, landuse_loc, strahler_loc, slope_loc)
		# Load rasters

#		roadRaster = gdal.Open(mergedRoad_loc)
#		road_array = np.array(roadRaster.GetRasterBand(1).ReadAsArray())
#		na_road = roadRaster.GetRasterBand(1).GetNoDataValue()
#		roadRaster = None

		landuseRaster = gdal.Open(landuse_loc)
		landuse_o = np.array(landuseRaster.GetRasterBand(1).ReadAsArray())
		landuse_speed = np.empty(landuse_o.shape, dtype=float)

		na_landuse = landuseRaster.GetRasterBand(1).GetNoDataValue()
		landuseRaster = None

#		landuse_speed = landuse_array
#		road_array = None

#		sys.exit("Memory test") ####################

		print("~~~ Populating speed map")
		print("~~~ Quadrant %s") %q

		slopeRaster = gdal.Open(slope_loc)

		print "Map dimensions"
		print "Landuse array: " + str(landuse_speed.shape)
#		print "Slope array: " + str(slope_array.shape)
#		print "River array: " + str(river_array.shape)
#		print "Road array: " + str(road_array.shape)

		slope_array = np.array(slopeRaster.GetRasterBand(1).ReadAsArray())
		na_slope = slopeRaster.GetRasterBand(1).GetNoDataValue()

		# Drop the heavy file:
		slopeRaster = None
		print "Slope array: " + str(slope_array.shape)

### NOTE: Slope generated by GDAL's gdaldem function, slope expressed as % (but as whole percent, not decimal)
### Slope function stems from Tobler's walking velocity function: W = 6e^(-3.5*(dH/dX + 0.05))
		print('~~~~~~~~Modifying landuse with elevation')
		for i in range(0, landuse_speed.shape[0]):
			# Visual ticker of progress in %
			percent_comp = ((i*1.0)/landuse_speed.shape[0])*100.0
			if(percent_comp%2==0):
				print".",
			if(percent_comp%10==0):
				print str(percent_comp)+str("%"),
			for j in range(0, landuse_speed.shape[1]):
				if(landuse_o[i,j]==na_landuse or landuse_o[i,j]==5 or slope_array[i,j]):
					landuse_speed[i,j] = -999
					continue
						# Skip the NAs and cloud covered DEM areas, maximize at 200 as slope
				if(slope_array[i,j]==-9999 or slope_array[i,j]==0 or slope_array[i,j] > 200):
					if(landuse_o[i,j]==1):
						landuse_speed[i,j] = 1.0
					elif(landuse_o[i,j]==2):
						landuse_speed[i,j] = 1.25
					elif(landuse_o[i,j]==3):
						landuse_speed[i,j] = 2.0
					elif(landuse_o[i,j]==4):
						landuse_speed[i,j] = 30.0
					else:
						landuse_speed[i,j] = -999
					print 'Empty array %s, %s' %(i,j)
				elif(slope_array[i,j]<1000): # There are weird cloud data still in the slope dataset
					modifier = float(math.exp(-3.5*((slope_array[i,j]/100.0)+0.05)))
					if modifier==0.0:
						print 'ZERO ALERT! Modifier = 0, slope %s' %slope_array[i,j]

					if(landuse_o[i,j])==1:
						landuse_speed[i,j] = float(1.0*modifier)
					elif(landuse_o[i,j]==2):
						landuse_speed[i,j] = float(1.25*modifier)
					elif(landuse_o[i,j]==3):
						landuse_speed[i,j] = float(2.0*modifier)
						if(i==0 and j==17):
							print 'LS is:' + str(landuse_speed[i,j])
							print str(modifier)
							print str(landuse_o[i,j])
							print str(float(2.0*modifier))
					elif(landuse_o[i,j]==4):
						landuse_speed[i,j] = 30.0
					else:
						landuse_speed[i,j] = 15.0
					if(float(landuse_speed[i,j])<0):
						print 'NEGATIVE ALERT! Modifier: %s, Base: %s (is readable? %s)' %(modifier, landuse_o[i,j], landuse_o[i,j]==3)
						print 'ZERO ALERT! Modifier * base = %s' %(2.0*modifier)
						sys.exit('%s, %s\nSlope: %s\nLanduse:%s' %(i,j,slope_array[i,j], landuse_o[i,j]))
		# Drop the heavy file:
		slope_array = None

		print 'Creating walking speed output'
		speed_map_out = '/home/mairin/deithvader/Documents/UBC/MScThesis/ResistanceData/EditedMapFiles/Testing/WalkingSlopes'

#		geotransform=

		output_raster = gdal.GetDriverByName('GTiff').Create(str(speed_map_out+".tif"),landuse_speed.shape[1], landuse_speed.shape[0],1,gdal.GDT_Float32)  # Open the file
		output_raster.SetGeoTransform(geotransform) # Specify its coordinates
		srs = osr.SpatialReference()              	# Establish its coordinate encoding
		srs.ImportFromEPSG(4326)                    # This one specifies WGS84 lat long.

		output_raster.SetProjection(srs.ExportToWkt())   # Exports the coordinate system
		                                                   # to the file
		output_raster.GetRasterBand(1).WriteArray(landuse_speed)   # Writes my array to the raster

		# Modify landscape with rivers generated from DEM
		riverRaster = gdal.Open(strahler_loc)
		river_array = np.array(riverRaster.GetRasterBand(1).ReadAsArray())
		na_river = riverRaster.GetRasterBand(1).GetNoDataValue()
		riverRaster = None

		print('~~~~~~~~Adding rivers')
		for i in range(0, landuse_speed.shape[0]):
			percent_comp = (i*1.0)/landuse_speed.shape[0]*100.0
			if(percent_comp%2==0):
				print".",
			if(percent_comp%10==0):
				print str(percent_comp)+str("%")

			for j in range(0, landuse_speed.shape[1]):
				if(river_array[i,j] == 0 or river_array[i,j]==na_river):
					continue
				elif(river_array[i,j] == 3 or river_array[i,j] == 4):
					landuse_speed[i,j] = -999
				elif(river_array[i,j] > 5 and river_array[i,j] < 20):
					landuse_speed[i,j] = 20
				else:
					continue

				# if(landuse_speed[i,j] < 0.0000001):
				#	print 'SMALL ALERT!' + str(i) + str(j) + '+' + str(landuse_speed[i,j])

		# Drop the heavy file
		river_array = None

		print 'Creating river speed output'
		speed_map_out = '/home/mairin/deithvader/Documents/UBC/MScThesis/ResistanceData/EditedMapFiles/Testing/RiverSpeeds'

		output_raster = gdal.GetDriverByName('GTiff').Create(str(speed_map_out+".tif"),landuse_speed.shape[1], landuse_speed.shape[0],1,gdal.GDT_Float32)  # Open the file
		output_raster.SetGeoTransform(geotransform) # Specify its coordinates
		srs = osr.SpatialReference()              	# Establish its coordinate encoding
		srs.ImportFromEPSG(4326)                    # This one specifies WGS84 lat long.

		output_raster.SetProjection(srs.ExportToWkt())   # Exports the coordinate system
		                                                   # to the file
		output_raster.GetRasterBand(1).WriteArray(landuse_speed)   # Writes my array to the raster


		# Final modification: overlay roads, these trump all other speeds
		roadRaster = gdal.Open(mergedRoad_loc)
		road_array = np.array(roadRaster.GetRasterBand(1).ReadAsArray())
		na_road = roadRaster.GetRasterBand(1).GetNoDataValue()

		# Drop the heavy file:
		roadRaster = None

		print('~~~~~~~~Overlaying roads')
		for i in range(0, landuse_speed.shape[0]):
			percent_comp = (i*1.0)/landuse_speed.shape[0]*100.0
			if(percent_comp%2==0):
				print".",
			if(percent_comp%10==0):
				print str(percent_comp)+str("%"),
			for j in range(0, landuse_speed.shape[1]):
				# Roads are the strongest speed modifier
				if(road_array[i,j]==na_road):
					continue
				elif(road_array[i,j]==1):
					# Major roads, 60km/h
					landuse_speed[i,j] = 60
					continue
				elif(road_array[i,j]==2):
					# Logging roads, 10km/h
					landuse_speed[i,j] = 10
				else:
					continue


		# Drop the heavy file:
		road_array = None

		print 'Converting speeds to resistances'

		resistance = np.empty(landuse_o.shape, dtype=float)

		for i in range(0, landuse_speed.shape[0]):
			for j in range(0, landuse_speed.shape[1]):
#				No values should be zero, negative, or greater than 60km/h
				if(landuse_speed[i,j]<=0 or landuse_speed[i,j]>60):
					resistance[i,j] = -999
				elif(landuse_speed[i,j]>=0.01 and landuse_speed[i,j]<=60):
					resistance[i,j] = 1/landuse_speed[i,j]

		print 'Creating final speed output'

		speed_map_out = outpath+'speed_map_out'
		resistance_map_out = outpath+'resistance_map_out'

# Procedure to create a geospatially explicit Tiff output (from http://gis.stackexchange.com/questions/37238/writing-numpy-array-to-raster-file)
		geotransform=(x_ll,pixel_size,0,y_ul,0, -pixel_size)
		# That's (top left x, w-e pixel resolution, rotation (0 if North is up),
		#         top left y, rotation (0 if North is up), n-s pixel resolution)


# Save raster for speed (ASCII):
		print("Saving final image for Q%s" %q)
		np.savetxt(str(speed_map_out+'.asc'), landuse_speed, delimiter=' ', newline='\n', comments='',
			header=("NCOLS "+str(landuse_speed.shape[1])+"\n"+
				"NROWS "+str(landuse_speed.shape[0])+"\n"+
				"XLL CORNER %s\nYLL CORNER %s\nCELLSIZE %s\nNODATA_value -999" %(x_ll, y_ll, pixel_size)))

# Save raster for resistance (ASCII):
		print("Saving final image for Q%s" %q)
		np.savetxt(str(resistance_map_out+'.asc'), resistance, delimiter=' ', newline='\n', comments='',
			header=("NCOLS "+str(resistance.shape[1])+"\n"+
				"NROWS "+str(resistance.shape[0])+"\n"+
				"XLL CORNER %s\nYLL CORNER %s\nCELLSIZE %s\nNODATA_value -999" %(x_ll, y_ll, pixel_size)))

if __name__ == "__main__":
    main()

#### Old unused functions:

#def cropResample(source_raster, template_raster, method, polygon, p_size):
#	if method != "average":
#		if method !="mode":
#			sys.exit("Resample method %s not recognized. Please input either 'mode' or 'average'.") %(str(method))
#
#	(src_path, src_name) = os.path.split(source_raster)
#	(tmp_path, tmp_name) = os.path.split(template_raster)
#
#	out_name = str(src_name) + "_cropped.tif"
#	print "Cropping raster " + str(src_name) + " to adminShape"
#
#	os.system('gdalwarp -cutline ' + str(polygon) + ' -crop_to_cutline ' + str(source_raster) + ' ' + str(src_path) + '/' + str(out_name))
#	# croppedFile = str(src_path + '/' + src_name + "_cropped.tif")
#	print "Resampling raster " + str(src_name) + " to match " + str(tmp_name)
#
#	out_name = str(src_name) + "_rs.tif"
#	os.system('gdalwarp -t_srs EPSG:32601 -te 108.9998611 -0.0001389 120.0001389 8.0001389 -tr %s %s -r %s %s %s/%s' %(p_size, p_size, method, source_raster, out_path, out_name))
#
#	return str(out_path) + '/' + str(out_name)
#
#def cropRaster(base_raster, polygon):
#	# Takes in the paths of the raster to be modified and the shapefile to be used as a cropping template
#	(in_path, in_name) = os.path.split(base_raster)
#	short_in_name = in_name[0:len(in_name)-4]
#	out_name = str(short_in_name + '_cropped.tif')
#	os.system('gdalwarp -cutline %s -crop_to_cutline %s %s/%s -overwrite' %(polygon, base_raster, out_path, out_name))
#	return str(out_path + '/' + out_name)
