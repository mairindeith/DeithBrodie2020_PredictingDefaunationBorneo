#! /bin/sh
# My first script, largely practice for batch processing GRASS

# TEST
g.version
COUNTER=1
THRESH=100
# cropfile=sabah_crop.shp

for file in $(ls DEM_clipped_1*); do
	echo COUNTER: $COUNTER
###	echo "Clipping total map with $file"
### 
###	gdalwarp -cutline "$file" -crop_to_cutline -of GTiff -dstnodata -999 ASTER_DEM_Mosaic_cropped.tif DEM_clipped_"$COUNTER".tif -co COMPRESS=LZW -co TILED=YES --config GDAL_CACHEMAX 2048 -multi
###
###	echo "Finished clipping with file $COUNTER"
	
#	echo "Reading in $file"
	shortname=`echo "${file%%.*}"`
	echo "Shortname $shortname created"
######
###	r.in.gdal -e input=$shortname output=$shortfile location="$shortname"_f --o
######	
######	# echo Input $file
	echo Extracting streams from $shortname, threshold=$THRESH

	r.stream.extract --overwrite elevation="$shortname" threshold=$THRESH stream_raster="$shortname"_streamid stream_vector="$shortname"_streamid direction="$shortname"_flowdir
	# Note: currently not attempting to limit memory

	echo Calculating stream order
	r.stream.order -m memory=2000 stream_rast="$shortname"_streamid direction="$shortname"_flowdir strahler="$shortname"_strahler.tif shreve="$shortname"_shreve.tif
	echo "Finished stream order, exporting $shortname to TIFF"

	r.out.gdal input=${shortname_strahler} 	output="$shorname"_strahler.tif format=GTiff
	r.out.gdal input=${shortname_shreve} 	output="$shorname"_shreve.tif format=GTiff
	
	COUNTER="$COUNTER"1

done

### run GRASS' cleanup
# $GISBASE/etc/clean_temp
# rm -rf /tmp/grass7-$USER-$GIS_LOCK
# 
# unset GISBASE
# unset GISRC
#### 1) PREPARATION
# First we generate a script which contains the command(s) to be executed:
# for convenience, we save the file in our HOME directory
## You may use also a text editor for this, here we use the "echo" shell command
#
#echo "export GRASS_MESSAGE_FORMAT=plain
## Extract streams from elevation tif
#r.stream.extract elevation=$file threshold=5 
#v.out.ogr input=mymap3000 output=mymap3000.shp"  # > $HOME/my_grassjob.sh
#
## verify the content of the file
#cat $HOME/my_grassjob.sh
#
## make it user executable (this is important, use 'chmod' or via file manager)
#chmod u+x $HOME/my_grassjob.sh
#
## create a directory (may be elsewhere) to hold the location used for processing
#mkdir -p $HOME/grassdata
#
## create new temporary location for the job, exit after creation of this location
#grass70 -c epsg:32632 $HOME/grassdata/mytemploc_utm32n -e
#
##### 2) USING THE BATCH JOB
## define job file as environmental variable
#export GRASS_BATCH_JOB="$HOME/my_grassjob.sh"
#
## now we can use this new location and run the job defined via GRASS_BATCH_JOB
#grass70 $HOME/grassdata/mytemploc_utm32n/PERMANENT
#
##### 3) CLEANUP
## switch back to interactive mode, for the next GRASS GIS session
#unset GRASS_BATCH_JOB
#
## delete temporary location (consider to export results first in your batch job)
#rm -rf $HOME/grassdata/mytemploc_utm32n
#
## Now you can use the resulting SHAPE file "mymap3000.shp" elsewhere.
