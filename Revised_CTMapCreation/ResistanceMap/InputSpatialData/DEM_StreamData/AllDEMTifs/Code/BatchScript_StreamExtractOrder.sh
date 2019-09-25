#! /bin/sh
# My first script, largely practice for batch processing GRASS

# TEST
g.version
COUNTER=0
THRESH=100
# cropfile=sabah_crop.shp

for file in $(g.list -r type=rast pattern=DEM_clipped*1); do
	echo COUNTER_$COUNTER
	echo "Processing $file"

	shortname=$(echo ${file##*/} | awk -F'[_.]' '{print $3}')_FULL
	# echo "Beginning $shortname processing"

	# gdalwarp -cutline $cropfile
	# r.in.gdal -e input=$file output=${file%.tif} location="$shortname"_run1 --o
	
	# echo Input $file
	echo Extracting streams from $file, threshold=$THRESH
	r.stream.extract --overwrite elevation=${file%.tif} threshold=$THRESH stream_raster="$shortname"_streamid stream_vector="$shortname"_streamid direction="$shortname"_flowdir
	# Note: currently not attempting to limit memory
	echo Calculating stream order
	r.stream.order -m memory=1200 stream_rast="$shortname"_streamid direction="$shortname"_flowdir strahler="$shortname"_strahler.tif
	echo "Finished stream order, exporting $shortname to TIFF"
#	r.out.gdal input=${shortname_strahler%} output="$shorname"_strahler.tif format=GTiff

done

file=DEM_clipped_0
shortname=DEM_clipped_0
echo Extracting streams from $file, threshold=$THRESH
	r.stream.extract --overwrite elevation=${file%.tif} threshold=$THRESH stream_raster="$shortname"_streamid stream_vector="$shortname"_streamid direction="$shortname"_flowdir
	# Note: currently not attempting to limit memory
	echo Calculating stream order
	r.stream.order -m memory=400 stream_rast="$shortname"_streamid direction="$shortname"_flowdir strahler="$shortname"_strahler.tif
	echo "Finished stream order, exporting $shortname to TIFF"
#	r.out.gdal input=${shortname_strahler%} output="$shorname"_strahler.tif format=GTiff

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