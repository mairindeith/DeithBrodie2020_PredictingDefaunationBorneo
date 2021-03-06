#!/bin/bash -x

# -----------------------
# Supplementary Material for Deith and Brodie 2020; “Predicting defaunation – accurately mapping bushmeat hunting pressure over large areas”
#    doi: 10.1098/rspb.2019-2677
#
# Code modified from Leonard et al's 2017 GFLOW: software for modelling circuit theory‐based connectivity at any scale. (doi: 10.1111/2041-210X.12689)
# -----------------------
#
# Gflow must be compiled locally before executing this script. Dependencies for GFlow include: openmpi, hypre, and petsc.
# If you would like to execute a random shuffle of all pairwise, please install 'coreutils'as well.

# Execution flags are commented below and demonstrate an example execution of GFlow. To execute script as is: type 'sh execute_example.sh'
# in Terminal.

which mpiexec

# Set and add PETSc to PATH (Please update this if you are using Linux or any installation proceedure that differs from the README)
export PETSC_DIR=/usr/lib/petscdir/3.7.6/x86_64-linux-gnu-real
export LD_LIBRARY_PATH=${PETSC_DIR}/lib:$LD_LIBRARY_PATH

# Set Output Directory: Default is Current Directory
OUTPUT_DIR=.

# DEBUG: Set number of random pairs to calculate from all possible pairs (Currently -n=5). Must be used with -node_pairs flag
# Allows exact number of test pairs. Remove or comment line for all pairwise. Requires 'coreutils' and 'all.tsv' from inputs.
	# gshuf all.tsv -n 5 > ${OUTPUT_DIR}/shuf.tsv

# Set the Clock
SECONDS=0
date

# REQUIRED flags to execute GFlow are below. Please read descriptions or see example arguments for use.
#	Set Number of processes (CPUs) and call Gflow. Currently -n = 4 below

	# -output_sum_density_filename
		# Set Output Path, file name, and format of cummulative current density.
		# For use see: https://github.com/Pbleonard/GFlow/issues/8
	# -ouput_prefix  (DEPRICATED)
		# Set Prefix for output files: Currently = 'local_'
	# -ouput_directory (DEPRICATED)
		# Set Output Directory: Set to default above
	# -habitat
		# Set Habitat Map or resistance surface (.asc)
	# -nodes
		# Set Focal Nodes or Source/Destination Points (.txt list of point pairs to calculate or .asc grid). Inputs must be points.
		# If using a .txt list, the point coorindates must be relative to the resistance surface grid. Please look at example inputs.
	# -output_format (DEPRICATED)
		# Set Desired Output formaat --'asc' or 'amps' -- Default = asc
	# -output_final_current_only (DEPRICATED)
		# Select output options. 1 = Only summation; 0 = Pairwise calculations + Summation

# OPTIONAL flags

	# -output_density_filename
		# Set Output Path, file name prefix, and format (i.e., *.asc, *.asc.gz) of pairwise calculations. Omitting this flag will
		# discard each pairwise solve output and assume you want the cummulative output only.
		# For use see: https://github.com/Pbleonard/GFlow/issues/8
	# -node_pairs
# Calculate only desired node pairs if input (e.g., '${OUTPUT_DIR}/shuf.tsv \' from gshuf above). Currently not used.
	# -converge_at
		# Set Convergence Factor to stop calculating. Typically used in place of 'node_pairs' or if all pairwise is too
		# computationlly time consuing. Acceptable formats include: '4N' or '.9999'. Set to '1N' Below.
	# -shuffle_node_pairs
		# Shuffles pairs for random selection. Input is binary. Currently set to shuffle below (= 1)
	# -effective_resistance
		# Print effective resistance to log file. Supply path for .csv

# Assigning Arguments to Flags for Execution:
final_output=cumulative_map.asc
info_textfile=../GFlowOutputs/05_03_2018/node_data.csv
echo Multiplier,nNodes,nIterations,InputFilename,OutputIntermediate > $info_textfile
loop=0
# for file in $(ls ../Nodes/NodesTXT/*.txt | tail -q | head -1); do
for file in $(ls ../Nodes/NodesTXT/*.txt); do
	# Extract the population density for that set of points, set as multiplier
	mult_int=$(echo $file | grep -o -E '[0-9]+' | sed -n '1p')
	mult_dec=$(echo $file | grep -o -E '[0-9]+' | sed -n '2p')
	count=$(echo $file | grep -o -E '[0-9]+' | sed -n '3p')
	multiplier=$(echo "$mult_int.$mult_dec")
	# echo $multiplier
	filename=nodes_${multiplier}_${count}
	echo ~~~~~~~~~~~~~~~~~~~~NODE: $multiplier, cost distance: $cost~~~~~~~~~~~~~~~~~~~~
	# Generate the file name without the path or extension
	mpiexec -n 8 GFlow-0.1.7-alpha/gflow.x \
		-habitat ../InputFiles/ButonResist_cs1.asc \
		-nodes ../Nodes/NodesTXT/$filename.txt \
		-node_pairs ../Nodes/NodesTSV/$filename.tsv \
		-effective_resistance ./R_eff.csv \
		-output_sum_density_filename "../GFlowOutputs/05_03_2018/Temp/output_temp_{iter}.asc" \
		-converge_at 3N \
		-shuffle_node_pairs 1 \

		: "walltime: $SECONDS seconds"

	# Notes:
	# These lines of code convert the output of GFlow, divide the current by the number of iterations run,
	# and then multiply it by the population modifier based on population density and the number of input nodes
	echo GFlow processing finished, applying map multipliers...
	output_file_full=$(basename $(ls ../GFlowOutputs/05_03_2018/Temp/output_temp*.asc))
	output_file=${output_file_full%.*}
	n_iterations=$(echo $output_file | grep -o -E '[0-9]+')
	sources_count=$(($(awk '{print $2}' ../Nodes/NodesTSV/$filename.tsv | head -n 1) - 1))
	log_percent=$(echo 'l('${multiplier}')/l(10)*-0.4536+1.0765' | bc -l)
	decimal_percent=$(echo 'e('$log_percent'*l(10))/100' | bc -l)
	adjusted_multiplier=$(echo $decimal_percent'*'$multiplier | bc -l)

	gdal_calc.py -A ../GFlowOutputs/05_03_2018/Temp/$output_file_full -B ../GFlowOutputs/05_03_2018/$final_output --outfile=../GFlowOutputs/05_03_2018/Temp/intermediate_cumulative.asc --calc="(A*$adjusted_multiplier*($sources_count/1)/$n_iterations)+B" --overwrite --NoDataValue=-9999 && mv ../GFlowOutputs/05_03_2018/Temp/$output_file_full ../GFlowOutputs/05_03_2018/IntermediateMaps/${output_file}_${multiplier}.asc && cp ../GFlowOutputs/05_03_2018/Temp/intermediate_cumulative.asc ../GFlowOutputs/05_03_2018/$final_output && mv ../Nodes/NodesTSV/$filename.tsv ../Nodes/NodesTSV/Processed/$filename.tsv && mv ../Nodes/NodesTXT/$filename.txt ../Nodes/NodesTXT/Processed/$filename.txt
	let "loop++"
done
echo Count: $loop
