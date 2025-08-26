#!/bin/bash
mkdir -p h5_files_pools_1_2
cd h5_files_pools_1_2

# Mapped data directory
#path0="/mnt/SCDC/Bumblebee/2023_NR4A1_Colon_Training_Data/Mapped_Data-2024-genome-update"
path0="/x/2023_NR4A1_Colon_Training_Data/Mapped_Data-2024-genome-update/"


# Custom pools
pools=("NR4A1-Pool1-2024"  "NR4A1-Pool2-2024" )

echo "Main path : ${path0}"
echo "Pools : ${pools[@]}"

# Expected subdirectory for pulling sample counts
subdir0="outs/per_sample_outs"

for name in "${pools[@]}" ; do 
	echo "Current path $name"
	#mkdir -p ${name}
	# Creating the whole path
	pathi="${path0}/${name}/${subdir0}"
	# Getting subdirectories in the current variable
	str=`ls ${pathi}`
	subdirs=($str) 
	for namesub in "${subdirs[@]}"; do
		echo "      ${namesub}"
		# Copying sample_filtered_features_bc_matrix to target
		#path1="${pathi}/${namesub}/count/sample_filtered_feature_bc_matrix"
		path1="${pathi}/${namesub}/count/sample_filtered_feature_bc_matrix.h5"
		# Target directory

		#target="./${name}/${namesub}"
		target="./${namesub}"
		mkdir -p ${target}
		cp -rv ${path1} ${target}/.
	done
done
cd ../
