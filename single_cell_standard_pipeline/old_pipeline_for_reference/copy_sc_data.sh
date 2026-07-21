#!/bin/bash

mkdir -p copy_sc_data
cd copy_sc_data
# Mapped data directory
path0="/mnt/matrix/tigger/2026_Safe_Study17_Sc/mapped-data"
# Custom pools
pools=("Safe_Pool53" "Safe_Pool54" "Safe_Pool55" "Safe_Pool56" )

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
		path1="${pathi}/${namesub}/sample_filtered_feature_bc_matrix.h5"
		path2="${pathi}/${namesub}/sample_raw_probe_bc_matrix.h5"
		# Target directory

		#target="./${name}/${namesub}"
		target="./${namesub}"
		mkdir -p ${target}
		cp -rv ${path1} ${target}/.
		cp -rv ${path2} ${target}/.
	done
done
cd ../

mv copy_sc_data h5_files
