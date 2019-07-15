#!/bin/bash

#flirt -interp nearestneighbour -in cc200_roi_atlas.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz  -out cc200_roi_atlas_2mm.nii.gz

HCPDirectory=/flush2/andek67/Data/HCP/FMRI/REST

MaximumThreads=10
threads=0

for i in ${HCPDirectory}/* ; do

    # Go to current directory
	cd $i
	# Get subject name
   	Subject=${PWD##*/}
    echo "Processing" $Subject
	# Go back to original directory
	cd $HCPDirectory

	3dROIstats -mask cc200_roi_atlas_2mm.nii.gz ${i}/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_hp2000_clean.nii.gz > ${HCPDirectory}/HCP_${Subject}_rois_cc200.1D &

    ((threads++))

    if [ $threads -eq "$MaximumThreads" ]; then
		wait
		threads=0
	fi

done


