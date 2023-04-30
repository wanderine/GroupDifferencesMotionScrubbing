#!/bin/bash

# Transform cc200 atlas to MNI, only needed once
#flirt -interp nearestneighbour -in cc200_roi_atlas.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz  -out cc200_roi_atlas_2mm.nii.gz

#HCPDirectory=/flush2/andek67/Data/HCP/FMRI/REST/ICAFIX
HCPDirectory=/flush2/andek67/Data/HCP/FMRI/REST/STANDARD

MaximumThreads=20
threads=0

# Loop over subjects
for i in ${HCPDirectory}/* ; do

    # Go to current directory
	cd $i
	# Get subject name
   	Subject=${PWD##*/}
    echo "Processing" $Subject
	# Go back to original directory
	cd $HCPDirectory

  # Get mean timeseries of all ROIs, ICAFIX
	#3dROIstats -mask cc200_roi_atlas_2mm.nii.gz ${i}/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_hp2000_clean.nii.gz > ${HCPDirectory}/HCP_${Subject}_rois_cc200.1D &

   # Get mean timeseries of all ROIs, STANDARD
    3dROIstats -mask cc200_roi_atlas_2mm.nii.gz ${i}/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR.nii.gz > ${HCPDirectory}/HCP_${Subject}_rois_cc200.1D &

    ((threads++))

    if [ $threads -eq "$MaximumThreads" ]; then
		wait
		threads=0
	fi

done
