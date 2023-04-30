#!/bin/bash

MaximumThreads=1 # Maximum number of parallel downloads to use

for site_ in 1 2; do

    if [ "$site_" -eq "1" ]; then
        site=NYU
    elif [ "$site_" -eq "2" ]; then
        site=UM_1
    fi

    for pipeline_ in 1 2 3 4; do
   
    	threads=0

	if [ "$pipeline_" -eq "1" ]; then
	    pipeline=ccs
	elif [ "$pipeline_" -eq "2" ]; then
	    pipeline=cpac
	elif [ "$pipeline_" -eq "3" ]; then
      	    pipeline=dparsf
	elif [ "$pipeline_" -eq "4" ]; then
	    pipeline=niak
	fi

	for preprocessing_ in 1 2 3 4; do

	    if [ "$preprocessing_" -eq "1" ]; then
	        preprocessing=filt_global
	    elif [ "$preprocessing_" -eq "2" ]; then
	        preprocessing=filt_noglobal
	    elif [ "$preprocessing_" -eq "3" ]; then
	        preprocessing=nofilt_global
	    elif [ "$preprocessing_" -eq "4" ]; then
	        preprocessing=nofilt_noglobal
	    fi

	    # Only downloads subjects with mean framewise displacement of less than 0.2
	    python2.7 download_abide_preproc.py -d rois_cc200 -p ${pipeline} -s ${preprocessing} -o rois_cc200 -t ${site} &

            # Download all subjects
	    #python2.7 download_abide_preproc_allSubjects.py -d rois_cc200 -p ${pipeline} -s ${preprocessing} -o rois_cc200_allSubjects -t ${site} &

	    ((threads++))
	
	    if [ $threads -eq "$MaximumThreads" ]; then
	         wait
	         threads=0
	    fi

        done
    done
done
