#!/bin/sh

## Subject identifiers - 8/7/2018 n= 368
bblid_path=/data/jux/BBL/projects/ballerDepHeterogen/results/csvs/subset_with_T1_NbackFC_and_dem_with_clusters_justBBLID.csv
scanid_path=$(cat /data/jux/BBL/projects/ballerDepHeterogen/results/csvs/subset_with_T1_NbackFC_and_dem_with_clusters_justSCANID.csv)

#scanid_path=8465
## Parcellation and Template Paths- these parcels are the NEW ones, using the same parcels as ted's jneurosci paper but with NEW data
parcels_path=/data/jux/BBL/projects/ballerDepHeterogen/data/neuroimaging/nback/parcels/frac2back_n951_parametric_27roi_thr20.nii.gz
pnc_2mm_path=/data/joy/BBL/studies/pnc/template/pnc_template_brain_2mm.nii.gz
parcels_2mm_path=/data/jux/BBL/projects/ballerDepHeterogen/data/neuroimaging/nback/parcels/frac2back_n951_parametric_27roi_thr20_2mm.nii.gz
affine_mat=/data/joy/BBL/studies/pnc/template/mni2pnc0GenericAffine.mat
mni2pnc_warp=/data/joy/BBL/studies/pnc/template/mni2pnc1Warp.nii.gz
atlasname="jneurosci_labels"
parcels_just_names=$(cat /data/jux/BBL/projects/ballerDepHeterogen/data/neuroimaging/nback/parcels/roi27_just_names.txt)
parcels_just_numbers=/data/jux/BBL/projects/ballerDepHeterogen/data/neuroimaging/nback/parcels/roi27_just_numbers.txt

#output paths
outdir=/data/jux/BBL/projects/ballerDepHeterogen/data/neuroimaging/nback/parcels/roi_vals
output_temp=${outdir}/output_temp
nback_parcel_output_csv=/data/jux/BBL/projects/ballerDepHeterogen/data/neuroimaging/nback/parcels/mean_activation_by_27roi_by_SCANID.csv

#output_path=/data/jux/BBL/projects/ballerDepHeterogen/data/neuroimaging/nback/parcels/quantify_atlas_output
#output_temp=/data/jux/BBL/projects/ballerDepHeterogen/data/neuroimaging/nback/parcels/quantify_atlas_output/output_temp.csv

#################################
## Create 2mm PNC_path ##
#################################
#antsApplyTransforms -e 3 -d 3 -i ${parcels_path} -r ${pnc_2mm_path} -o ${parcels_2mm_path} -n MultiLabel
echo "applying antsApplyTransforms"
antsApplyTransforms -e 3 -d 3 -i ${parcels_path} -o ${parcels_2mm_path} -r ${pnc_2mm_path} -t ${affine_mat} -t ${mni2pnc_warp} -n MultiLabel
echo "transform complete"

## Define output directory for ROI values
mkdir -p ${outdir}
nback_paths=${outdir}/n368_2back0backStd_paths.txt

#remove old files and make room for new ones!
rm ${nback_paths}
touch ${nback_paths}

rm ${nback_parcel_output_csv}
touch ${nback_parcel_output_csv}

#set up the headers of the output
echo "scanid	"$parcels_just_names > ${nback_parcel_output_csv}
#############################################################################
## Extract 2back-0back ROI activation in Parcel from New Jneurosci parcels ##
#############################################################################
counter=1
for scanid in ${scanid_path}; do
	echo $scanid
	echo $counter
	# Define path to 2back-0back activation map in PNC template space
	img=$(ls /data/jux/BBL/projects/ballerDepHeterogen/data/neuroimaging/nback/nbackGlmBlockDesign/n1601_voxelwiseMaps_cope/*"${scanid}"_contrast4_2back0backStd.nii.gz) 
	echo "nback path: "${img}
	echo $img >> ${nback_paths}
	
	#remove output_temp and make new one
	rm ${output_temp}
	touch ${output_temp}
	
	#3droistats to get the mean per subject
	3dROIstats -mask ${parcels_2mm_path} ${img} > ${output_temp}

	#if it is the first run through, set the output file's first line to have the first line of output temp
	#if (($counter == 1))
	#then
	#	title=`more ${output_temp} | head -1`
	#	echo $title"	scanid" > ${nback_parcel_output_csv}	
		#more ${output_temp} | head -1 > $nback_parcel_output_csv
	#fi
		
	#pull out the last line (which has all of the mean data) and add to running file
	to_add=`more ${output_temp} | tail -1 | awk '{$1=$2=""; print $0}'`
	echo ${scanid}"	"$to_add >> ${nback_parcel_output_csv}
	#more ${output_temp} | tail -1 >> $nback_parcel_output_csv
	
	#extract mean values for subject within each roi and store it in file temp.csv
	#$XCPEDIR/utils/quantifyAtlas -v ${img} -s mean -o ${output_path}/${scanid} -a ${parcels_2mm_path} -n ${atlasname} -r ${parcels_just_numbers} -p ${scanid} -w 0
	#rm ${output_path}/*.1D  

	#append info in temp.csv into master output
	#echo "cat temp.csv into "$nback_parcel_output_csv
	#cat temp.csv >> ${nback_parcel_output_csv}
	((counter++))
	
done

#the next part is to combine everything, and relabel the numbers 1-27 with parcels_just_names

#run script
#${XCPEDIR}/utils/quantifyAtlas -v $inputfile  -s mean  -o ${ouput directory with subject indfentifiers}.csv  -a atlas.nii.gz -n atlasname  -r  atlasregioname.txt(orcsv) -p id0,id1 


