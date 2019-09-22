#!/bin/sh

## Subject identifiers - 8/7/2018 n= 368
bblid_path=/data/jux/BBL/projects/ballerDepHeterogen/results/csvs/subset_with_T1_NbackFC_and_dem_with_clusters_justBBLID.csv
scanid_path=$(cat /data/jux/BBL/projects/ballerDepHeterogen/results/csvs/subset_with_T1_NbackFC_and_dem_with_clusters_justSCANID.csv)

## Parcellation and Template Paths
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


#################################
## Create 2mm PNC_path ##
#################################
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


	#pull out the last line (which has all of the mean data) and add to running file
	to_add=`more ${output_temp} | tail -1 | awk '{$1=$2=""; print $0}'`
	echo ${scanid}"	"$to_add >> ${nback_parcel_output_csv}
	
	((counter++))
	
done


