#!/bin/bash 
set -e

# Requirements for this script
#  installed versions of: FSL5.0.1 or higher 
#  environment: FSLDIR

#Ely version 7/16/14 to generate gradient- and EPI-distortion corrected BOLD in native T1w space (rather than MNI space)

################################################ SUPPORT FUNCTIONS ##################################################

Usage() {
  echo "`basename $0`: Script to combine warps and affine transforms together and do a single resampling, with specified output resolution"
  echo " "
  echo "Usage: `basename $0` --workingdir=<working dir>"
  echo "             --infmri=<input fMRI 4D image>"
  echo "             --t1=<input native T1w restored image>"
  echo "             --fmriresout=<output resolution for images, typically the fmri resolution>"
  echo "             --fmrifolder=<fMRI processing folder>"
  echo "             --fmri2structin=<input fMRI to T1w warp>"
  echo "             --owarp=<output fMRI to T1w warp>"
  echo "             --oiwarp=<output T1w to fMRI warp>"
  echo "             --motionmatdir=<input motion correcton matrix directory>"
  echo "             --motionmatprefix=<input motion correcton matrix filename prefix>"
  echo "             --ofmri=<output fMRI 4D image>"
  echo "             --freesurferbrainmask=<input FreeSurfer brain mask, nifti format in native T1w space>"
  echo "             --biasfield=<input biasfield image, in native T1w space>"
  echo "             --gdfield=<input warpfield for gradient non-linearity correction>"
  echo "             --scoutin=<input scout image (EPI pre-sat, before gradient non-linearity distortion correction)>"
  echo "             --scoutgdcin=<input scout gradient nonlinearity distortion corrected image (EPI pre-sat)>"
  echo "             --oscout=<output transformed + distortion corrected scout image>"
  echo "             --jacobianin=<input Jacobian image>"
  echo "             --ojacobian=<output transformed + distortion corrected Jacobian image>"
}

# function for parsing options
getopt1() {
    sopt="$1"
    shift 1
    for fn in $@ ; do
	if [ `echo $fn | grep -- "^${sopt}=" | wc -w` -gt 0 ] ; then
	    echo $fn | sed "s/^${sopt}=//"
	    return 0
	fi
    done
}

defaultopt() {
    echo $1
}

################################################### OUTPUT FILES #####################################################

# Outputs (in $WD): 
#         NB: all these images are in native T1w space 
#             but at the specified resolution (to match the fMRI - i.e. low-res)
#     ${T1wImageFile}.${FinalfMRIResolution}  
#     ${FreeSurferBrainMaskFile}.${FinalfMRIResolution}
#     ${BiasFieldFile}.${FinalfMRIResolution}  
#     Scout_gdc_native_warp     : a warpfield from original (distorted) scout to low-res native T1w space
#
# Outputs (not in either of the above):
#     ${OutputTransform}  : the warpfield from fMRI to native T1w (low-res)
#     ${OutputfMRI}       
#     ${JacobianOut}
#     ${ScoutOutput}
#          NB: last three images are all in low-res native T1w space

################################################## OPTION PARSING #####################################################

# Just give usage if no arguments specified
if [ $# -eq 0 ] ; then Usage; exit 0; fi
# check for correct options

# parse arguments
WD=`getopt1 "--workingdir" $@`  # "$1"
InputfMRI=`getopt1 "--infmri" $@`  # "$2"
T1wImage=`getopt1 "--t1" $@`  # "$3"
FinalfMRIResolution=`getopt1 "--fmriresout" $@`  # "$4"
fMRIFolder=`getopt1 "--fmrifolder" $@`
fMRIToStructuralInput=`getopt1 "--fmri2structin" $@`  # "$6"
OutputTransform=`getopt1 "--owarp" $@`  # "$8"
OutputInvTransform=`getopt1 "--oiwarp" $@`
MotionMatrixFolder=`getopt1 "--motionmatdir" $@`  # "$9"
MotionMatrixPrefix=`getopt1 "--motionmatprefix" $@`  # "${10}"
OutputfMRI=`getopt1 "--ofmri" $@`  # "${11}"
FreeSurferBrainMask=`getopt1 "--freesurferbrainmask" $@`  # "${12}"
BiasField=`getopt1 "--biasfield" $@`  # "${13}"
GradientDistortionField=`getopt1 "--gdfield" $@`  # "${14}"
ScoutInput=`getopt1 "--scoutin" $@`  # "${15}"
ScoutInputgdc=`getopt1 "--scoutgdcin" $@`  # "${15}"
ScoutOutput=`getopt1 "--oscout" $@`  # "${16}"
JacobianIn=`getopt1 "--jacobianin" $@`  # "${17}"
JacobianOut=`getopt1 "--ojacobian" $@`  # "${18}"

BiasFieldFile=`basename "$BiasField"`
T1wImageFile=`basename $T1wImage`
FreeSurferBrainMaskFile=`basename "$FreeSurferBrainMask"`

echo " "
echo " START: OneStepResampling_Native"

mkdir -p $WD

# Record the input options in a log file
echo "$0 $@" >> $WD/log.txt
echo "PWD = `pwd`" >> $WD/log.txt
echo "date: `date`" >> $WD/log.txt
echo " " >> $WD/log.txt


########################################## DO WORK ########################################## 

#Save TR for later
TR_vol=`${FSLDIR}/bin/fslval ${InputfMRI} pixdim4 | cut -d " " -f 1`
NumFrames=`${FSLDIR}/bin/fslval ${InputfMRI} dim4`

#Ely tweak 9/5/14 to resume when aborted during OneStepResampling
if [ -e ${WD}/postvols ] ; then
	#Figure out the last volume processed and overwrite in case incomplete
	CompletedVols=(${WD}/postvols/*_mask.nii.gz)
	k=$((${#CompletedVols[@]}-1))
	echo "Resuming from volume $k"

	#Repopulate the strings that will be used to concatenate the processed volumes
	FrameMergeSTRING=""
	FrameMergeSTRINGII=""
	for ((j=0; j<$k; ++j)) ; do
		FrameMergeSTRING="${FrameMergeSTRING}${WD}/postvols/vol${j}.nii.gz " 
		FrameMergeSTRINGII="${FrameMergeSTRINGII}${WD}/postvols/vol${j}_mask.nii.gz " 
	done
else

# Create fMRI resolution native space files for T1w image, wmparc, and brain mask
${FSLDIR}/bin/flirt -interp spline -in ${T1wImage} -ref ${T1wImage} -applyisoxfm $FinalfMRIResolution -out ${WD}/${T1wImageFile}.${FinalfMRIResolution}
ResampRefIm=${WD}/${T1wImageFile}.${FinalfMRIResolution} 
${FSLDIR}/bin/applywarp --rel --interp=spline -i ${T1wImage} -r ${ResampRefIm} --premat=$FSLDIR/etc/flirtsch/ident.mat -o ${WD}/${T1wImageFile}.${FinalfMRIResolution}

# Create brain masks in this space from the FreeSurfer output (changing resolution)
${FSLDIR}/bin/applywarp --rel --interp=nn -i ${FreeSurferBrainMask}.nii.gz -r ${WD}/${T1wImageFile}.${FinalfMRIResolution} --premat=$FSLDIR/etc/flirtsch/ident.mat -o ${WD}/${FreeSurferBrainMaskFile}.${FinalfMRIResolution}.nii.gz

# Create versions of the biasfield (changing resolution)
${FSLDIR}/bin/applywarp --rel --interp=spline -i ${BiasField} -r ${WD}/${FreeSurferBrainMaskFile}.${FinalfMRIResolution}.nii.gz --premat=$FSLDIR/etc/flirtsch/ident.mat -o ${WD}/${BiasFieldFile}.${FinalfMRIResolution}
${FSLDIR}/bin/fslmaths ${WD}/${BiasFieldFile}.${FinalfMRIResolution} -thr 0.1 ${WD}/${BiasFieldFile}.${FinalfMRIResolution}

# Downsample warpfield (fMRI to native T1w) to increase speed 
#   NB: warpfield resolution is 10mm, so 1mm to fMRIres downsample loses no precision
${FSLDIR}/bin/convertwarp --relout --rel --warp1=${fMRIToStructuralInput} --ref=${WD}/${T1wImageFile}.${FinalfMRIResolution} --out=${OutputTransform}

# Stuff for RMS
invwarp -w ${OutputTransform} -o ${OutputInvTransform} -r ${ScoutInputgdc}
applywarp --rel --interp=nn -i ${FreeSurferBrainMask}.nii.gz -r ${ScoutInputgdc} -w ${OutputInvTransform} -o ${ScoutInputgdc}_mask_native.nii.gz
if [ -e ${fMRIFolder}/Movement_RelativeRMS_native.txt ] ; then
	/bin/rm -v ${fMRIFolder}/Movement_RelativeRMS_native.txt
fi
if [ -e ${fMRIFolder}/Movement_AbsoluteRMS_native.txt ] ; then
	/bin/rm -v ${fMRIFolder}/Movement_AbsoluteRMS_native.txt
fi
if [ -e ${fMRIFolder}/Movement_RelativeRMS_mean_native.txt ] ; then
	/bin/rm -v ${fMRIFolder}/Movement_RelativeRMS_mean_native.txt
fi
if [ -e ${fMRIFolder}/Movement_AbsoluteRMS_mean_native.txt ] ; then
	/bin/rm -v ${fMRIFolder}/Movement_AbsoluteRMS_mean_native.txt
fi

${FSLDIR}/bin/imcp ${WD}/${T1wImageFile}.${FinalfMRIResolution} ${fMRIFolder}/${T1wImageFile}.${FinalfMRIResolution}
${FSLDIR}/bin/imcp ${WD}/${FreeSurferBrainMaskFile}.${FinalfMRIResolution} ${fMRIFolder}/${FreeSurferBrainMaskFile}_native.${FinalfMRIResolution}
${FSLDIR}/bin/imcp ${WD}/${BiasFieldFile}.${FinalfMRIResolution} ${fMRIFolder}/${BiasFieldFile}.${FinalfMRIResolution}

mkdir -p ${WD}/prevols
mkdir -p ${WD}/postvols

# Apply combined transformations to fMRI (combines gradient non-linearity distortion, motion correction, and registration to T1w space, but keeping fMRI resolution)
${FSLDIR}/bin/fslsplit ${InputfMRI} ${WD}/prevols/vol -t
FrameMergeSTRING=""
FrameMergeSTRINGII=""
k=0
fi

while [ $k -lt $NumFrames ] ; do
  vnum=`${FSLDIR}/bin/zeropad $k 4`
  rmsdiff ${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum} ${MotionMatrixFolder}/${MotionMatrixPrefix}0000 ${ScoutInputgdc} ${ScoutInputgdc}_mask_native.nii.gz | tail -n 1 >> ${fMRIFolder}/Movement_AbsoluteRMS_native.txt
  prevmatrix="${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}" #Ely moved from below subsequent if statement as was causing error with resuming 9/5/14
  if [ $k -eq 0 ] ; then
    echo "0" >> ${fMRIFolder}/Movement_RelativeRMS_native.txt
  else
    rmsdiff ${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum} $prevmatrix ${ScoutInputgdc} ${ScoutInputgdc}_mask_native.nii.gz | tail -n 1 >> ${fMRIFolder}/Movement_RelativeRMS_native.txt
  fi
  ${FSLDIR}/bin/convertwarp --relout --rel --ref=${WD}/prevols/vol${vnum}.nii.gz --warp1=${GradientDistortionField} --postmat=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum} --out=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_native_gdc_warp.nii.gz
  ${FSLDIR}/bin/convertwarp --relout --rel --ref=${WD}/${T1wImageFile}.${FinalfMRIResolution} --warp1=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_native_gdc_warp.nii.gz --warp2=${OutputTransform} --out=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_native_all_warp.nii.gz
  ${FSLDIR}/bin/fslmaths ${WD}/prevols/vol${vnum}.nii.gz -mul 0 -add 1 ${WD}/prevols/vol${vnum}_mask.nii.gz
  ${FSLDIR}/bin/applywarp --rel --interp=spline --in=${WD}/prevols/vol${vnum}.nii.gz --warp=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_native_all_warp.nii.gz --ref=${WD}/${T1wImageFile}.${FinalfMRIResolution} --out=${WD}/postvols/vol${k}.nii.gz
  ${FSLDIR}/bin/applywarp --rel --interp=nn --in=${WD}/prevols/vol${vnum}_mask.nii.gz --warp=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_native_all_warp.nii.gz --ref=${WD}/${T1wImageFile}.${FinalfMRIResolution} --out=${WD}/postvols/vol${k}_mask.nii.gz
  FrameMergeSTRING="${FrameMergeSTRING}${WD}/postvols/vol${k}.nii.gz " 
  FrameMergeSTRINGII="${FrameMergeSTRINGII}${WD}/postvols/vol${k}_mask.nii.gz " 
  k=`echo "$k + 1" | bc`
  echo "Volume ${k} of ${NumFrames}complete" #Ely added 9/5/14
done

# Merge together results and restore the TR (saved beforehand)
${FSLDIR}/bin/fslmerge -tr ${OutputfMRI} $FrameMergeSTRING $TR_vol
${FSLDIR}/bin/fslmerge -tr ${OutputfMRI}_mask $FrameMergeSTRINGII $TR_vol
fslmaths ${OutputfMRI}_mask -Tmin ${OutputfMRI}_mask

# Combine transformations: gradient non-linearity distortion + fMRI_dc to native T1w space
${FSLDIR}/bin/convertwarp --relout --rel --ref=${WD}/${T1wImageFile}.${FinalfMRIResolution} --warp1=${GradientDistortionField} --warp2=${OutputTransform} --out=${WD}/Scout_gdc_native_warp.nii.gz
${FSLDIR}/bin/applywarp --rel --interp=spline --in=${ScoutInput} -w ${WD}/Scout_gdc_native_warp.nii.gz -r ${WD}/${T1wImageFile}.${FinalfMRIResolution} -o ${ScoutOutput}

# Create spline interpolated version of Jacobian  (native T1w space, fMRI resolution)
${FSLDIR}/bin/applywarp --rel --interp=spline -i ${JacobianIn} -r ${WD}/${T1wImageFile}.${FinalfMRIResolution} --premat=$FSLDIR/etc/flirtsch/ident.mat -o ${JacobianOut}

# Concatenate RMS parameters
cat ${fMRIFolder}/Movement_RelativeRMS_native.txt | awk '{ sum += $1} END { print sum / NR }' >> ${fMRIFolder}/Movement_RelativeRMS_mean_native.txt
cat ${fMRIFolder}/Movement_AbsoluteRMS_native.txt | awk '{ sum += $1} END { print sum / NR }' >> ${fMRIFolder}/Movement_AbsoluteRMS_mean_native.txt

echo " "
echo "END: OneStepResampling_Native"
echo " END: `date`" >> $WD/log.txt

########################################## QA STUFF ########################################## 

if [ -e $WD/qa.txt ] ; then rm -f $WD/qa.txt ; fi
echo "cd `pwd`" >> $WD/qa.txt
echo "# Check registrations to low-res standard space" >> $WD/qa.txt
echo "fslview ${WD}/${T1wImageFile}.${FinalfMRIResolution} ${WD}/${FreeSurferBrainMaskFile}.${FinalfMRIResolution} ${WD}/${BiasFieldFile}.${FinalfMRIResolution} ${OutputfMRI}" >> $WD/qa.txt

##############################################################################################


