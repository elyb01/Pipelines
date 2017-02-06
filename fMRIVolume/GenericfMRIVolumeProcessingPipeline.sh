#!/bin/bash 
set -e

# Requirements for this script
#  installed versions of: FSL (version 5.0.6), FreeSurfer (version 5.3.0-HCP) , gradunwarp (HCP version 1.0.2) 
#  environment: use SetUpHCPPipeline.sh  (or individually set FSLDIR, FREESURFER_HOME, HCPPIPEDIR, PATH - for gradient_unwarp.py)

########################################## PIPELINE OVERVIEW ########################################## 

# TODO

########################################## OUTPUT DIRECTORIES ########################################## 

# TODO

# --------------------------------------------------------------------------------
#  Load Function Libraries
# --------------------------------------------------------------------------------

source $HCPPIPEDIR/global/scripts/log.shlib  # Logging related functions
source $HCPPIPEDIR/global/scripts/opts.shlib # Command line option functions

################################################ SUPPORT FUNCTIONS ##################################################

# --------------------------------------------------------------------------------
#  Usage Description Function
# --------------------------------------------------------------------------------

show_usage() {
    echo "Usage information To Be Written"
    exit 1
}

# --------------------------------------------------------------------------------
#   Establish tool name for logging
# --------------------------------------------------------------------------------
log_SetToolName "GenericfMRIVolumeProcessingPipelinePlusNativeBOLD_AutoResume.sh"

################################################## OPTION PARSING #####################################################

opts_ShowVersionIfRequested $@

if opts_CheckForHelpRequest $@; then
    show_usage
fi

log_Msg "Parsing Command Line Options"

# parse arguments
Path=`opts_GetOpt1 "--path" $@`
log_Msg "Path: ${Path}"

Subject=`opts_GetOpt1 "--subject" $@`
log_Msg "Subject: ${Subject}"

NameOffMRI=`opts_GetOpt1 "--fmriname" $@`
log_Msg "NameOffMRI: ${NameOffMRI}"

fMRITimeSeries=`opts_GetOpt1 "--fmritcs" $@`
log_Msg "fMRITimeSeries: ${fMRITimeSeries}"

fMRIScout=`opts_GetOpt1 "--fmriscout" $@`
log_Msg "fMRIScout: ${fMRIScout}"

SpinEchoPhaseEncodeNegative=`opts_GetOpt1 "--SEPhaseNeg" $@`
log_Msg "SpinEchoPhaseEncodeNegative: ${SpinEchoPhaseEncodeNegative}"

SpinEchoPhaseEncodePositive=`opts_GetOpt1 "--SEPhasePos" $@`
log_Msg "SpinEchoPhaseEncodePositive: ${SpinEchoPhaseEncodePositive}"

MagnitudeInputName=`opts_GetOpt1 "--fmapmag" $@`  # Expects 4D volume with two 3D timepoints
log_Msg "MagnitudeInputName: ${MagnitudeInputName}"

PhaseInputName=`opts_GetOpt1 "--fmapphase" $@`  
log_Msg "PhaseInputName: ${PhaseInputName}"

GEB0InputName=`opts_GetOpt1 "--fmapgeneralelectric" $@`
log_Msg "GEB0InputName: ${GEB0InputName}"

DwellTime=`opts_GetOpt1 "--echospacing" $@`  
log_Msg "DwellTime: ${DwellTime}"

deltaTE=`opts_GetOpt1 "--echodiff" $@`  
log_Msg "deltaTE: ${deltaTE}"

UnwarpDir=`opts_GetOpt1 "--unwarpdir" $@`  
log_Msg "UnwarpDir: ${UnwarpDir}"

FinalfMRIResolution=`opts_GetOpt1 "--fmrires" $@`  
log_Msg "FinalfMRIResolution: ${FinalfMRIResolution}"

# FIELDMAP, SiemensFieldMap, GeneralElectricFieldMap, or TOPUP
# Note: FIELDMAP and SiemensFieldMap are equivalent
DistortionCorrection=`opts_GetOpt1 "--dcmethod" $@`
log_Msg "DistortionCorrection: ${DistortionCorrection}"

BiasCorrection=`opts_GetOpt1 "--biascorrection" $@`
# Convert BiasCorrection value to all UPPERCASE (to allow the user the flexibility to use NONE, None, none, legacy, Legacy, etc.)
BiasCorrection="$(echo ${BiasCorrection} | tr '[:lower:]' '[:upper:]')"
log_Msg "BiasCorrection: ${BiasCorrection}"

GradientDistortionCoeffs=`opts_GetOpt1 "--gdcoeffs" $@`  
log_Msg "GradientDistortionCoeffs: ${GradientDistortionCoeffs}"

TopupConfig=`opts_GetOpt1 "--topupconfig" $@`  # NONE if Topup is not being used
log_Msg "TopupConfig: ${TopupConfig}"

dof=`opts_GetOpt1 "--dof" $@`
dof=`opts_DefaultOpt $dof 6`
log_Msg "dof: ${dof}"

RUN=`opts_GetOpt1 "--printcom" $@`  # use ="echo" for just printing everything and not running the commands (default is to run)
log_Msg "RUN: ${RUN}"

#NOTE: the jacobian option only applies the jacobian of the distortion corrections to the fMRI data, and NOT from the nonlinear T1 to template registration
UseJacobian=`opts_GetOpt1 "--usejacobian" $@`
# Convert UseJacobian value to all lowercase (to allow the user the flexibility to use True, true, TRUE, False, False, false, etc.)
UseJacobian="$(echo ${UseJacobian} | tr '[:upper:]' '[:lower:]')"
log_Msg "UseJacobian: ${UseJacobian}"

MotionCorrectionType=`opts_GetOpt1 "--mctype" $@`  # use = "FLIRT" to run FLIRT-based mcflirt_acc.sh, or "MCFLIRT" to run MCFLIRT-based mcflirt.sh
MotionCorrectionType=`opts_DefaultOpt $MotionCorrectionType MCFLIRT` #use mcflirt by default

#error check
case "$MotionCorrectionType" in
    MCFLIRT|FLIRT)
        #nothing
    ;;
    
    *)
        log_Msg "ERROR: --mctype must be 'MCFLIRT' (default) or 'FLIRT'"
        exit 1
    ;;
esac

JacobianDefault="true"
if [[ $DistortionCorrection != "TOPUP" ]]
then
    #because the measured fieldmap can cause the warpfield to fold over, default to doing nothing about any jacobians
    JacobianDefault="false"
    #warn if the user specified it
    if [[ $UseJacobian == "true" ]]
    then
        log_Msg "WARNING: using --jacobian=true with --dcmethod other than TOPUP is not recommended, as the distortion warpfield is less stable than TOPUP"
    fi
fi
log_Msg "JacobianDefault: ${JacobianDefault}"

UseJacobian=`opts_DefaultOpt $UseJacobian $JacobianDefault`
log_Msg "After taking default value if necessary, UseJacobian: ${UseJacobian}"

if [[ -n $HCPPIPEDEBUG ]]
then
    set -x
fi

#sanity check the jacobian option
if [[ "$UseJacobian" != "true" && "$UseJacobian" != "false" ]]
then
    log_Msg "the --usejacobian option must be 'true' or 'false'"
    exit 1
fi

# Setup PATHS
PipelineScripts=${HCPPIPEDIR_fMRIVol}
GlobalScripts=${HCPPIPEDIR_Global}

#Naming Conventions
T1wImage="T1w_acpc_dc"
T1wRestoreImage="T1w_acpc_dc_restore"
T1wRestoreImageBrain="T1w_acpc_dc_restore_brain"
T1wFolder="T1w" #Location of T1w images
AtlasSpaceFolder="MNINonLinear"
ResultsFolder="Results"
BiasField="BiasField_acpc_dc"
BiasFieldMNI="BiasField"
T1wAtlasName="T1w_restore"
T1wNativeName="T1w_acpc_dc_restore"
MovementRegressor="Movement_Regressors" #No extension, .txt appended
MotionMatrixFolder="MotionMatrices"
MotionMatrixPrefix="MAT_"
FieldMapOutputName="FieldMap"
MagnitudeOutputName="Magnitude"
MagnitudeBrainOutputName="Magnitude_brain"
ScoutName="Scout"
OrigScoutName="${ScoutName}_orig"
OrigTCSName="${NameOffMRI}_orig"
FreeSurferBrainMask="brainmask_fs"
fMRI2strOutputTransform="${NameOffMRI}2str"
RegOutput="Scout2T1w"
AtlasTransform="acpc_dc2standard"
OutputfMRI2StandardTransform="${NameOffMRI}2standard"
Standard2OutputfMRITransform="standard2${NameOffMRI}"
OutputfMRI2NativeTransform="${NameOffMRI}2native" #Native
Native2OutputfMRITransform="native2${NameOffMRI}" #Native
QAImage="T1wMulEPI"
JacobianOut="Jacobian"
SubjectFolder="$Path"/"$Subject"
#note, this file doesn't exist yet, gets created by ComputeSpinEchoBiasField.sh during DistortionCorrectionAnd...
sebasedBiasFieldMNI="$SubjectFolder/$AtlasSpaceFolder/Results/$NameOffMRI/${NameOffMRI}_sebased_bias.nii.gz"

fMRIFolder="$Path"/"$Subject"/"$NameOffMRI"

#error check bias correction opt
case "$BiasCorrection" in
    NONE)
        UseBiasFieldMNI=""
    ;;
    LEGACY)
        UseBiasFieldMNI="${fMRIFolder}/${BiasFieldMNI}.${FinalfMRIResolution}"
    ;;
    
    SEBASED)
        echo "CAUTION: SEBASED BIAS CORRECTION NOT CURRENTLY SET UP FOR NATIVE FMRI PREPROCESSING!"
		if [[ "$DistortionCorrection" != "TOPUP" ]]
        then
            log_Msg "SEBASED bias correction is only available with --dcmethod=TOPUP"
            exit 1
        fi
        UseBiasFieldMNI="$sebasedBiasFieldMNI"
    ;;
    
    "")
        log_Msg "--biascorrection option not specified"
        exit 1
    ;;
    
    *)
        log_Msg "unrecognized value for bias correction: $BiasCorrection"
    exit 1
esac


########################################## DO WORK ########################################## 

T1wFolder="$Path"/"$Subject"/"$T1wFolder"
AtlasSpaceFolder="$Path"/"$Subject"/"$AtlasSpaceFolder"
ResultsFolder="$AtlasSpaceFolder"/"$ResultsFolder"/"$NameOffMRI"
NativeResultsFolder="$T1wFolder"/"$ResultsFolder"/"$NameOffMRI" #Native
fMRIFolder="$Path"/"$Subject"/"$NameOffMRI"

if [ ! -e "$fMRIFolder" ] ; then
  log_Msg "mkdir ${fMRIFolder}"
  mkdir "$fMRIFolder"
fi

if [ ! -e  "$fMRIFolder"/"$OrigTCSName".nii.gz ] ; then
	cp "$fMRITimeSeries" "$fMRIFolder"/"$OrigTCSName".nii.gz
fi

if [ ! -e "$fMRIFolder"/"$OrigScoutName".nii.gz ] ; then
	#Create fake "Scout" if it doesn't exist
	if [ $fMRIScout = "NONE" ] ; then
		${RUN} ${FSLDIR}/bin/fslroi "$fMRIFolder"/"$OrigTCSName" "$fMRIFolder"/"$OrigScoutName" 0 1
	else
		cp "$fMRIScout" "$fMRIFolder"/"$OrigScoutName".nii.gz
	fi
fi

if [ ! -e "${fMRIFolder}/MotionCorrection_FLIRTbased" ] ; then
	#Gradient Distortion Correction of fMRI
	log_Msg "Gradient Distortion Correction of fMRI"
	if [ ! $GradientDistortionCoeffs = "NONE" ] ; then
		log_Msg "mkdir -p ${fMRIFolder}/GradientDistortionUnwarp"
		mkdir -p "$fMRIFolder"/GradientDistortionUnwarp
    	${RUN} "$GlobalScripts"/GradientDistortionUnwarp.sh \
		--workingdir="$fMRIFolder"/GradientDistortionUnwarp \
		--coeffs="$GradientDistortionCoeffs" \
		--in="$fMRIFolder"/"$OrigTCSName" \
		--out="$fMRIFolder"/"$NameOffMRI"_gdc \
		--owarp="$fMRIFolder"/"$NameOffMRI"_gdc_warp
		
		log_Msg "mkdir -p ${fMRIFolder}/${ScoutName}_GradientDistortionUnwarp"	
     	mkdir -p "$fMRIFolder"/"$ScoutName"_GradientDistortionUnwarp
     	${RUN} "$GlobalScripts"/GradientDistortionUnwarp.sh \
	 	--workingdir="$fMRIFolder"/"$ScoutName"_GradientDistortionUnwarp \
		--coeffs="$GradientDistortionCoeffs" \
		--in="$fMRIFolder"/"$OrigScoutName" \
		--out="$fMRIFolder"/"$ScoutName"_gdc \
		--owarp="$fMRIFolder"/"$ScoutName"_gdc_warp
	else
		log_Msg "NOT PERFORMING GRADIENT DISTORTION CORRECTION"
		${RUN} ${FSLDIR}/bin/imcp "$fMRIFolder"/"$OrigTCSName" "$fMRIFolder"/"$NameOffMRI"_gdc
		${RUN} ${FSLDIR}/bin/fslroi "$fMRIFolder"/"$NameOffMRI"_gdc "$fMRIFolder"/"$NameOffMRI"_gdc_warp 0 3
		${RUN} ${FSLDIR}/bin/fslmaths "$fMRIFolder"/"$NameOffMRI"_gdc_warp -mul 0 "$fMRIFolder"/"$NameOffMRI"_gdc_warp
		${RUN} ${FSLDIR}/bin/imcp "$fMRIFolder"/"$OrigScoutName" "$fMRIFolder"/"$ScoutName"_gdc
	fi
else
	echo "\nGradient Distortion Correction appears complete. Skipping to Motion Correction.\n"
fi

if [ ! -e "${fMRIFolder}/DistortionCorrectionAndEPIToT1wReg_FLIRTBBRAndFreeSurferBBRbased" ] ; then
	log_Msg "mkdir -p ${fMRIFolder}/MotionCorrection_FLIRTbased"
	mkdir -p "$fMRIFolder"/MotionCorrection_FLIRTbased
	${RUN} "$PipelineScripts"/MotionCorrection_FLIRTbased.sh \
		"$fMRIFolder"/MotionCorrection_FLIRTbased \
		"$fMRIFolder"/"$NameOffMRI"_gdc \
    	"$fMRIFolder"/"$ScoutName"_gdc \
    	"$fMRIFolder"/"$NameOffMRI"_mc \
    	"$fMRIFolder"/"$MovementRegressor" \
    	"$fMRIFolder"/"$MotionMatrixFolder" \
    	"$MotionMatrixPrefix" 
else
	echo "\nMotion Correction appears complete. Skipping to EPI Correction.\n"
fi

if [ ! -e "${fMRIFolder}/OneStepResampling" ] ; then
	#EPI Distortion Correction and EPI to T1w Registration
	log_Msg "EPI Distortion Correction and EPI to T1w Registration"
	if [ -e ${fMRIFolder}/DistortionCorrectionAndEPIToT1wReg_FLIRTBBRAndFreeSurferBBRbased ] ; then
		rm -r ${fMRIFolder}/DistortionCorrectionAndEPIToT1wReg_FLIRTBBRAndFreeSurferBBRbased
	fi
	log_Msg "mkdir -p ${fMRIFolder}/DistortionCorrectionAndEPIToT1wReg_FLIRTBBRAndFreeSurferBBRbased"
	mkdir -p ${fMRIFolder}/DistortionCorrectionAndEPIToT1wReg_FLIRTBBRAndFreeSurferBBRbased
	
	${RUN} ${PipelineScripts}/DistortionCorrectionAndEPIToT1wReg_FLIRTBBRAndFreeSurferBBRbased.sh \
		--workingdir=${fMRIFolder}/DistortionCorrectionAndEPIToT1wReg_FLIRTBBRAndFreeSurferBBRbased \
    	--scoutin=${fMRIFolder}/${ScoutName}_gdc \
    	--t1=${T1wFolder}/${T1wImage} \
    	--t1restore=${T1wFolder}/${T1wRestoreImage} \
    	--t1brain=${T1wFolder}/${T1wRestoreImageBrain} \
    	--fmapmag=${MagnitudeInputName} \
    	--fmapphase=${PhaseInputName} \
    	--echodiff=${deltaTE} \
    	--SEPhaseNeg=${SpinEchoPhaseEncodeNegative} \
   		--SEPhasePos=${SpinEchoPhaseEncodePositive} \
    	--echospacing=${DwellTime} \
    	--unwarpdir=${UnwarpDir} \
    	--owarp=${T1wFolder}/xfms/${fMRI2strOutputTransform} \
    	--biasfield=${T1wFolder}/${BiasField} \
    	--oregim=${fMRIFolder}/${RegOutput} \
    	--freesurferfolder=${T1wFolder} \
    	--freesurfersubjectid=${Subject} \
    	--gdcoeffs=${GradientDistortionCoeffs} \
    	--qaimage=${fMRIFolder}/${QAImage} \
    	--method=${DistortionCorrection} \
    	--topupconfig=${TopupConfig} \
    	--ojacobian=${fMRIFolder}/${JacobianOut} 
else
	echo "\nEPI Correction appears complete. Skipping to One Step Resampling.\n"
fi

if [ ! -e "${fMRIFolder}/${NameOffMRI}_nonlin.nii.gz" ] ; then
	#One Step Resampling
	log_Msg "One Step Resampling"
	log_Msg "mkdir -p ${fMRIFolder}/OneStepResampling"
	mkdir -p ${fMRIFolder}/OneStepResampling
	${RUN} ${PipelineScripts}/OneStepResampling.sh \
		--workingdir=${fMRIFolder}/OneStepResampling \
    	--infmri=${fMRIFolder}/${OrigTCSName}.nii.gz \
    	--t1=${AtlasSpaceFolder}/${T1wAtlasName} \
    	--fmriresout=${FinalfMRIResolution} \
    	--fmrifolder=${fMRIFolder} \
    	--fmri2structin=${T1wFolder}/xfms/${fMRI2strOutputTransform} \
    	--struct2std=${AtlasSpaceFolder}/xfms/${AtlasTransform} \
    	--owarp=${AtlasSpaceFolder}/xfms/${OutputfMRI2StandardTransform} \
    	--oiwarp=${AtlasSpaceFolder}/xfms/${Standard2OutputfMRITransform} \
    	--motionmatdir=${fMRIFolder}/${MotionMatrixFolder} \
    	--motionmatprefix=${MotionMatrixPrefix} \
    	--ofmri=${fMRIFolder}/${NameOffMRI}_nonlin \
    	--freesurferbrainmask=${AtlasSpaceFolder}/${FreeSurferBrainMask} \
    	--biasfield=${AtlasSpaceFolder}/${BiasFieldMNI} \
    	--gdfield=${fMRIFolder}/${NameOffMRI}_gdc_warp \
    	--scoutin=${fMRIFolder}/${OrigScoutName} \
    	--scoutgdcin=${fMRIFolder}/${ScoutName}_gdc \
    	--oscout=${fMRIFolder}/${NameOffMRI}_SBRef_nonlin \
    	--jacobianin=${fMRIFolder}/${JacobianOut} \
    	--ojacobian=${fMRIFolder}/${JacobianOut}_MNI.${FinalfMRIResolution}
else
	echo "\nOne Step Resampling appears complete. Skipping to Intensity Normalization.\n"
fi

if [ ! -e "${fMRIFolder}/${NameOffMRI}_nonlin_norm.nii.gz" ] ; then
	#Intensity Normalization and Bias Removal
	log_Msg "Intensity Normalization and Bias Removal"
	${RUN} ${PipelineScripts}/IntensityNormalization.sh \
		--infmri=${fMRIFolder}/${NameOffMRI}_nonlin \
    	--biasfield=${fMRIFolder}/${BiasFieldMNI}.${FinalfMRIResolution} \
    	--jacobian=${fMRIFolder}/${JacobianOut}_MNI.${FinalfMRIResolution} \
    	--brainmask=${fMRIFolder}/${FreeSurferBrainMask}.${FinalfMRIResolution} \
    	--ofmri=${fMRIFolder}/${NameOffMRI}_nonlin_norm \
    	--inscout=${fMRIFolder}/${NameOffMRI}_SBRef_nonlin \
    	--oscout=${fMRIFolder}/${NameOffMRI}_SBRef_nonlin_norm \
    	--usejacobian=false
	
	log_Msg "mkdir -p ${ResultsFolder}"
	mkdir -p ${ResultsFolder}
	# MJ QUERY: WHY THE -r OPTIONS BELOW?
	# TBr Response: Since the copy operations are specifying individual files
	# to be copied and not directories, the recursive copy options (-r) to the
	# cp calls below definitely seem unnecessary. They should be removed in 
	# a code clean up phase when tests are in place to verify that removing them
	# has no unexpected bad side-effect.
	${RUN} cp -r ${fMRIFolder}/${NameOffMRI}_nonlin_norm.nii.gz ${ResultsFolder}/${NameOffMRI}.nii.gz
	${RUN} cp -r ${fMRIFolder}/${MovementRegressor}.txt ${ResultsFolder}/${MovementRegressor}.txt
	${RUN} cp -r ${fMRIFolder}/${MovementRegressor}_dt.txt ${ResultsFolder}/${MovementRegressor}_dt.txt
	${RUN} cp -r ${fMRIFolder}/${NameOffMRI}_SBRef_nonlin_norm.nii.gz ${ResultsFolder}/${NameOffMRI}_SBRef.nii.gz
	${RUN} cp -r ${fMRIFolder}/${JacobianOut}_MNI.${FinalfMRIResolution}.nii.gz ${ResultsFolder}/${NameOffMRI}_${JacobianOut}.nii.gz
	${RUN} cp -r ${fMRIFolder}/${FreeSurferBrainMask}.${FinalfMRIResolution}.nii.gz ${ResultsFolder}
	###Add stuff for RMS###
	${RUN} cp -r ${fMRIFolder}/Movement_RelativeRMS.txt ${ResultsFolder}/Movement_RelativeRMS.txt
	${RUN} cp -r ${fMRIFolder}/Movement_AbsoluteRMS.txt ${ResultsFolder}/Movement_AbsoluteRMS.txt
	${RUN} cp -r ${fMRIFolder}/Movement_RelativeRMS_mean.txt ${ResultsFolder}/Movement_RelativeRMS_mean.txt
	${RUN} cp -r ${fMRIFolder}/Movement_AbsoluteRMS_mean.txt ${ResultsFolder}/Movement_AbsoluteRMS_mean.txt
	###Add stuff for RMS###
else
	echo "\nIntensity Normalization appears complete. Skipping to One Step Resampling (Native).\n"
fi

if [ ! -e "${fMRIFolder}/${NameOffMRI}_gcd_mc_dc_native.nii.gz" ] ; then
	#One Step Resampling for Native Space BOLD Output
	log_Msg "One Step Resampling Native"
	log_Msg "mkdir -p ${fMRIFolder}/OneStepResampling_native"
	mkdir -p ${fMRIFolder}/OneStepResampling_native
	${RUN} ${PipelineScripts}/OneStepResampling_Native.sh \
		--workingdir=${fMRIFolder}/OneStepResampling_native \
    	--infmri=${fMRIFolder}/${OrigTCSName}.nii.gz \
    	--t1=${T1wFolder}/${T1wNativeName} \
    	--fmriresout=${FinalfMRIResolution} \
    	--fmrifolder=${fMRIFolder} \
    	--fmri2structin=${T1wFolder}/xfms/${fMRI2strOutputTransform} \
    	--owarp=${T1wFolder}/xfms/${OutputfMRI2NativeTransform} \
    	--oiwarp=${T1wFolder}/xfms/${Native2OutputfMRITransform} \
    	--motionmatdir=${fMRIFolder}/${MotionMatrixFolder} \
    	--motionmatprefix=${MotionMatrixPrefix} \
    	--ofmri=${fMRIFolder}/${NameOffMRI}_gcd_mc_dc_native \
    	--freesurferbrainmask=${T1wFolder}/${FreeSurferBrainMask} \
    	--biasfield=${T1wFolder}/${BiasField} \
    	--gdfield=${fMRIFolder}/${NameOffMRI}_gdc_warp \
    	--scoutin=${fMRIFolder}/${OrigScoutName} \
    	--scoutgdcin=${fMRIFolder}/${ScoutName}_gdc \
    	--oscout=${fMRIFolder}/${NameOffMRI}_SBRef_gdc_mc_dc_native \
    	--jacobianin=${fMRIFolder}/${JacobianOut} \
    	--ojacobian=${fMRIFolder}/${JacobianOut}_native.${FinalfMRIResolution}
else
	echo "\nOne Step Resampling (Native) appears complete. Skipping to Intensity Normalization (Native).\n"
fi

if [ ! -e "${fMRIFolder}/${NameOffMRI}_gcd_mc_dc_native_norm.nii.gz" ] ; then
	#Intensity Normalization and Bias Removal (Native)
	log_Msg "Intensity Normalization and Bias Removal (Native)"
	${RUN} ${PipelineScripts}/IntensityNormalization.sh \
		--infmri=${fMRIFolder}/${NameOffMRI}_gcd_mc_dc_native \
    	--biasfield=${fMRIFolder}/${BiasField}.${FinalfMRIResolution} \
    	--jacobian=${fMRIFolder}/${JacobianOut}_native.${FinalfMRIResolution} \
    	--brainmask=${fMRIFolder}/${FreeSurferBrainMask}_native.${FinalfMRIResolution} \
    	--ofmri=${fMRIFolder}/${NameOffMRI}_gcd_mc_dc_native_norm \
    	--inscout=${fMRIFolder}/${NameOffMRI}_SBRef_gdc_mc_dc_native \
    	--oscout=${fMRIFolder}/${NameOffMRI}_SBRef_gdc_mc_dc_native_norm \
    	--usejacobian=false
	
	log_Msg "mkdir -p ${NativeResultsFolder}"
	mkdir -p ${NativeResultsFolder}
	${RUN} cp ${fMRIFolder}/${NameOffMRI}_gcd_mc_dc_native_norm.nii.gz ${NativeResultsFolder}/${NameOffMRI}.nii.gz
	${RUN} cp ${fMRIFolder}/${MovementRegressor}.txt ${NativeResultsFolder}/${MovementRegressor}.txt
	${RUN} cp ${fMRIFolder}/${MovementRegressor}_dt.txt ${NativeResultsFolder}/${MovementRegressor}_dt.txt
	${RUN} cp ${fMRIFolder}/${NameOffMRI}_SBRef_gdc_mc_dc_native_norm.nii.gz ${NativeResultsFolder}/${NameOffMRI}_SBRef.nii.gz
	${RUN} cp ${fMRIFolder}/${JacobianOut}_native.${FinalfMRIResolution}.nii.gz ${NativeResultsFolder}/${NameOffMRI}_${JacobianOut}.nii.gz
	${RUN} cp ${fMRIFolder}/${FreeSurferBrainMask}_native.${FinalfMRIResolution}.nii.gz ${NativeResultsFolder}/${FreeSurferBrainMask}.${FinalfMRIResolution}.nii.gz
	${RUN} cp ${fMRIFolder}/Movement_RelativeRMS_native.txt ${NativeResultsFolder}/Movement_RelativeRMS.txt
	${RUN} cp ${fMRIFolder}/Movement_AbsoluteRMS_native.txt ${NativeResultsFolder}/Movement_AbsoluteRMS.txt
	${RUN} cp ${fMRIFolder}/Movement_RelativeRMS_mean_native.txt ${NativeResultsFolder}/Movement_RelativeRMS_mean.txt
	${RUN} cp ${fMRIFolder}/Movement_AbsoluteRMS_mean_native.txt ${NativeResultsFolder}/Movement_AbsoluteRMS_mean.txt
else
	echo "Um, looks like you're pretty much finished with this subject already, buddy?"
fi

log_Msg "Completed"
