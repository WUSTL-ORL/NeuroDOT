#!/bin/bash
set -x

StudyFolder="$1" #FreeSurfer "subjects" folder
Subject="$2" #Subject name
T1wFolder="$3" #FreeSurfer subject "mri" folder
AtlasSpaceFolder="$4" #standard_mesh_atlases folder
NativeFolder="$5" #copy of FreeSurfer mri folder that you can write to
FreeSurferFolder="$6" #copy of FreeSurfer subject dir that you can write to
FreeSurferInput="$7" #FreeSurfer folder for your subject <subjects/Subject>
T1wImage="$8" #path to T1w image
SurfaceAtlasDIR="${9}" #standard_mesh_atlases folder
HighResMesh="${10}" #164 (resolution of the high res mesh)
LowResMeshes="${11}" #32 (resolution of the low res  mesh)
CARET7DIR="${12}" #path to workbench command (bin_rh_linux64)

MatrixX=`mri_info --cras $FreeSurferInput/mri/orig.mgz | cut -f1 -d' '`
MatrixY=`mri_info --cras $FreeSurferInput/mri/orig.mgz | cut -f2 -d' '`
MatrixZ=`mri_info --cras $FreeSurferInput/mri/orig.mgz | cut -f3 -d' '`
Matrix1= echo "1 0 0 $MatrixX" > "$FreeSurferInput"/c_ras.mat
Matrix2= echo "0 1 0 $MatrixY" >> "$FreeSurferInput"/c_ras.mat
Matrix3= echo "0 0 1 $MatrixZ" >> "$FreeSurferInput"/c_ras.mat
Matrix4= echo "0 0 0 1">> "$FreeSurferInput"/c_ras.mat
MatrixCRAS=`echo "$Matrix1"" ""$Matrix2"" ""$Matrix3"" ""$Matrix4"`>> "$FreeSurferInput"/c_ras.mat
if [ ! -e "$AtlasSpaceFolder"/fsaverage_LR"$LowResMeshes"k ] ; then
            mkdir "$AtlasSpaceFolder"/fsaverage_LR"$LowResMeshes"k 
        fi
if [ ! -e "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k ] ; then
            mkdir "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k 
        fi
if [ ! -e "$AtlasSpaceFolder"/fsaverage ] ; then
            mkdir "$AtlasSpaceFolder"/fsaverage 
        fi
#Loop through left and right hemispheres
for Hemisphere in L R ; do
    #Set a bunch of different ways of saying left and right
    if [ $Hemisphere = "L" ] ; then
        hemisphere="l"
        Structure="CORTEX_LEFT"
    elif [ $Hemisphere = "R" ] ; then
        hemisphere="r"
        Structure="CORTEX_RIGHT"
    fi

    #native Mesh Processing
    #Convert and volumetrically register white and pial surfaces makign linear and nonlinear copies, add each to the appropriate spec file
    Types="ANATOMICAL@GRAY_WHITE ANATOMICAL@PIAL"
    i=1
    for Surface in white pial ; do
        Type=$(echo "$Types" | cut -d " " -f $i)
        Secondary=$(echo "$Type" | cut -d "@" -f 2)
        Type=$(echo "$Type" | cut -d "@" -f 1)
        if [ ! $Secondary = $Type ] ; then
            Secondary=$(echo " -surface-secondary-type ""$Secondary")
        else
            Secondary=""
        fi
        mris_convert "$FreeSurferFolder"/surf/"$hemisphere"h."$Surface" "$T1wFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
        ${CARET7DIR}/wb_command -set-structure "$T1wFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${Structure} -surface-type $Type$Secondary
        #${CARET7DIR}/wb_command -surface-apply-affine "$T1wFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii "$FreeSurferFolder"/mri/c_ras.mat "$T1wFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
        ${CARET7DIR}/wb_command -add-to-spec-file "$T1wFolder"/"$Subject".native.wb.spec $Structure "$T1wFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
        #${CARET7DIR}/wb_command -surface-apply-warpfield "$T1wFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii "$InverseAtlasTransform".nii.gz "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii -fnirt "$AtlasTransform".nii.gz
        ${CARET7DIR}/wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$Subject".native.wb.spec $Structure "$T1wFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
        i=$(( i+1 ))
    done

    #Create midthickness by averaging white and pial surfaces and use it to make inflated surfacess
    for Folder in "$T1wFolder" "$AtlasSpaceFolder" ; do
        ${CARET7DIR}/wb_command -surface-average "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii -surf "$T1wFolder"/"$Subject"."$Hemisphere".white.native.surf.gii -surf "$T1wFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii
        ${CARET7DIR}/wb_command -set-structure "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii ${Structure} -surface-type ANATOMICAL -surface-secondary-type MIDTHICKNESS
        ${CARET7DIR}/wb_command -add-to-spec-file "$T1wFolder"/"$Subject".native.wb.spec $Structure "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii

        #get number of vertices from native file
        NativeVerts=$(${CARET7DIR}/wb_command -file-information "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii | grep 'Number of Vertices:' | cut -f2 -d: | tr -d '[:space:]')

        #HCP fsaverage_LR32k used -iterations-scale 0.75. Compute new param value for native mesh density
        NativeInflationScale=$(echo "scale=4; $InflateExtraScale * 0.75 * $NativeVerts / 32492" | bc -l)

        ${CARET7DIR}/wb_command -surface-generate-inflated "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMeshes"k/"$Subject"."$Hemisphere".inflated.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMeshes"k/"$Subject"."$Hemisphere".very_inflated.native.surf.gii -iterations-scale 2.5
        ${CARET7DIR}/wb_command -add-to-spec-file "$T1wFolder"/"$Subject".native.wb.spec $Structure "$AtlasSpaceFolder"/fsaverage_LR"$LowResMeshes"k/"$Subject"."$Hemisphere".inflated.native.surf.gii
        ${CARET7DIR}/wb_command -add-to-spec-file "$T1wFolder"/"$Subject".native.wb.spec $Structure "$AtlasSpaceFolder"/fsaverage_LR"$LowResMeshes"k/"$Subject"."$Hemisphere".very_inflated.native.surf.gii
    done

    #Convert original and registered spherical surfaces and add them to the nonlinear spec file
    for Surface in sphere.reg sphere ; do
        mris_convert "$FreeSurferFolder"/surf/"$hemisphere"h."$Surface" "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii
        ${CARET7DIR}/wb_command -set-structure "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${Structure} -surface-type SPHERICAL
    done
    ${CARET7DIR}/wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$Subject".native.wb.spec $Structure "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.native.surf.gii

    #Add more files to the spec file and convert other FreeSurfer surface data to metric/GIFTI including sulc, curv, and thickness.
    for Map in sulc@sulc@Sulc thickness@thickness@Thickness curv@curvature@Curvature ; do
        fsname=$(echo $Map | cut -d "@" -f 1)
        wbname=$(echo $Map | cut -d "@" -f 2)
        mapname=$(echo $Map | cut -d "@" -f 3)
        mris_convert -c "$FreeSurferFolder"/surf/"$hemisphere"h."$fsname" "$FreeSurferFolder"/surf/"$hemisphere"h.white "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii
        ${CARET7DIR}/wb_command -set-structure "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii ${Structure}
        ${CARET7DIR}/wb_command -metric-math "var * -1" "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -var var "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii
        ${CARET7DIR}/wb_command -set-map-names "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -map 1 "$Subject"_"$Hemisphere"_"$mapname"
        ${CARET7DIR}/wb_command -metric-palette "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true
    done
    #Thickness specific operations
    ${CARET7DIR}/wb_command -metric-math "abs(thickness)" "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii -var thickness "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii
    ${CARET7DIR}/wb_command -metric-palette "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false
    ${CARET7DIR}/wb_command -metric-math "thickness > 0" "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii -var thickness "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii
    ${CARET7DIR}/wb_command -metric-fill-holes "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
    ${CARET7DIR}/wb_command -metric-remove-islands "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
    ${CARET7DIR}/wb_command -set-map-names "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii -map 1 "$Subject"_"$Hemisphere"_ROI
    ${CARET7DIR}/wb_command -metric-dilate "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii 10 "$T1wFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii -nearest
    ${CARET7DIR}/wb_command -metric-dilate "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".curvature.native.shape.gii "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii 10 "$T1wFolder"/"$Subject"."$Hemisphere".curvature.native.shape.gii -nearest

    #Label operations
    for Map in aparc aparc.a2009s ; do #Remove BA because it doesn't convert properly
        if [ -e "$FreeSurferFolder"/label/"$hemisphere"h."$Map".annot ] ; then
            mris_convert --annot "$FreeSurferFolder"/label/"$hemisphere"h."$Map".annot "$FreeSurferFolder"/surf/"$hemisphere"h.white "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii
            ${CARET7DIR}/wb_command -set-structure "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii $Structure
            ${CARET7DIR}/wb_command -set-map-names "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii -map 1 "$Subject"_"$Hemisphere"_"$Map"
            ${CARET7DIR}/wb_command -gifti-label-add-prefix "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii "${Hemisphere}_" "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii
        fi
    done
    #End main native mesh processing

    #Copy Atlas Files
    cp "$SurfaceAtlasDIR"/fs_"$Hemisphere"/fsaverage."$Hemisphere".sphere."$HighResMesh"k_fs_"$Hemisphere".surf.gii "$AtlasSpaceFolder"/fsaverage/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_"$Hemisphere".surf.gii
    cp "$SurfaceAtlasDIR"/fs_"$Hemisphere"/fs_"$Hemisphere"-to-fs_LR_fsaverage."$Hemisphere"_LR.spherical_std."$HighResMesh"k_fs_"$Hemisphere".surf.gii "$AtlasSpaceFolder"/fsaverage/"$Subject"."$Hemisphere".def_sphere."$HighResMesh"k_fs_"$Hemisphere".surf.gii
    cp "$SurfaceAtlasDIR"/fsaverage."$Hemisphere"_LR.spherical_std."$HighResMesh"k_fs_LR.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii
    ${CARET7DIR}/wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$Subject"."$HighResMesh"k_fs_LR.wb.spec $Structure "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii
    cp "$SurfaceAtlasDIR"/"$Hemisphere".atlasroi."$HighResMesh"k_fs_LR.shape.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".atlasroi."$HighResMesh"k_fs_LR.shape.gii
    cp "$SurfaceAtlasDIR"/"$Hemisphere".refsulc."$HighResMesh"k_fs_LR.shape.gii "$AtlasSpaceFolder"/${Subject}.${Hemisphere}.refsulc."$HighResMesh"k_fs_LR.shape.gii
    if [ -e "$SurfaceAtlasDIR"/colin.cerebral."$Hemisphere".flat."$HighResMesh"k_fs_LR.surf.gii ] ; then
        cp "$SurfaceAtlasDIR"/colin.cerebral."$Hemisphere".flat."$HighResMesh"k_fs_LR.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".flat."$HighResMesh"k_fs_LR.surf.gii
        ${CARET7DIR}/wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$Subject"."$HighResMesh"k_fs_LR.wb.spec $Structure "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".flat."$HighResMesh"k_fs_LR.surf.gii
    fi

    #Concatenate FS registration to FS --> FS_LR registration
    ${CARET7DIR}/wb_command -surface-sphere-project-unproject "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.reg.native.surf.gii "$AtlasSpaceFolder"/fsaverage/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_"$Hemisphere".surf.gii "$AtlasSpaceFolder"/fsaverage/"$Subject"."$Hemisphere".def_sphere."$HighResMesh"k_fs_"$Hemisphere".surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.reg.reg_LR.native.surf.gii

    #Make FreeSurfer Registration Areal Distortion Maps
    ${CARET7DIR}/wb_command -surface-vertex-areas "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.native.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.native.shape.gii
    ${CARET7DIR}/wb_command -surface-vertex-areas "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.reg.reg_LR.native.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.reg.reg_LR.native.shape.gii
    ${CARET7DIR}/wb_command -metric-math "ln(spherereg / sphere) / ln(2)" "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".ArealDistortion_FS.native.shape.gii -var sphere "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.native.shape.gii -var spherereg "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.reg.reg_LR.native.shape.gii
    rm "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.native.shape.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.reg.reg_LR.native.shape.gii
    ${CARET7DIR}/wb_command -set-map-names "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".ArealDistortion_FS.native.shape.gii -map 1 "$Subject"_"$Hemisphere"_Areal_Distortion_FS
    ${CARET7DIR}/wb_command -metric-palette "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".ArealDistortion_FS.native.shape.gii MODE_AUTO_SCALE -palette-name ROY-BIG-BL -thresholding THRESHOLD_TYPE_NORMAL THRESHOLD_TEST_SHOW_OUTSIDE -1 1

    ${CARET7DIR}/wb_command -surface-distortion "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.native.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.reg.reg_LR.native.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".EdgeDistortion_FS.native.shape.gii -edge-method

    ${CARET7DIR}/wb_command -surface-distortion "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.native.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.reg.reg_LR.native.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".Strain_FS.native.shape.gii -local-affine-method
    ${CARET7DIR}/wb_command -metric-merge "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainJ_FS.native.shape.gii -metric "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".Strain_FS.native.shape.gii -column 1
    ${CARET7DIR}/wb_command -metric-merge "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainR_FS.native.shape.gii -metric "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".Strain_FS.native.shape.gii -column 2
    ${CARET7DIR}/wb_command -metric-math "ln(var) / ln (2)" "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainJ_FS.native.shape.gii -var var "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainJ_FS.native.shape.gii
    ${CARET7DIR}/wb_command -metric-math "ln(var) / ln (2)" "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainR_FS.native.shape.gii -var var "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainR_FS.native.shape.gii
    rm "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".Strain_FS.native.shape.gii

    #If desired, run MSMSulc folding-based registration to FS_LR initialized with FS affine
    if [ ${RegName} = "MSMSulc" ] ; then
        #Calculate Affine Transform and Apply
        if [ ! -e "$AtlasSpaceFolder"/MSMSulc ] ; then
            mkdir "$AtlasSpaceFolder"/MSMSulc
        fi
        ${CARET7DIR}/wb_command -surface-affine-regression "$AtlasSpaceFolder"/${Subject}.${Hemisphere}.sphere.native.surf.gii "$AtlasSpaceFolder"/${Subject}.${Hemisphere}.sphere.reg.reg_LR.native.surf.gii "$AtlasSpaceFolder"/MSMSulc/${Hemisphere}.mat
        ${CARET7DIR}/wb_command -surface-apply-affine "$AtlasSpaceFolder"/${Subject}.${Hemisphere}.sphere.native.surf.gii "$AtlasSpaceFolder"/MSMSulc/${Hemisphere}.mat "$AtlasSpaceFolder"/MSMSulc/${Hemisphere}.sphere_rot.surf.gii
        ${CARET7DIR}/wb_command -surface-modify-sphere "$AtlasSpaceFolder"/MSMSulc/${Hemisphere}.sphere_rot.surf.gii 100 "$AtlasSpaceFolder"/MSMSulc/${Hemisphere}.sphere_rot.surf.gii
        cp "$AtlasSpaceFolder"/MSMSulc/${Hemisphere}.sphere_rot.surf.gii "$AtlasSpaceFolder"/${Subject}.${Hemisphere}.sphere.rot.native.surf.gii
        DIR=$(pwd)
        cd "$AtlasSpaceFolder"/MSMSulc
        #Register using FreeSurfer Sulc Folding Map Using MSM Algorithm Configured for Reduced Distortion
        #${MSMBINDIR}/msm --version
        #${MSMBINDIR}/msm --levels=4 --conf=${MSMCONFIGDIR}/allparameterssulcDRconf --inmesh="$AtlasSpaceFolder"/${Subject}.${Hemisphere}.sphere.rot.native.surf.gii --trans="$AtlasSpaceFolder"/${Subject}.${Hemisphere}.sphere.rot.native.surf.gii --refmesh="$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii --indata="$AtlasSpaceFolder"/${Subject}.${Hemisphere}.sulc.native.shape.gii --refdata="$AtlasSpaceFolder"/${Subject}.${Hemisphere}.refsulc."$HighResMesh"k_fs_LR.shape.gii --out="$AtlasSpaceFolder"/MSMSulc/${Hemisphere}. --verbose
        ${MSMBINDIR}/msm --conf=${MSMCONFIGDIR}/MSMSulcStrainFinalconf --inmesh="$AtlasSpaceFolder"/${Subject}.${Hemisphere}.sphere.rot.native.surf.gii --refmesh="$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii --indata="$AtlasSpaceFolder"/${Subject}.${Hemisphere}.sulc.native.shape.gii --refdata="$AtlasSpaceFolder"/${Subject}.${Hemisphere}.refsulc."$HighResMesh"k_fs_LR.shape.gii --out="$AtlasSpaceFolder"/MSMSulc/${Hemisphere}. --verbose
        cp ${MSMCONFIGDIR}/MSMSulcStrainFinalconf "$AtlasSpaceFolder"/MSMSulc/${Hemisphere}.logdir/conf
        cd $DIR
        #cp "$AtlasSpaceFolder"/MSMSulc/${Hemisphere}.HIGHRES_transformed.surf.gii "$AtlasSpaceFolder"/${Subject}.${Hemisphere}.sphere.MSMSulc.native.surf.gii
        cp "$AtlasSpaceFolder"/MSMSulc/${Hemisphere}.sphere.reg.surf.gii "$AtlasSpaceFolder"/${Subject}.${Hemisphere}.sphere.MSMSulc.native.surf.gii
        ${CARET7DIR}/wb_command -set-structure "$AtlasSpaceFolder"/${Subject}.${Hemisphere}.sphere.MSMSulc.native.surf.gii ${Structure}

        #Make MSMSulc Registration Areal Distortion Maps
        ${CARET7DIR}/wb_command -surface-vertex-areas "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.native.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.native.shape.gii
        ${CARET7DIR}/wb_command -surface-vertex-areas "$AtlasSpaceFolder"/${Subject}.${Hemisphere}.sphere.MSMSulc.native.surf.gii "$AtlasSpaceFolder"/${Subject}.${Hemisphere}.sphere.MSMSulc.native.shape.gii
        ${CARET7DIR}/wb_command -metric-math "ln(spherereg / sphere) / ln(2)" "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".ArealDistortion_MSMSulc.native.shape.gii -var sphere "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.native.shape.gii -var spherereg "$AtlasSpaceFolder"/${Subject}.${Hemisphere}.sphere.MSMSulc.native.shape.gii
        rm "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.native.shape.gii "$AtlasSpaceFolder"/${Subject}.${Hemisphere}.sphere.MSMSulc.native.shape.gii
        ${CARET7DIR}/wb_command -set-map-names "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".ArealDistortion_MSMSulc.native.shape.gii -map 1 "$Subject"_"$Hemisphere"_Areal_Distortion_MSMSulc
        ${CARET7DIR}/wb_command -metric-palette "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".ArealDistortion_MSMSulc.native.shape.gii MODE_AUTO_SCALE -palette-name ROY-BIG-BL -thresholding THRESHOLD_TYPE_NORMAL THRESHOLD_TEST_SHOW_OUTSIDE -1 1

        ${CARET7DIR}/wb_command -surface-distortion "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.native.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.MSMSulc.native.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".EdgeDistortion_MSMSulc.native.shape.gii -edge-method

        ${CARET7DIR}/wb_command -surface-distortion "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.native.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere.MSMSulc.native.surf.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".Strain_MSMSulc.native.shape.gii -local-affine-method
        ${CARET7DIR}/wb_command -metric-merge "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainJ_MSMSulc.native.shape.gii -metric "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".Strain_MSMSulc.native.shape.gii -column 1
        ${CARET7DIR}/wb_command -metric-merge "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainR_MSMSulc.native.shape.gii -metric "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".Strain_MSMSulc.native.shape.gii -column 2
        ${CARET7DIR}/wb_command -metric-math "ln(var) / ln (2)" "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainJ_MSMSulc.native.shape.gii -var var "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainJ_MSMSulc.native.shape.gii
        ${CARET7DIR}/wb_command -metric-math "ln(var) / ln (2)" "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainR_MSMSulc.native.shape.gii -var var "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainR_MSMSulc.native.shape.gii
        rm "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".Strain_MSMSulc.native.shape.gii

        RegSphere="${AtlasSpaceFolder}/${Subject}.${Hemisphere}.sphere.MSMSulc.native.surf.gii"
    else
        RegSphere="${AtlasSpaceFolder}/${Subject}.${Hemisphere}.sphere.reg.reg_LR.native.surf.gii"
    fi

    #Ensure no zeros in atlas medial wall ROI
    ${CARET7DIR}/wb_command -metric-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".atlasroi."$HighResMesh"k_fs_LR.shape.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii ${RegSphere} BARYCENTRIC "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".atlasroi.native.shape.gii -largest
    ${CARET7DIR}/wb_command -metric-math "(atlas + individual) > 0" "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii -var atlas "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".atlasroi.native.shape.gii -var individual "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
    ${CARET7DIR}/wb_command -metric-mask "$T1wFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".thickness.native.shape.gii
    ${CARET7DIR}/wb_command -metric-mask "$T1wFolder"/"$Subject"."$Hemisphere".curvature.native.shape.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".curvature.native.shape.gii


    #Populate Highres fs_LR spec file.  Deform surfaces and other data according to native to folding-based registration selected above.  Regenerate inflated surfaces.
    for Surface in white midthickness pial ; do
        ${CARET7DIR}/wb_command -surface-resample "$T1wFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${RegSphere} "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii BARYCENTRIC "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere"."$Surface"."$HighResMesh"k_fs_LR.surf.gii
        ${CARET7DIR}/wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$Subject"."$HighResMesh"k_fs_LR.wb.spec $Structure "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere"."$Surface"."$HighResMesh"k_fs_LR.surf.gii
    done

    #HCP fsaverage_LR32k used -iterations-scale 0.75. Compute new param value for high res mesh density
    HighResInflationScale=$(echo "scale=4; $InflateExtraScale * 0.75 * $HighResMesh / 32" | bc -l)

    ${CARET7DIR}/wb_command -surface-generate-inflated "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".inflated."$HighResMesh"k_fs_LR.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".very_inflated."$HighResMesh"k_fs_LR.surf.gii -iterations-scale 2.5
    ${CARET7DIR}/wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$Subject"."$HighResMesh"k_fs_LR.wb.spec $Structure "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".inflated."$HighResMesh"k_fs_LR.surf.gii
    ${CARET7DIR}/wb_command -add-to-spec-file "$AtlasSpaceFolder"/"$Subject"."$HighResMesh"k_fs_LR.wb.spec $Structure "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".very_inflated."$HighResMesh"k_fs_LR.surf.gii

    for Map in thickness curvature ; do
        ${CARET7DIR}/wb_command -metric-resample "$T1wFolder"/"$Subject"."$Hemisphere"."$Map".native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$Map"."$HighResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii -current-roi "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
        ${CARET7DIR}/wb_command -metric-mask "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$Map"."$HighResMesh"k_fs_LR.shape.gii "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".atlasroi."$HighResMesh"k_fs_LR.shape.gii "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere"."$Map"."$HighResMesh"k_fs_LR.shape.gii
    done
    ${CARET7DIR}/wb_command -metric-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".ArealDistortion_FS.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".ArealDistortion_FS."$HighResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii
    ${CARET7DIR}/wb_command -metric-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".EdgeDistortion_FS.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".EdgeDistortion_FS."$HighResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii
    ${CARET7DIR}/wb_command -metric-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainJ_FS.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainJ_FS."$HighResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii
    ${CARET7DIR}/wb_command -metric-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainR_FS.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainR_FS."$HighResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii
    if [ ${RegName} = "MSMSulc" ] ; then
        ${CARET7DIR}/wb_command -metric-resample "$T1wFolder"/"$Subject"."$Hemisphere".ArealDistortion_MSMSulc.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".ArealDistortion_MSMSulc."$HighResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii
        ${CARET7DIR}/wb_command -metric-resample "$T1wFolder"/"$Subject"."$Hemisphere".EdgeDistortion_MSMSulc.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".EdgeDistortion_MSMSulc."$HighResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii
        ${CARET7DIR}/wb_command -metric-resample "$T1wFolder"/"$Subject"."$Hemisphere".StrainJ_MSMSulc.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainJ_MSMSulc."$HighResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii
        ${CARET7DIR}/wb_command -metric-resample "$T1wFolder"/"$Subject"."$Hemisphere".StrainR_MSMSulc.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainR_MSMSulc."$HighResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii
    fi
    ${CARET7DIR}/wb_command -metric-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sulc.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sulc."$HighResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$HighResMesh"k/"$Subject"."$Hemisphere".midthickness."$HighResMesh"k_fs_LR.surf.gii

    for Map in aparc aparc.a2009s ; do #Remove BA because it doesn't convert properly
        if [ -e "$FreeSurferFolder"/label/"$hemisphere"h."$Map".annot ] ; then
            ${CARET7DIR}/wb_command -label-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii ${RegSphere} "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sphere."$HighResMesh"k_fs_LR.surf.gii BARYCENTRIC "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$Map"."$HighResMesh"k_fs_LR.label.gii -largest
        fi
    done

    for LowResMesh in ${LowResMeshes} ; do
        #Copy Atlas Files
        cp "$SurfaceAtlasDIR"/"$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii
        ${CARET7DIR}/wb_command -add-to-spec-file "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$LowResMesh"k_fs_LR.wb.spec $Structure "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii
        cp "$AtlasSpaceFolder"/"$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii
        if [ -e "$SurfaceAtlasDIR"/colin.cerebral."$Hemisphere".flat."$LowResMesh"k_fs_LR.surf.gii ] ; then
            cp "$SurfaceAtlasDIR"/colin.cerebral."$Hemisphere".flat."$LowResMesh"k_fs_LR.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".flat."$LowResMesh"k_fs_LR.surf.gii
            ${CARET7DIR}/wb_command -add-to-spec-file "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$LowResMesh"k_fs_LR.wb.spec $Structure "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".flat."$LowResMesh"k_fs_LR.surf.gii
        fi

        #Create downsampled fs_LR spec files.
        for Surface in white midthickness pial ; do
            ${CARET7DIR}/wb_command -surface-resample "$T1wFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${RegSphere} "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii BARYCENTRIC "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere"."$Surface"."$LowResMesh"k_fs_LR.surf.gii
            ${CARET7DIR}/wb_command -add-to-spec-file "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$LowResMesh"k_fs_LR.wb.spec $Structure "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere"."$Surface"."$LowResMesh"k_fs_LR.surf.gii
        done

        #HCP fsaverage_LR32k used -iterations-scale 0.75. Recalculate in case using a different mesh
        LowResInflationScale=$(echo "scale=4; $InflateExtraScale * 0.75 * $LowResMesh / 32" | bc -l)

        ${CARET7DIR}/wb_command -surface-generate-inflated "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".inflated."$LowResMesh"k_fs_LR.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".very_inflated."$LowResMesh"k_fs_LR.surf.gii -iterations-scale 2.5
        ${CARET7DIR}/wb_command -add-to-spec-file "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$LowResMesh"k_fs_LR.wb.spec $Structure "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".inflated."$LowResMesh"k_fs_LR.surf.gii
        ${CARET7DIR}/wb_command -add-to-spec-file "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$LowResMesh"k_fs_LR.wb.spec $Structure "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".very_inflated."$LowResMesh"k_fs_LR.surf.gii

        for Map in sulc thickness curvature ; do
            ${CARET7DIR}/wb_command -metric-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$Map".native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere"."$Map"."$LowResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii -current-roi "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii
            ${CARET7DIR}/wb_command -metric-mask "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere"."$Map"."$LowResMesh"k_fs_LR.shape.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere"."$Map"."$LowResMesh"k_fs_LR.shape.gii
        done
        ${CARET7DIR}/wb_command -metric-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".ArealDistortion_FS.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".ArealDistortion_FS."$LowResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii
        ${CARET7DIR}/wb_command -metric-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".EdgeDistortion_FS.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".EdgeDistortion_FS."$LowResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii
        ${CARET7DIR}/wb_command -metric-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainJ_FS.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".StrainJ_FS."$LowResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii
        ${CARET7DIR}/wb_command -metric-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainR_FS.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".StrainR_FS."$LowResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii
        if [ ${RegName} = "MSMSulc" ] ; then
            ${CARET7DIR}/wb_command -metric-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".ArealDistortion_MSMSulc.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".ArealDistortion_MSMSulc."$LowResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii
            ${CARET7DIR}/wb_command -metric-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".EdgeDistortion_MSMSulc.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".EdgeDistortion_MSMSulc."$LowResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii
            ${CARET7DIR}/wb_command -metric-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainJ_MSMSulc.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".StrainJ_MSMSulc."$LowResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii
            ${CARET7DIR}/wb_command -metric-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".StrainR_MSMSulc.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".StrainR_MSMSulc."$LowResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii
        fi
        ${CARET7DIR}/wb_command -metric-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere".sulc.native.shape.gii ${RegSphere} "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii ADAP_BARY_AREA "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sulc."$LowResMesh"k_fs_LR.shape.gii -area-surfs "$T1wFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii

        for Map in aparc aparc.a2009s ; do #Remove BA because it doesn't convert properly
            if [ -e "$FreeSurferFolder"/label/"$hemisphere"h."$Map".annot ] ; then
                ${CARET7DIR}/wb_command -label-resample "$AtlasSpaceFolder"/"$Subject"."$Hemisphere"."$Map".native.label.gii ${RegSphere} "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii BARYCENTRIC "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere"."$Map"."$LowResMesh"k_fs_LR.label.gii -largest
            fi
        done

        #Create downsampled fs_LR spec file in structural space.
        for Surface in white midthickness pial ; do
            ${CARET7DIR}/wb_command -surface-resample "$T1wFolder"/"$Subject"."$Hemisphere"."$Surface".native.surf.gii ${RegSphere} "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii BARYCENTRIC "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere"."$Surface"."$LowResMesh"k_fs_LR.surf.gii
            ${CARET7DIR}/wb_command -add-to-spec-file "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$LowResMesh"k_fs_LR.wb.spec $Structure "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere"."$Surface"."$LowResMesh"k_fs_LR.surf.gii
        done

        #HCP fsaverage_LR32k used -iterations-scale 0.75. Recalculate in case using a different mesh
        LowResInflationScale=$(echo "scale=4; $InflateExtraScale * 0.75 * $LowResMesh / 32" | bc -l)

        ${CARET7DIR}/wb_command -surface-generate-inflated "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".inflated."$LowResMesh"k_fs_LR.surf.gii "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".very_inflated."$LowResMesh"k_fs_LR.surf.gii -iterations-scale 2.5
        ${CARET7DIR}/wb_command -add-to-spec-file "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$LowResMesh"k_fs_LR.wb.spec $Structure "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".inflated."$LowResMesh"k_fs_LR.surf.gii
        ${CARET7DIR}/wb_command -add-to-spec-file "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$LowResMesh"k_fs_LR.wb.spec $Structure "$AtlasSpaceFolder"/fsaverage_LR"$LowResMesh"k/"$Subject"."$Hemisphere".very_inflated."$LowResMesh"k_fs_LR.surf.gii
    done
done

STRINGII=""
for LowResMesh in ${LowResMeshes} ; do
    STRINGII=$(echo "${STRINGII}${AtlasSpaceFolder}/fsaverage_LR${LowResMesh}k@${LowResMesh}k_fs_LR@atlasroi ")
done

#Create CIFTI Files
for STRING in "$AtlasSpaceFolder"/@native@roi "$AtlasSpaceFolder"@"$HighResMesh"k_fs_LR@atlasroi ${STRINGII} ; do
    Folder=$(echo $STRING | cut -d "@" -f 1)
    Mesh=$(echo $STRING | cut -d "@" -f 2)
    ROI=$(echo $STRING | cut -d "@" -f 3)

    ${CARET7DIR}/wb_command -cifti-create-dense-scalar "$AtlasSpaceFolder"/"$Subject".sulc."$Mesh".dscalar.nii -left-metric "$T1wFolder"/"$Subject".L.sulc."$Mesh".shape.gii -right-metric "$AtlasSpaceFolder"/"$Subject".R.sulc."$Mesh".shape.gii
    ${CARET7DIR}/wb_command -set-map-names "$AtlasSpaceFolder"/"$Subject".sulc."$Mesh".dscalar.nii -map 1 "${Subject}_Sulc"
    ${CARET7DIR}/wb_command -cifti-palette "$AtlasSpaceFolder"/"$Subject".sulc."$Mesh".dscalar.nii MODE_AUTO_SCALE_PERCENTAGE "$AtlasSpaceFolder"/"$Subject".sulc."$Mesh".dscalar.nii -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true

    ${CARET7DIR}/wb_command -cifti-create-dense-scalar "$AtlasSpaceFolder"/"$Subject".curvature."$Mesh".dscalar.nii -left-metric "$AtlasSpaceFolder"/"$Subject".L.curvature."$Mesh".shape.gii -roi-left "$AtlasSpaceFolder"/"$Subject".L."$ROI"."$Mesh".shape.gii -right-metric "$AtlasSpaceFolder"/"$Subject".R.curvature."$Mesh".shape.gii -roi-right "$AtlasSpaceFolder"/"$Subject".R."$ROI"."$Mesh".shape.gii
    ${CARET7DIR}/wb_command -set-map-names "$AtlasSpaceFolder"/"$Subject".curvature."$Mesh".dscalar.nii -map 1 "${Subject}_Curvature"
    ${CARET7DIR}/wb_command -cifti-palette "$AtlasSpaceFolder"/"$Subject".curvature."$Mesh".dscalar.nii MODE_AUTO_SCALE_PERCENTAGE "$AtlasSpaceFolder"/"$Subject".curvature."$Mesh".dscalar.nii -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true

    ${CARET7DIR}/wb_command -cifti-create-dense-scalar "$AtlasSpaceFolder"/"$Subject".thickness."$Mesh".dscalar.nii -left-metric "$AtlasSpaceFolder"/"$Subject".L.thickness."$Mesh".shape.gii -roi-left "$AtlasSpaceFolder"/"$Subject".L."$ROI"."$Mesh".shape.gii -right-metric "$AtlasSpaceFolder"/"$Subject".R.thickness."$Mesh".shape.gii -roi-right "$AtlasSpaceFolder"/"$Subject".R."$ROI"."$Mesh".shape.gii
    ${CARET7DIR}/wb_command -set-map-names "$AtlasSpaceFolder"/"$Subject".thickness."$Mesh".dscalar.nii -map 1 "${Subject}_Thickness"
    ${CARET7DIR}/wb_command -cifti-palette "$AtlasSpaceFolder"/"$Subject".thickness."$Mesh".dscalar.nii MODE_AUTO_SCALE_PERCENTAGE "$AtlasSpaceFolder"/"$Subject".thickness."$Mesh".dscalar.nii -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false

    ${CARET7DIR}/wb_command -cifti-create-dense-scalar "$AtlasSpaceFolder"/"$Subject".ArealDistortion_FS."$Mesh".dscalar.nii -left-metric "$AtlasSpaceFolder"/"$Subject".L.ArealDistortion_FS."$Mesh".shape.gii -right-metric "$AtlasSpaceFolder"/"$Subject".R.ArealDistortion_FS."$Mesh".shape.gii
    ${CARET7DIR}/wb_command -set-map-names "$AtlasSpaceFolder"/"$Subject".ArealDistortion_FS."$Mesh".dscalar.nii -map 1 "${Subject}_ArealDistortion_FS"
    ${CARET7DIR}/wb_command -cifti-palette "$AtlasSpaceFolder"/"$Subject".ArealDistortion_FS."$Mesh".dscalar.nii MODE_USER_SCALE "$AtlasSpaceFolder"/"$Subject".ArealDistortion_FS."$Mesh".dscalar.nii -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false

    ${CARET7DIR}/wb_command -cifti-create-dense-scalar "$AtlasSpaceFolder"/"$Subject".EdgeDistortion_FS."$Mesh".dscalar.nii -left-metric "$AtlasSpaceFolder"/"$Subject".L.EdgeDistortion_FS."$Mesh".shape.gii -right-metric "$AtlasSpaceFolder"/"$Subject".R.EdgeDistortion_FS."$Mesh".shape.gii
    ${CARET7DIR}/wb_command -set-map-names "$AtlasSpaceFolder"/"$Subject".EdgeDistortion_FS."$Mesh".dscalar.nii -map 1 "${Subject}_EdgeDistortion_FS"
    ${CARET7DIR}/wb_command -cifti-palette "$AtlasSpaceFolder"/"$Subject".EdgeDistortion_FS."$Mesh".dscalar.nii MODE_USER_SCALE "$AtlasSpaceFolder"/"$Subject".EdgeDistortion_FS."$Mesh".dscalar.nii -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false

    ${CARET7DIR}/wb_command -cifti-create-dense-scalar "$AtlasSpaceFolder"/"$Subject".StrainJ_FS."$Mesh".dscalar.nii -left-metric "$AtlasSpaceFolder"/"$Subject".L.StrainJ_FS."$Mesh".shape.gii -right-metric "$AtlasSpaceFolder"/"$Subject".R.StrainJ_FS."$Mesh".shape.gii
    ${CARET7DIR}/wb_command -set-map-names "$AtlasSpaceFolder"/"$Subject".StrainJ_FS."$Mesh".dscalar.nii -map 1 "${Subject}_StrainJ_FS"
    ${CARET7DIR}/wb_command -cifti-palette "$AtlasSpaceFolder"/"$Subject".StrainJ_FS."$Mesh".dscalar.nii MODE_USER_SCALE "$AtlasSpaceFolder"/"$Subject".StrainJ_FS."$Mesh".dscalar.nii -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false

    ${CARET7DIR}/wb_command -cifti-create-dense-scalar "$AtlasSpaceFolder"/"$Subject".StrainR_FS."$Mesh".dscalar.nii -left-metric "$AtlasSpaceFolder"/"$Subject".L.StrainR_FS."$Mesh".shape.gii -right-metric "$AtlasSpaceFolder"/"$Subject".R.StrainR_FS."$Mesh".shape.gii
    ${CARET7DIR}/wb_command -set-map-names "$AtlasSpaceFolder"/"$Subject".StrainR_FS."$Mesh".dscalar.nii -map 1 "${Subject}_StrainR_FS"
    ${CARET7DIR}/wb_command -cifti-palette "$AtlasSpaceFolder"/"$Subject".StrainR_FS."$Mesh".dscalar.nii MODE_USER_SCALE "$AtlasSpaceFolder"/"$Subject".StrainR_FS."$Mesh".dscalar.nii -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false

    if [ ${RegName} = "MSMSulc" ] ; then
        ${CARET7DIR}/wb_command -cifti-create-dense-scalar "$AtlasSpaceFolder"/"$Subject".ArealDistortion_MSMSulc."$Mesh".dscalar.nii -left-metric "$AtlasSpaceFolder"/"$Subject".L.ArealDistortion_MSMSulc."$Mesh".shape.gii -right-metric "$AtlasSpaceFolder"/"$Subject".R.ArealDistortion_MSMSulc."$Mesh".shape.gii
        ${CARET7DIR}/wb_command -set-map-names "$AtlasSpaceFolder"/"$Subject".ArealDistortion_MSMSulc."$Mesh".dscalar.nii -map 1 "${Subject}_ArealDistortion_MSMSulc"
        ${CARET7DIR}/wb_command -cifti-palette "$AtlasSpaceFolder"/"$Subject".ArealDistortion_MSMSulc."$Mesh".dscalar.nii MODE_USER_SCALE "$AtlasSpaceFolder"/"$Subject".ArealDistortion_MSMSulc."$Mesh".dscalar.nii -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false

        ${CARET7DIR}/wb_command -cifti-create-dense-scalar "$AtlasSpaceFolder"/"$Subject".EdgeDistortion_MSMSulc."$Mesh".dscalar.nii -left-metric "$AtlasSpaceFolder"/"$Subject".L.EdgeDistortion_MSMSulc."$Mesh".shape.gii -right-metric "$AtlasSpaceFolder"/"$Subject".R.EdgeDistortion_MSMSulc."$Mesh".shape.gii
        ${CARET7DIR}/wb_command -set-map-names "$AtlasSpaceFolder"/"$Subject".EdgeDistortion_MSMSulc."$Mesh".dscalar.nii -map 1 "${Subject}_EdgeDistortion_MSMSulc"
        ${CARET7DIR}/wb_command -cifti-palette "$AtlasSpaceFolder"/"$Subject".EdgeDistortion_MSMSulc."$Mesh".dscalar.nii MODE_USER_SCALE "$AtlasSpaceFolder"/"$Subject".EdgeDistortion_MSMSulc."$Mesh".dscalar.nii -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false

        ${CARET7DIR}/wb_command -cifti-create-dense-scalar "$AtlasSpaceFolder"/"$Subject".StrainJ_MSMSulc."$Mesh".dscalar.nii -left-metric "$AtlasSpaceFolder"/"$Subject".L.StrainJ_MSMSulc."$Mesh".shape.gii -right-metric "$AtlasSpaceFolder"/"$Subject".R.StrainJ_MSMSulc."$Mesh".shape.gii
        ${CARET7DIR}/wb_command -set-map-names "$AtlasSpaceFolder"/"$Subject".StrainJ_MSMSulc."$Mesh".dscalar.nii -map 1 "${Subject}_StrainJ_MSMSulc"
        ${CARET7DIR}/wb_command -cifti-palette "$AtlasSpaceFolder"/"$Subject".StrainJ_MSMSulc."$Mesh".dscalar.nii MODE_USER_SCALE "$AtlasSpaceFolder"/"$Subject".StrainJ_MSMSulc."$Mesh".dscalar.nii -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false

        ${CARET7DIR}/wb_command -cifti-create-dense-scalar "$AtlasSpaceFolder"/"$Subject".StrainR_MSMSulc."$Mesh".dscalar.nii -left-metric "$AtlasSpaceFolder"/"$Subject".L.StrainR_MSMSulc."$Mesh".shape.gii -right-metric "$AtlasSpaceFolder"/"$Subject".R.StrainR_MSMSulc."$Mesh".shape.gii
        ${CARET7DIR}/wb_command -set-map-names "$AtlasSpaceFolder"/"$Subject".StrainR_MSMSulc."$Mesh".dscalar.nii -map 1 "${Subject}_StrainR_MSMSulc"
        ${CARET7DIR}/wb_command -cifti-palette "$AtlasSpaceFolder"/"$Subject".StrainR_MSMSulc."$Mesh".dscalar.nii MODE_USER_SCALE "$AtlasSpaceFolder"/"$Subject".StrainR_MSMSulc."$Mesh".dscalar.nii -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false
    fi

    for Map in aparc aparc.a2009s ; do #Remove BA because it doesn't convert properly
        if [ -e "$T1wFolder"/"$Subject".L.${Map}."$Mesh".label.gii ] ; then
            ${CARET7DIR}/wb_command -cifti-create-label "$T1wFolder"/"$Subject".${Map}."$Mesh".dlabel.nii -left-label "$T1wFolder"/"$Subject".L.${Map}."$Mesh".label.gii -roi-left "$T1wFolder"/"$Subject".L."$ROI"."$Mesh".shape.gii -right-label "$T1wFolder"/"$Subject".R.${Map}."$Mesh".label.gii -roi-right "$T1wFolder"/"$Subject".R."$ROI"."$Mesh".shape.gii
            ${CARET7DIR}/wb_command -set-map-names "$T1wFolder"/"$Subject".${Map}."$Mesh".dlabel.nii -map 1 "$Subject"_${Map}
        fi
    done
done

STRINGII=""
for LowResMesh in ${LowResMeshes} ; do
    STRINGII=$(echo "${STRINGII}${AtlasSpaceFolder}/fsaverage_LR${LowResMesh}k@${AtlasSpaceFolder}/fsaverage_LR${LowResMesh}k@${LowResMesh}k_fs_LR ${T1wFolder}/fsaverage_LR${LowResMesh}k@${AtlasSpaceFolder}/fsaverage_LR${LowResMesh}k@${LowResMesh}k_fs_LR ")
done

#Add CIFTI Maps to Spec Files
for STRING in "$T1wFolder"/"$NativeFolder"@"$AtlasSpaceFolder"/"$NativeFolder"@native "$AtlasSpaceFolder"/"$NativeFolder"@"$AtlasSpaceFolder"/"$NativeFolder"@native "$AtlasSpaceFolder"@"$AtlasSpaceFolder"@"$HighResMesh"k_fs_LR ${STRINGII} ; do
    FolderI=$(echo $STRING | cut -d "@" -f 1)
    FolderII=$(echo $STRING | cut -d "@" -f 2)
    Mesh=$(echo $STRING | cut -d "@" -f 3)
    for STRINGII in sulc@dscalar thickness@dscalar curvature@dscalar aparc@dlabel aparc.a2009s@dlabel ; do #Remove BA@dlabel because it doesn't convert properly
        Map=$(echo $STRINGII | cut -d "@" -f 1)
        Ext=$(echo $STRINGII | cut -d "@" -f 2)
        if [ -e "$FolderII"/"$Subject"."$Map"."$Mesh"."$Ext".nii ] ; then
            ${CARET7DIR}/wb_command -add-to-spec-file "$FolderI"/"$Subject"."$Mesh".wb.spec INVALID "$FolderII"/"$Subject"."$Map"."$Mesh"."$Ext".nii
        fi
    done
done
