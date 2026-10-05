#!/bin/bash
set -e
set -x

OUTPUTROOTDIR=fs_LR_output_directory
mkdir -p "$OUTPUTROOTDIR"
InitialMeshDirectory=InitialMesh
InputFreesurferSubjectDirectory=$1
Subject=$2
FREESURFERDIR="${3}" # path to FreeSurfer installation (e.g. /usr/local/freesurfer)
AtlasSpaceFolder="$4" # standard_mesh_atlases folder
WorkbenchDir="$5"

ReplaceOutputSubjectDirectory=true
Species="Human"

# Clean previous output directory for this subject if requested
if [ "$ReplaceOutputSubjectDirectory" = true ] ; then
  rm -rf "$OUTPUTROOTDIR/$Subject"
fi

# Create directory structure
mkdir -p "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory"
mkdir -p "$OUTPUTROOTDIR/$Subject/fsaverage"

# ------------------------------------------------------------------------------
# 1. Volume Processing (FreeSurfer mri_convert)
# ------------------------------------------------------------------------------

# Axialize orig to LPI orientation
"$FREESURFERDIR/bin/mri_convert" "$InputFreesurferSubjectDirectory/mri/orig.mgz" \
  --reslice_like "$AtlasSpaceFolder/grid_lpi.nii.gz" \
  -rt nearest "$OUTPUTROOTDIR/$Subject/orig_lpi.nii.gz"

# Affine transform orig to MNI space
"$FREESURFERDIR/bin/mri_convert" "$OUTPUTROOTDIR/$Subject/orig_lpi.nii.gz" \
  --apply_transform "$InputFreesurferSubjectDirectory/mri/transforms/talairach.xfm" \
  -oc 0 0 0 "$OUTPUTROOTDIR/$Subject/orig_mni.nii.gz"

# Generate c_ras 4x4 translation matrix for surface alignment
MatrixX=$(mri_info --cras "$InputFreesurferSubjectDirectory/mri/orig.mgz" | cut -f1 -d' ')
MatrixY=$(mri_info --cras "$InputFreesurferSubjectDirectory/mri/orig.mgz" | cut -f2 -d' ')
MatrixZ=$(mri_info --cras "$InputFreesurferSubjectDirectory/mri/orig.mgz" | cut -f3 -d' ')

CRAS_MAT="$OUTPUTROOTDIR/$Subject/cras.mat"
echo "1 0 0 $MatrixX" > "$CRAS_MAT"
echo "0 1 0 $MatrixY" >> "$CRAS_MAT"
echo "0 0 1 $MatrixZ" >> "$CRAS_MAT"
echo "0 0 0 1" >> "$CRAS_MAT"

# Extract talairach affine 4x4 matrix
MNI_MAT="$OUTPUTROOTDIR/$Subject/mni.mat"
tail -3 "$InputFreesurferSubjectDirectory/mri/transforms/talairach.xfm" | tr -d ';' > "$MNI_MAT"
echo "0.0 0.0 0.0 1.0" >> "$MNI_MAT"

# ------------------------------------------------------------------------------
# 2. Loop Through Hemispheres
# ------------------------------------------------------------------------------

for Hemisphere in L R ; do
  if [ "$Hemisphere" = "L" ] ; then 
    hemisphere="l"
    CaretStructure="CORTEX_LEFT"
  else 
    hemisphere="r"
    CaretStructure="CORTEX_RIGHT"
  fi

  # --- A. Convert Native FS Surfaces to GIFTI ---
  mris_convert "$InputFreesurferSubjectDirectory/surf/${hemisphere}h.white" \
    "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.white_orig.initial_mesh.surf.gii"

  mris_convert "$InputFreesurferSubjectDirectory/surf/${hemisphere}h.pial.T1" \
    "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.pial_orig.initial_mesh.surf.gii"

  $WorkbenchDir/wb_command -set-structure "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.white_orig.initial_mesh.surf.gii" "$CaretStructure"
  $WorkbenchDir/wb_command -set-structure "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.pial_orig.initial_mesh.surf.gii" "$CaretStructure"

  # Apply c_ras offset
  $WorkbenchDir/wb_command -surface-apply-affine \
    "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.white_orig.initial_mesh.surf.gii" \
    "$CRAS_MAT" \
    "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.white_orig.initial_mesh.surf.gii"

  $WorkbenchDir/wb_command -surface-apply-affine \
    "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.pial_orig.initial_mesh.surf.gii" \
    "$CRAS_MAT" \
    "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.pial_orig.initial_mesh.surf.gii"

  # Generate midthickness surface
  $WorkbenchDir/wb_command -surface-average \
    "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.midthickness_orig.initial_mesh.surf.gii" \
    -surf "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.pial_orig.initial_mesh.surf.gii" \
    -surf "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.white_orig.initial_mesh.surf.gii"

  # Apply MNI affine transformation
  for Surface in white midthickness pial ; do
    $WorkbenchDir/wb_command -surface-apply-affine \
      "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.${Surface}_orig.initial_mesh.surf.gii" \
      "$MNI_MAT" \
      "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.${Surface}_mni.initial_mesh.surf.gii"
  done

  # --- B. Convert Spherical Surfaces ---
  for Surface in sphere.reg sphere ; do
    mris_convert "$InputFreesurferSubjectDirectory/surf/${hemisphere}h.${Surface}" \
      "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.${Surface}.initial_mesh.surf.gii"
    
    $WorkbenchDir/wb_command -set-structure "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.${Surface}.initial_mesh.surf.gii" "$CaretStructure"
    
    $WorkbenchDir/wb_command -surface-modify-sphere \
      "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.${Surface}.initial_mesh.surf.gii" \
      100 \
      "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.${Surface}.initial_mesh.surf.gii" \
      -recenter
  done

  # --- C. Convert Scalar Features (sulc, thickness, curv) ---
  mris_convert -c "$InputFreesurferSubjectDirectory/surf/${hemisphere}h.sulc" \
    "$InputFreesurferSubjectDirectory/surf/${hemisphere}h.white" \
    "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.sulc.initial_mesh.shape.gii"
  $WorkbenchDir/wb_command -set-structure "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.sulc.initial_mesh.shape.gii" "$CaretStructure"

  mris_convert -c "$InputFreesurferSubjectDirectory/surf/${hemisphere}h.thickness" \
    "$InputFreesurferSubjectDirectory/surf/${hemisphere}h.white" \
    "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.thickness.initial_mesh.shape.gii"
  $WorkbenchDir/wb_command -set-structure "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.thickness.initial_mesh.shape.gii" "$CaretStructure"
  
  # Ensure non-negative thickness
  $WorkbenchDir/wb_command -metric-math "abs(x)" \
    "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.thickness.initial_mesh.shape.gii" \
    -var x "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.thickness.initial_mesh.shape.gii"

  mris_convert -c "$InputFreesurferSubjectDirectory/surf/${hemisphere}h.curv" \
    "$InputFreesurferSubjectDirectory/surf/${hemisphere}h.white" \
    "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.curvature.initial_mesh.shape.gii"
  $WorkbenchDir/wb_command -set-structure "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.curvature.initial_mesh.shape.gii" "$CaretStructure"

  # --- D. Generate Inflated Surfaces ---
  $WorkbenchDir/wb_command -surface-generate-inflated \
    "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.midthickness_mni.initial_mesh.surf.gii" \
    "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.inflated.initial_mesh.surf.gii" \
    "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.very_inflated.initial_mesh.surf.gii" \
    -iterations-scale 2.5

  # --- E. Prepare 164k Atlas Spheres ---
  cp "$AtlasSpaceFolder/fs_$Hemisphere/fsaverage.$Hemisphere.sphere.164k_fs_$Hemisphere.surf.gii" \
     "$OUTPUTROOTDIR/$Subject/fsaverage/$Subject.$Hemisphere.sphere.164k_fs_$Hemisphere.surf.gii"

  cp "$AtlasSpaceFolder/fs_$Hemisphere/fs_$Hemisphere-to-fs_LR_fsaverage.${Hemisphere}_LR.spherical_std.164k_fs_$Hemisphere.surf.gii" \
     "$OUTPUTROOTDIR/$Subject/fsaverage/$Subject.$Hemisphere.def_sphere.164k_fs_$Hemisphere.surf.gii"

  $WorkbenchDir/wb_command -surface-modify-sphere \
    "$OUTPUTROOTDIR/$Subject/fsaverage/$Subject.$Hemisphere.def_sphere.164k_fs_$Hemisphere.surf.gii" \
    100 \
    "$OUTPUTROOTDIR/$Subject/fsaverage/$Subject.$Hemisphere.def_sphere.164k_fs_$Hemisphere.surf.gii" \
    -recenter

  cp "$AtlasSpaceFolder/fsaverage.${Hemisphere}_LR.spherical_std.164k_fs_LR.surf.gii" \
     "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.sphere.164k_fs_LR.surf.gii"

  $WorkbenchDir/wb_command -surface-modify-sphere \
    "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.sphere.164k_fs_LR.surf.gii" \
    100 \
    "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.sphere.164k_fs_LR.surf.gii" \
    -recenter

  # --- F. Project/Unproject Spherical Surface ---
  $WorkbenchDir/wb_command -surface-sphere-project-unproject \
    "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.sphere.reg.initial_mesh.surf.gii" \
    "$OUTPUTROOTDIR/$Subject/fsaverage/$Subject.$Hemisphere.sphere.164k_fs_$Hemisphere.surf.gii" \
    "$OUTPUTROOTDIR/$Subject/fsaverage/$Subject.$Hemisphere.def_sphere.164k_fs_$Hemisphere.surf.gii" \
    "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.sphere.reg.reg_LR.initial_mesh.surf.gii"

  # --- G. Resample Native Data to 164k_fs_LR Space ---
  for Space in orig mni ; do
    for Surface in white midthickness pial ; do
      $WorkbenchDir/wb_command -surface-resample \
        "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.${Surface}_${Space}.initial_mesh.surf.gii" \
        "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.sphere.reg.reg_LR.initial_mesh.surf.gii" \
        "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.sphere.164k_fs_LR.surf.gii" \
        BARYCENTRIC \
        "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.${Surface}_${Space}.164k_fs_LR.surf.gii"
    done
  done

  for Feature in curvature sulc thickness ; do
    $WorkbenchDir/wb_command -metric-resample \
      "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.$Feature.initial_mesh.shape.gii" \
      "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.sphere.reg.reg_LR.initial_mesh.surf.gii" \
      "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.sphere.164k_fs_LR.surf.gii" \
      ADAP_BARY_AREA \
      "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.$Feature.164k_fs_LR.shape.gii" \
      -area-surfs \
        "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.midthickness_orig.initial_mesh.surf.gii" \
        "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.midthickness_orig.164k_fs_LR.surf.gii"
  done

  for Surface in inflated very_inflated ; do
    $WorkbenchDir/wb_command -surface-resample \
      "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.${Surface}.initial_mesh.surf.gii" \
      "$OUTPUTROOTDIR/$Subject/$InitialMeshDirectory/$Subject.$Hemisphere.sphere.reg.reg_LR.initial_mesh.surf.gii" \
      "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.sphere.164k_fs_LR.surf.gii" \
      BARYCENTRIC \
      "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.${Surface}.164k_fs_LR.surf.gii"
  done

  # ------------------------------------------------------------------------------
  # --- H. Downsample from 164k_fs_LR to 32k_fs_LR ---
  # ------------------------------------------------------------------------------

  Sphere32kAtlas="$AtlasSpaceFolder/${Hemisphere}.sphere.32k_fs_LR.surf.gii"

  if [ ! -f "$Sphere32kAtlas" ] ; then
    echo "ERROR: Expected atlas sphere file not found at: $Sphere32kAtlas" >&2
    exit 1
  fi

  TargetSphere32k="$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.sphere.32k_fs_LR.surf.gii"
  cp "$Sphere32kAtlas" "$TargetSphere32k"

  $WorkbenchDir/wb_command -set-structure "$TargetSphere32k" "$CaretStructure"

  $WorkbenchDir/wb_command -surface-modify-sphere \
    "$TargetSphere32k" \
    100 \
    "$TargetSphere32k" \
    -recenter

  Sphere164k="$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.sphere.164k_fs_LR.surf.gii"

  # Resample anatomical surfaces (164k -> 32k)
  for Space in orig mni ; do
    for Surface in white midthickness pial ; do
      $WorkbenchDir/wb_command -surface-resample \
        "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.${Surface}_${Space}.164k_fs_LR.surf.gii" \
        "$Sphere164k" \
        "$TargetSphere32k" \
        BARYCENTRIC \
        "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.${Surface}_${Space}.32k_fs_LR.surf.gii"
    done
  done

  # Resample scalar metrics (164k -> 32k)
  for Feature in curvature sulc thickness ; do
    $WorkbenchDir/wb_command -metric-resample \
      "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.$Feature.164k_fs_LR.shape.gii" \
      "$Sphere164k" \
      "$TargetSphere32k" \
      ADAP_BARY_AREA \
      "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.$Feature.32k_fs_LR.shape.gii" \
      -area-surfs \
        "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.midthickness_orig.164k_fs_LR.surf.gii" \
        "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.midthickness_orig.32k_fs_LR.surf.gii"
  done

  # Resample inflated surfaces (164k -> 32k)
  for Surface in inflated very_inflated ; do
    $WorkbenchDir/wb_command -surface-resample \
      "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.${Surface}.164k_fs_LR.surf.gii" \
      "$Sphere164k" \
      "$TargetSphere32k" \
      BARYCENTRIC \
      "$OUTPUTROOTDIR/$Subject/$Subject.$Hemisphere.${Surface}.32k_fs_LR.surf.gii"
  done

done

# Clean temporary matrices
rm -f "$CRAS_MAT" "$MNI_MAT"

exit 0