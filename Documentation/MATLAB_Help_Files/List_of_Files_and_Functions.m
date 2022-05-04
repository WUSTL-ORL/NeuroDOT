%% NeuroDOT 2.2.0 Base Edition
% *List of Files and Functions*

%% Tutorials
% <html><table border=0>
% <tr><td><a href='DOT_Processing_Tutorial.html'>DOT Processing
% Tutorial</a></td><td>Covers the DOT Processing pipeline.</td></tr>
% </table></html>

%% Data
% <html><table border=0>
% <tr><td>/NeuroDOT_Base_AC_Sample.mat</td><td>Based on the data set from
% subject 1409, date 161110, acquisition AC002.</td></tr>
% <tr><td>/NeuroDOT_Base_CCW_Sample.mat</td><td>Based on the data set from
% subject 0813, date 101019, acquisition CCW.</td></tr>
% <tr><td>/NeuroDOT_Base_CW_Sample.mat</td><td>Based on the data set from
% subject 0813, date 101019, acquisition CW.</td></tr>
% <tr><td>/NeuroDOT_Base_GV_Sample.mat</td><td>Based on the data set from
% subject 1106, date 121023, acquisition GV.</td></tr>
% <tr><td>/NeuroDOT_Base_HW_Sample_1.mat</td><td>Based on the data set from
% subject 1106, date 120702, acquisition HW.</td></tr>
% <tr><td>/NeuroDOT_Base_HW_Sample_2.mat</td><td>Based on the data set from
% subject 1205, date 120626, acquisition HW.</td></tr>
% <tr><td>/NeuroDOT_Base_HW_Sample_Noisy.mat</td><td>Based on the data set
% from subject 1531, date 150810, acquisition HW001.</td></tr>
% <tr><td>/NeuroDOT_Base_IN_Sample.mat</td><td>Based on the data set from
% subject 0813, date 101019, acquisition IN.</td></tr>
% <tr><td>/NeuroDOT_Base_OUT_Sample.mat</td><td>Based on the data set from
% subject 0813, date 101019, acquisition OUT.</td></tr>
% <tr><td>/NeuroDOT_Base_RW_Sample.mat</td><td>Based on the data set from
% subject 1106, date 120702, acquisition RW.</td></tr>
% <tr><td>/NeuroDOT_Base_Converter_Sample.mat</td><td>This sample is
% provided for use in the File IO tutorial appendix, based on HW Sample
% 1.</td></tr>
% <tr><td>/test_homer2.nirs</td><td>This sample is provided for use in the
% File IO tutorial appendix. It was derived from the introductory sample
% used in the HOMER2 beginners’ guide.</td></tr>
% </table></html>

%% Functions
%%% _Analysis_
% <html><table border=0>
% <tr><td><a
% href='BlockAverage_help.html'>BlockAverage</a></td><td>Averaging of
% stimulus blocks.</td></tr>
% <tr><td><a href='FindGoodMeas_help.html'>FindGoodMeas</a></td><td>Filter
% out noisy channels.</td></tr>
% </table></html>

%%% _File IO_
% <html><table border=0>
% <tr><td><a
% href='Check4MissingData_help.html'>Check4MissingData</a></td><td>
% Validates framesynch file against key file info.</td></tr>
% <tr><td><a
% href='converter_data_help.html'>converter_data</a></td><td>Reshapes data
% from ND1 to ND2 format.</td></tr>
% <tr><td><a
% href='converter_HOMER_to_ND2_help.html'>converter_HOMER_to_ND2</a></td><td>Converts
% workspace from HOMER to ND2.</td></tr>
% <tr><td><a
% href='converter_info_help.html'>converter_info</a></td><td>Converts the
% ND1 info structure to ND2 format.</td></tr>
% <tr><td><a
% href='converter_ND1_to_ND2_help.html'>converter_ND1_to_ND2</a></td><td>Converts
% both data and info from ND1 format to ND2.</td></tr>
% <tr><td><a
% href='converter_ND2_to_HOMER_help.html'>converter_ND2_to_HOMER</a></td><td>Converts
% workspace from ND2 to HOMER.</td></tr>
% <tr><td><a
% href='converter_ND2_to_ND1_help.html'>converter_ND2_to_ND1</a></td><td>Converts
% both data and info from ND2 format to ND1.</td></tr>
% <tr><td><a href='Crop2Synch_help.html'>Crop2Synch</a></td><td>Crops data
% to synch points.</td></tr>
% <tr><td><a
% href='InterpretPulses_help.html'>InterpretPulses</a></td><td>Reads the
% synch pulses into an experiment paradigm.</td></tr>
% <tr><td><a
% href='InterpretStimSynch_help.html'>InterpretStimSynch</a></td><td>Reads
% stimulus from frame synch file.</td></tr>
% <tr><td><a
% href='InterpretSynchBeeps_help.html'>InterpretSynchBeeps</a></td><td>Automated
% synch file analysis.</td></tr>
% <tr><td><a
% href='Load_AcqDecode_Data_help.html'>Load_AcqDecode_Data</a></td><td>Loads
% a single AcqDecode file.</td></tr>
% <tr><td><a href='Load_HOMER_help.html'>Load_HOMER</a></td><td>Converts
% and loads a HOMER2 ".nirs" file into ND2.</td></tr>
% <tr><td><a
% href='LoadMulti_AcqDecode_Data_help.html'>LoadMulti_AcqDecode_Data</a></td><td>Loads
% and combines data from multiple AcqDecode files for a single
% scan.</td></tr>
% <tr><td><a
% href='LoadVolumetricData_help.html'>LoadVolumetricData</a></td><td>Loads
% a volumetric data file.</td></tr>
% <tr><td><a
% href='Make_NativeSpace_4dfp_help.html'>Make_NativeSpace_4dfp</a></td><td>Calculates
% native space of incomplete 4dfp header.</td></tr>
% <tr><td><a
% href='Read_4dfp_Header_help.html'>Read_4dfp_Header</a></td><td>Read the
% .ifh header of a 4dfp file.</td></tr>
% <tr><td><a
% href='Read_AcqDecode_Header_help.html'>Read_AcqDecode_Header</a></td><td>Interprets
% the header of an AcqDecode generated data file.</td></tr>
% <tr><td><a
% href='Read_NIFTI_Header_help.html'>Read_NIFTI_Header</a></td><td>Reads a
% NIFTI file header.</td></tr>
% <tr><td><a href='ReadAux_help.html'>ReadAux</a></td><td>Reads auxiliary
% channel files.</td></tr>
% <tr><td><a href='ReadInfoTxt_help.html'>ReadInfoTxt</a></td><td>Reads a
% key file.</td></tr>
% <tr><td><a href='Save_HOMER_help.html'>Save_HOMER</a></td><td>Converts an
% ND2 workspace to HOMER2 ".nirs" format to save.</td></tr>
% <tr><td><a
% href='SaveVolumetricData_help.html'>SaveVolumetricData</a></td><td>Saves
% volumetric data in supported formats.</td></tr>
% <tr><td><a
% href='Write_4dfp_Header_help.html'>Write_4dfp_Header</a></td><td>Writes a
% 4dfp header.</td></tr>
% <tr><td><a href='Write_NIFTI_help.html'>Write_NIFTI</a></td><td>Prepares
% a "nii" structure to be written to file.</td></tr>
% </table></html>

%%% _Reconstruction_
% <html><table border=0>
% <tr><td><a
% href='reconstruct_img_help.html'>reconstruct_img</a></td><td>Performs
% image reconstruction by wavelength using inverted A-matrix.</td></tr>
% <tr><td><a href='smooth_Amat_help.html'>smooth_Amat</a></td><td>Performs
% Gaussian smoothing on a sensitivity matrix.</td></tr>
% <tr><td><a
% href='spectroscopy_img_help.html'>spectroscopy_img</a></td><td>Completes
% the Beer-Lambert law from a reconstructed image.</td></tr>
% <tr><td><a
% href='Tikhonov_invert_Amat_help.html'>Tikhonov_invert_Amat</a></td><td>Inverts
% a sensitivity matrix.</td></tr>
% </table></html>

%%% _Spatial Transforms_
% <html><table border=0>
% <tr><td><a
% href='affine3d_img_help.html'>affine3d_img</a></td><td>Transforms a 3D
% data set to a new space.</td></tr>
% <tr><td><a
% href='change_space_coords_help.html'>change_space_coords</a></td><td>Applies
% a look up to change 3D coordinates into a new space.</td></tr>
% <tr><td><a href='Good_Vox2vol_help.html'>Good_Vox2vol</a></td><td>Turns a
% Good Voxels data stream into a volume.</td></tr>
% <tr><td><a href='rotate_cap_help.html'>rotate_cap</a></td><td>Rotates the
% cap in space.</td></tr>
% <tr><td><a
% href='rotation_matrix_help.html'>rotation_matrix</a></td><td>Creates a
% rotation matrix.</td></tr>
% </table></html>

%%% _Support_
% <html><table border=0>
% <tr><td><a
% href='CheckOrientation_help.html'>CheckOrientation</a></td><td>Interprets
% 3D orientation information.</td></tr>
% <tr><td><a href='istablevar_help.html'>istablevar</a></td><td>Same
% operation as "isfield" for tables.</td></tr>
% </table></html>

%%% _Temporal Transforms_
% <html><table border=0>
% <tr><td><a href='detrend_tts_help.html'>detrend_tts</a></td><td>Performs
% linear detrending.</td></td>
% <tr><td><a href='fft_tts_help.html'>fft_tts</a></td><td>Fast Fourier
% Transform function with normalization and multiple outputs.</td></tr>
% <tr><td><a href='gethem_help.html'>gethem</a></td><td>Generate
% hemodynamics for SSR.</td></tr>
% <tr><td><a href='highpass_help.html'>highpass</a></td><td>High pass
% filter.</td></tr>
% <tr><td><a href='logmean_help.html'>logmean</a></td><td>Negative ln of
% each channel divided by its mean.</td></tr>
% <tr><td><a href='lowpass_help.html'>lowpass</a></td><td>Low pass
% filter.</td></tr>
% <tr><td><a
% href='normalize2range_tts_help.html'>normalize2range_tts</a></td><td>Normalizes
% time traces to average of a range.</td></tr>
% <tr><td><a href='regcorr_help.html'>regcorr</a></td><td>SSR to remove
% hemodynamics from data.</td></tr>
% <tr><td><a
% href='resample_tts_help.html'>resample_tts</a></td><td>Resample data
% while maintaining linear signal component.</td></tr>
% </table></html>

%%% _Visualization_
% <html><table border=0>
% <tr><td><a
% href='adjust_brain_pos_help.html'>adjust_brain_pos</a></td><td>Repositions
% mesh orientations for display.</td></tr>
% <tr><td><a href='applycmap_help.html'>applycmap</a></td><td>Performs
% color mapping and fuses images with anatomical models.</td></tr>
% <tr><td><a href='PlotCap_help.html'>PlotCap</a></td><td>A visualization
% of the cap grid.</td></tr>
% <tr><td><a href='PlotCapData_help.html'>PlotCapData</a></td><td>A basic
% plotting function for generating and labeling cap grids.</td></tr>
% <tr><td><a
% href='PlotCapGoodMeas_help.html'>PlotCapGoodMeas</a></td><td>A Good
% Measurements visualization overlaid on a cap grid.</td></tr>
% <tr><td><a href='PlotCapMeanLL_help.html'>PlotCapMeanLL</a></td><td>A
% visualization of mean light levels overlaid on a cap grid.</td></tr>
% <tr><td><a href='PlotFalloffData_help.html'>PlotFalloffData</a></td><td>A
% basic falloff plotting function.</td></tr>
% <tr><td><a href='PlotFalloffLL_help.html'>PlotFalloffLL</a></td><td>A
% light-level falloff visualization.</td></tr>
% <tr><td><a href='PlotGray_help.html'>PlotGray</a></td><td>A visualization
% of data as a scaled image.</td></tr>
% <tr><td><a href='PlotGrayData_help.html'>PlotGrayData</a></td><td>A basic
% scaled image plotting function.</td></tr>
% <tr><td><a
% href='PlotHistogramData_help.html'>PlotHistogramData</a></td><td>A basic
% histogram plotting function.</td></tr>
% <tr><td><a
% href='PlotHistogramSTD_help.html'>PlotHistogramSTD</a></td><td>A
% histogram of channel standard deviations.</td></tr>
% <tr><td><a
% href='PlotInterpSurfMesh_help.html'>PlotInterpSurfMesh</a></td><td>Interpolates
% volumetric data onto hemispheric meshes for display.</td></tr>
% <tr><td><a href='PlotLRMeshes_help.html'>PlotLRMeshes</a></td><td>Renders
% a pair of hemispheric meshes.</td></tr>
% <tr><td><a
% href='PlotMeshSurface_help.html'>PlotMeshSurface</a></td><td>Shows a head
% mesh model's outer surface.</td></tr>
% <tr><td><a
% href='PlotPowerSpectrumAllMeas_help.html'>PlotPowerSpectrumAllMeas</a></td><td>An
% all-measurements frequency-domain time trace visualization.</td></tr>
% <tr><td><a
% href='PlotPowerSpectrumData_help.html'>PlotPowerSpectrumData</a></td><td>A
% basic frequency-domain time traces plotting function.</td></tr>
% <tr><td><a
% href='PlotPowerSpectrumMean_help.html'>PlotPowerSpectrumMean</a></td><td>A
% mean frequency-domain time trace visualization.</td></tr>
% <tr><td><a href='PlotSlices_help.html'>PlotSlices</a></td><td>Creates an
% interactive three-slice plot.</td></tr>
% <tr><td><a
% href='PlotSlicesMov_help.html'>PlotSlicesMov</a></td><td>Creates a video
% file from a volumetric data set.</td></tr>
% <tr><td><a
% href='PlotSlicesTimeTrace_help.html'>PlotSlicesTimeTrace</a></td><td>Creates
% an interactive 4D plot.</td></tr>
% <tr><td><a
% href='PlotTimeTraceAllMeas_help.html'>PlotTimeTraceAllMeas</a></td><td>An
% all-measurements time trace visualization.</td></tr>
% <tr><td><a
% href='PlotTimeTraceData_help.html'>PlotTimeTraceData</a></td><td>A basic
% time traces plotting function.</td></tr>
% <tr><td><a
% href='PlotTimeTraceMean_help.html'>PlotTimeTraceMean</a></td><td>A mean
% time traces visualization.</td></tr>
% <tr><td><a
% href='vol2surf_mesh_help.html'>vol2surf_mesh</a></td><td>Interpolates
% volumetric data onto a surface mesh.</td></tr>
% </table></html>

%% Support Files
% <html><table border=0>
% <tr><td>/Raw_Support_Files/Spectroscopy/E.mat</td><td>Contains an “E”
% matrix of spectroscopy coefficients that relate the light intensity to
% HbO and HbR concentrations via the modified Beer-Lambert law.</td></tr>
% <tr><td>/A_Adult_96x92.mat</td><td>Contains several variables required
% for image reconstruction with the ORL’s adult scanning cap: an “A” matrix
% for cap sensitivity, a “dim” structure describing the DOT space, and an
% “infoA” structure describing the atlas space.</td></tr>
% <tr><td>/Atlas_MNI152nl_T1_on_111.mat</td><td>This file contains a
% version of the MNI 152 non-linear T1 111 atlas that has been transformed
% into the DOT space. It is saved as a variable “atlas”.</td></tr>
% <tr><td>/Atlas_MNI152nl_T1_on_333.mat</td><td>This file contains a
% version of the MNI 152 non-linear T1 333 atlas that has been transformed
% into the DOT space. It is saved as a variable “atlas”.</td></tr>
% <tr><td>/Atlas_Segmented_MNI152nl_T1_on_111.mat</td><td>This file
% contains a version of the MNI 152 non-linear T1 111 atlas segmentation
% that has been transformed into the DOT space. It is saved as a variable
% “atlas”.</td></tr>
% <tr><td>/Atlas_Segmented_MNI152nl_T1_on_333.mat</td><td>This file
% contains a version of the MNI 152 non-linear T1 333 atlas segmentation
% that has been transformed into the DOT space. It is saved as a variable
% “atlas”.</td></tr>
% <tr><td>/Cap_Fitter_Mesh_Sample.mat</td><td>A full-head sample mesh
% provided for use with the Cap_Fitter GUI.</td></tr>
% <tr><td>/LR_Meshes_MNI_164k.mat</td><td>Contains two meshes, "MNIl" and
% "MNIr", corresponding to the left and right halves of the segmented MNI
% 152 non-linear T1 111 atlas.</td></tr>
% <tr><td>/Pad_Adult_96x92.mat</td><td>Contains an “info” structure with
% the “optodes”, “pairs”, and “tissue” substructures corresponding to the
% ORL’s adult scanning cap, with 96 sources and 92 detectors.</td></tr>
% <tr><td>/Pad_AdultV24x28.mat</td><td>Contains an “info” structure with
% the “optodes”, “pairs”, and “tissue” substructures corresponding to the
% ORL’s adult visual cortex cap, with 24 sources and 28
% detectors.</td></tr>
% <tr><td>/Pad_Baby32x34.mat</td><td>Contains an “info” structure with the
% “optodes”, “pairs”, and “tissue” substructures corresponding to the ORL’s
% infant scanning cap, with 32 sources and 34 detectors.</td></tr>
% </table></html>
