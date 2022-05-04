%% SaveVolumetricData
% Saves a volumetric data file.
%
%% Description
% |SaveVolumetricData(volume, header, filename, pn, file_type)| saves
% volumetric data defined by |volume| and |header| into a file specified by
% |filename|, path |pn|, and |file_type|.
%
% |SaveVolumetricData(volume, header, filename)| supports a full filename
% input, as long as the extension is included in the file name and matches
% a supported file type.
% 
% Supported File Types/Extensions: |'.4dfp'| 4dfp, |'.nii'| NIFTI.
% 
% NOTE: This function uses the NIFTI_Reader toolbox available on MATLAB
% Central. This toolbox has been included with NeuroDOT 2.
% 
%% Dependencies
% <Write_4dfp_Header_help.html Write_4dfp_Header> | <Write_NIFTI_help.html
% Write_NIFTI>
%
%% See Also
% <LoadVolumetricData_help.html LoadVolumetricData> |
% <Make_NativeSpace_4dfp_help.html Make_NativeSpace_4dfp>