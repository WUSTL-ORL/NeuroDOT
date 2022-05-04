%% LoadVolumetricData
% Loads a volumetric data file.
%
%% Description
% |[volume, header] = LoadVolumetricData(filename, pn, file_type)| loads a
% file specified by |filename|, path |pn|, and |file_type|, and returns it
% in two parts: the raw data file |volume|, and the header |header|,
% containing a number of key-value pairs in 4dfp format.
%
% |[volume, header] = LoadVolumetricData(filename)| supports a full
% filename input, as long as the extension is included in the file name and
% matches a supported file type.
% 
% Supported File Types/Extensions: |'.4dfp'| 4dfp, |'.nii'| NIFTI.
% 
% NOTE: This function uses the NIFTI_Reader toolbox available on MATLAB
% Central. This toolbox has been included with NeuroDOT 2.
% 
%% Dependencies
% <Read_4dfp_Header_help.html Read_4dfp_Header> |
% <Read_NIFTI_Header_help.html Read_NIFTI_Header> |
% <Make_NativeSpace_4dfp_help.html Make_NativeSpace_4dfp>
%
%% See Also
% <SaveVolumetricData_help.html SaveVolumetricData>