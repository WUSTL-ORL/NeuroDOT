%% Read_NIFTI_Header
% Reads header of |nii| structure into 4dfp.
%
%% Description
% |header = READ_NIFTI_HEADER(nii)| reads the structure |nii| and converts
% it to a 4dfp format header |header|.
% 
% NOTE: The math in this function to calculate the |center| field for 4dfp
% is not straightforward, and is based on source code shared by Tim Coalson
% from his |nifti-4dfp| toolbox.
%
%% See Also
% <LoadVolumetricData_help.html LoadVolumetricData>