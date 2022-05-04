%% Write_NIFTI
% Assembles an "nii" structure from 4dfp data.
% 
%% Description
% |nii = Write_NIFTI(volume, header, filename)| converts the volumetric
% data |volume| described by |header| from 4dfp to NIFTI format, output as
% |nii|.
% 
% NOTE: The math in this function to calculate the |center| field for 4dfp
% is not straightforward, and is based on source code shared by Tim Coalson
% from his |nifti-4dfp| toolbox.
%
%% See Also
% <SaveVolumetricData_help.html SaveVolumetricData>
