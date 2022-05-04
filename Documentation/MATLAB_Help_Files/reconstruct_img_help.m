%% reconstruct_img
% Performs image reconstruction by wavelength using inverted A-matrix.
%
%% Description
% |img = reconstruct_img(data, iA)| takes the inverted VOX x MEAS
% sensitivity matrix |iA| and right-multiplies it by the preprocessed MEAS
% x TIME light-level matrix |data| to reconstruct an image in voxel space.
% The image is output in a VOX x TIME matrix |img|.
%
%% See Also
% <Tikhonov_invert_Amat_help.html Tikhonov_invert_Amat> |
% <smooth_Amat_help.html smooth_Amat> | <spectroscopy_img_help.html
% spectroscopy_img> | <FindGoodMeas_help.html FindGoodMeas>