%% smooth_Amat
% Performs Gaussian smoothing on a sensitivity matrix.
%
%% Description
%   |iA_out = SMOOTH_AMAT(iA_in, dim, gsigma)| takes the inverted VOX x MEAS
%   sensitivity matrix |iA_in| and performs Gaussian smoothing in the 3D
%   voxel space on each of the concatenated measurement matrices within,
%   returning it as |iA_out|. The user specifies the Gaussian filter
%   width "gsigma".
%
% 
%% See Also
% <Tikhonov_invert_Amat_help.html Tikhonov_invert_Amat> |
% <reconstruct_img_help.html reconstruct_img> | <FindGoodMeas_help.html
% FindGoodMeas>