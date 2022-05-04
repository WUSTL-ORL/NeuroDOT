%% spectroscopy_img
% Completes the Beer-Lambert law from a reconstructed image.
%
%% Description
% |img_out = spectroscopy_img(img_in, E)| takes the reconstructed VOX x
% TIME x WL mua image |img_in| and multiplies it by the inverse of
% the extinction coefficient matrix |E| to create an output VOX x TIME x HB
% matrix |img_out|, where HB 1 and 2 are the voxel-space time series images
% for HbR and HbO, respectively.
%
%% See Also
% <reconstruct_img_help.html reconstruct_img> | <affine3d_img_help.html
% affine3d_img>