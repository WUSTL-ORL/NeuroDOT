%% affine3d_img
% Transforms a 3D data set to a new space.
% 
%% Description
% |imgB = affine3d_img(imgA, infoA, infoB, affine)| takes a reconstructed,
% VOX x TIME image |imgA| and transforms it from its initial voxel space
% defined by the structure |infoA| into a target voxel space defined by the
% structure |infoB| and using the transform matrix |affine|. The output is
% a VOX x TIME matrix |imgB| in the target voxel space.
% 
% |imgB = affine3d_img(imgA, infoA, infoB, affine, interp_type)| allows the
% user to specify an interpolation method for the |interp3| function that
% |affine3d_img| uses. Other methods that can be used (input as strings)
% are |'nearest'|, |'spline'|, and |'cubic'|. The default value is
% |'linear'|.
% 
%% See Also
% <spectroscopy_img_help.html spectroscopy_img> |
% <change_space_coords_help.html change_space_coords> |
% <matlab:doc('interp3') interp3>