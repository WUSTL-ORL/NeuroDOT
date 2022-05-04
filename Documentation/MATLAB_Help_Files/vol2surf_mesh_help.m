%% vol2surf_mesh
% Interpolates volumetric data onto a surface mesh.
%
%% Description
% |Smesh = vol2surf_mesh(Smesh, volume, dim)| takes the surface mesh |Smesh|
% and interpolates the values of the volumetric data |volume| at the mesh's
% nodes, using the spatial information in |dim|. These values are
% output as |Smesh|.
%
% |Smesh = vol2surf_mesh(Smesh, volume, dim, params)| allows the user
% to specify parameters for plot creation.
%
%% Visualization Parameters
% |params| fields that apply to this function (and their defaults):
%
% <html>
% <table border = 1>
% <tr><td>Name</td><td>Default</td><td>Effect</td></tr>
% <tr><td>OL</td><td>0</td><td>If "overlap" data is presented (OL==1), this
% sets the interpolation method to "nearest". Default is
% "linear".</td></tr>
% </table></html>
%
%% See Also
% <PlotInterpSurfMesh_help.html PlotInterpSurfMesh> |
% <Good_Vox2vol_help.html Good_Vox2vol> | <affine3d_img_help.html
% affine3d_img>