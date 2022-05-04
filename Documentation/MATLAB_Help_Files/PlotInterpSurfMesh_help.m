%% PlotInterpSurfMesh
% Interpolates volumetric data onto hemispheric meshes for display.
% 
%% Description
% |PlotInterpSurfMesh(volume, meshL, meshR, dim)| takes functional imaging
% data |volume| and interpolates it onto the left and right hemispheric
% surface meshes given in |meshL| and |meshR|, respectively, using the
% spatial information in |dim|. The result is a bilateral sagittal view of
% the activations overlain onto the surface of the brain represented by the
% meshes.
%
% |PlotInterpSurfMesh(volume, meshL, meshR, dim, params)| allows the user
% to specify parameters for plot creation.
% 
%% Visualization Parameters
% |params| fields that apply to this function (and their defaults):
%
% <html>
% <table border = 1>
% <tr><td>Name</td><td>Default</td><td>Effect</td></tr>
% <tr><td>Scale</td><td>(90% max)</td><td>Maximum value to which image is
% scaled.</td><tr>
% <tr><td>Th</td><td></td><td>Thresholds.</td></tr>
% <tr><td>.P</td><td>(25% max)</td><td>Value of min threshold to display
% positive data values.</td></tr>
% <tr><td>.N</td><td>(-Th.P)</td><td>Value of max threshold to display
% negative data values.</td></tr>
% </table></html>
% 
%% Dependencies
% <vol2surf_mesh_help.html vol2surf_mesh> | <PlotLRMeshes_help.html
% PlotLRMeshes> | <adjust_brain_pos_help.html adjust_brain_pos> |
% <rotate_cap_help.html rotate_cap> | <rotation_matrix_help.html
% rotation_matrix>
%
%% See Also
% <PlotSlices_help.html PlotSlices>