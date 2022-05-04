%% adjust_brain_pos
% Repositions mesh orientations for display.
%
%% Description
% |[Lnodes, Rnodes] = adjust_brain_pos(meshL, meshR)| takes the left and
% right hemispheric meshes |meshL| and |meshR|, respectively, and
% repositions them to the proper perspective for display.
%
% |[Lnodes, Rnodes] = adjust_brain_pos(meshL, meshR, params)| allows the
% user to specify parameters for plot creation.
%
%% Visualization Parameters
% |params| fields that apply to this function (and their defaults):
%
% <html>
% <table border = 1>
% <tr><td>Name</td><td>Default</td><td>Effect</td></tr>
% <tr><td>ctx</td><td>'std'</td><td>Defines inflation of mesh.</td></tr>
% <tr><td>orientation</td><td>'t'</td><td>Select orientation of volume. 't'
% for transverse, 's' for sagittal.</td></tr>
% <tr><td>view</td><td>'lat'</td><td>Sets the view perspective.</td></tr>
% </table></html>
% 
%% Dependencies
% <rotate_cap_help.html rotate_cap> | <rotation_matrix_help.html
% rotation_matrix>
% 
%% See Also
% <PlotLRMeshes_help.html PlotLRMeshes>