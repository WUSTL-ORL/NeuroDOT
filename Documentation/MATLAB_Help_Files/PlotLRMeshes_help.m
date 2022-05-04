%% PlotLRMeshes
% Renders a pair of hemispheric meshes.
%
%% Description
% |PlotLRMeshes(meshL, meshR)| renders the data in a pair of left and right
% hemispheric meshes |meshL.data| and |meshR.data|, respectively, and
% applies full color mapping to them. If no data is present for either
% mesh, a default gray mesh will be plotted.
%
% |PlotLRMeshes(meshL, meshR, params)| allows the user to specify
% parameters for plot creation.
% 
%% Visualization Parameters
% |params| fields that apply to this function (and their defaults):
%
% <html>
% <table border = 1>
% <tr><td>Name</td><td>Default</td><td>Effect</td></tr>
% <tr><td>fig_size</td><td>[20, 200, 960, 420]</td><td>Default figure
% position vector.</td></tr>
% <tr><td>fig_handle</td></td>(none)</td><td>Specifies a figure to
% target. If empty, spawns a new figure.</td></tr>
% <tr><td>Scale</td><td>(90% max)</td><td>Maximum value to which image is
% scaled.</td></tr>
% <tr><td>PD</td><td>0</td><td>Declares that input image is positive
% definite.</td></tr>
% <tr><td>cblabels</td><td>([-90% max, 90% max])</td><td>Colorbar axis
% labels. Min defaults to 0 if PD==1, both default to +/- Scale if
% supplied.</td></tr>
% <tr><td>cbticks</td><td>(none)</td><td>Specifies positions of tick marks
% on colorbar axis.</td></tr>
% <tr><td>alpha</td><td>1</td><td>Transparency of mesh.</td></tr>
% <tr><td>view</td><td>'lat'</td><td>Sets the view perspective.</td></tr>
% </table></html>
%
%% Dependencies
% <adjust_brain_pos_help.html adjust_brain_pos> | <rotate_cap_help.html
% rotate_cap> | <rotation_matrix_help.html rotation_matrix> |
% <applycmap_help.html applycmap> 
% 
%% See Also
% <PlotInterpSurfMesh_help.html PlotInterpSurfMesh> |
% <vol2surf_mesh_help.html vol2surf_mesh>