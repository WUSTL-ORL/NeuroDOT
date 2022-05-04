%% PlotMeshSurface
% Creates a 3D surface mesh visualization.
%
%% Description
% |PlotMeshSurface(mesh)| creates a 3D visualization of the surface mesh
% |mesh|. If no region data is provided in |mesh.region|, all nodes will
%  be assumed to form a single region. If a field |data| is provided as
%  part of the |mesh| structure, that data will be used to color the
%  visualization. If both |data| and |region| are present, the |region|
%  values are used as an underlay for the colormapping.
%   
% |PlotMeshSurface(mesh, params)| allows the user to specify
% parameters for plot creation.
%
%% Visualization Parameters
% |params| fields that apply to this function (and their defaults):
%
% <html>
% <table border = 1>
% <tr><td>Name</td><td>Default</td><td>Effect</td></tr>
% <tr><td>fig_size</td><td>[20, 200, 1240, 420]</td><td>Default figure
% position vector.</td></tr>
% <tr><td>fig_handle</td><td>(none)</td><td>Specifies a figure to
% target.</td></tr>
% <tr><td>TC</td><td>1</td><td>Direct map integer data values to defined
% color map ("True Color").</td></tr>
% <tr><td>Cmap.P</td><td>'jet'</td><td>Colormap for positive data
% values.</td></tr>
% <tr><td>BG</td><td>[0.8, 0.8, 0.8]</td><td>Background color, as an RGB
% triplet.</td></tr>
% <tr><td>orientation</td><td>'t'</td><td>Select orientation of volume. 't'
% for transverse, 's' for sagittal.</td></tr>
% </table>
% </html>
%
% Note: |applycmap| has further options for using |params| to specify
% parameters for the fusion, scaling, and colormapping process.
%
%% Dependencies
% <applycmap_help.html applycmap>
%
%% See Also
% <PlotSlices_help.html PlotSlices> | <PlotCap_help.html PlotCap> |
% <Cap_Fitter_help.html Cap_Fitter>