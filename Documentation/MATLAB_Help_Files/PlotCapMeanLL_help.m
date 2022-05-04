%% PlotCapMeanLL
% A visualization of mean light levels overlaid on a cap grid.
%
%% Description
% |PlotCapMeanLL(data, info)| plots an intensity map of the mean light
% levels for specified measurement groupings of each optode on the cap and
% arranges them based on the metadata in |info.optodes|.
%
% |PlotCapMeanLL(data, info, params)| allows the user to specify parameters
% for plot creation.
%
%% Visualization Parameters
% |params| fields that apply to this function (and their defaults):
%
% <html>
% <table border = 1>
% <tr><td>Name</td><td>Default</td><td>Effect</td></tr>
% <tr><td>fig_size</td><td>[20, 200, 1240, 420]</td><td>Default figure
% position vector.</td></tr>
% <tr><td>fig_handle</td><td>(none)</td><td>Specifies a figure to target.
% If empty, spawns a new figure.</td></tr>
% <tr><td>dimension</td><td>'2D'</td><td>Specifies either a 2D or 3D plot
% rendering.</td></tr>
% <tr><td>rlimits</td><td>(all R2D)</td><td>Limits of pair radii
% displayed.</td></tr>
% <tr><td>Nnns</td><td>(all NNs)</td><td>Number of NNs displayed.</td></tr>
% <tr><td>Nwls</td><td>(all WLs)</td><td>Number of WLs displayed.</td></tr>
% <tr><td>useGM</td><td>0</td><td>Use Good Measurements.</td></tr>
% <tr><td>Cmap.P</td><td>'hot'</td><td>Default color mapping.</td></tr>
% </table></html>
%
%% Dependencies
% <PlotCapData_help.html PlotCapData> | <istablevar_help.html istablevar> |
% <applycmap_help.html applycmap>
% 
%% See Also
% <PlotCap_help.html PlotCap> | <PlotCapGoodMeas_help.html PlotCapGoodMeas>