%% PlotGrayData
% A basic scaled image plotting function.
%
%% Description
% |PlotGrayData(data)| plots the input |data| as a scaled image.
%
% |PlotGrayData(data, params)| allows the user to specify parameters for
% plot creation.
%
%% Visualization Parameters
% |params| fields that apply to this function (and their defaults):
%
% <html>
% <table border = 1>
% <tr><td>Name</td><td>Default</td><td>Effect</td></tr>
% <tr><td>fig_size</td><td>[200, 200, 560, 420]</td><td>Default figure
% position vector.</td></tr>
% <tr><td>fig_handle</td><td>(none)</td><td>Specifies a figure to target.
% If empty, spawns a new figure.</td></tr>
% <tr><td>climits</td><td>[data min and max]</td><td>Limits of color
% axis.</td></tr>
% </table></html>
%
%% See Also
% <PlotGray_help.html PlotGray>