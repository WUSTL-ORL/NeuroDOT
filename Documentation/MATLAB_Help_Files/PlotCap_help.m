%% PlotCap
% A visualization of the cap grid.
% 
%% Description
% |PlotCap(info)| plots a basic diagram of the cap based on the cap
% metadata stored in |info.optodes|. All optodes are numbered; sources are
% colored light red, and detectors light blue.
% 
% |PlotCap(info, params)| allows the user to specify parameters for plot
% creation.
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
% </table></html>
% 
%% Dependencies
% <PlotCapData_help.html PlotCapData>
%
%% See Also
% <PlotCapGoodMeas_help.html PlotCapGoodMeas> | <PlotCapMeanLL_help.html
% PlotCapMeanLL>