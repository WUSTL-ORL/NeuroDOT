%% PlotCapData
% A basic plotting function for generating and labeling cap grids.
% 
%% Description
% |PlotCapData(SrcRGB, DetRGB, info)| plots the input RGB information in
% one of three modes:
% 
% |'text'| - Optode numbers are arranged in a cap grid and colored with the
% RGB input.
% 
% |'patch'| - Optodes are plotted as patches and colored with the RGB
% input.
% 
% |'textpatch'| - Optodes are plotted as patches and colored with the RGB
% input, with optode numbers overlain in white.
%
% |PlotCapData(SrcRGB, DetRGB, info, params)| allows the user to specify
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
% <tr><td>fig_handle</td><td>(none)</td><td>Specifies a figure to target.
% If empty, spawns a new figure.</td></tr>
% <tr><td>dimension</td><td>'2D'</td><td>Specifies either a 2D or 3D plot
% rendering.</td></tr>
% <tr><td>mode</td><td>'textpatch'</td><td>Display mode.</td></tr>
% </table></html>
% 
%% Dependencies
% <PlotCapData_help.html PlotCapData>
%
%% See Also
% <PlotCapGoodMeas_help.html PlotCapGoodMeas> | <PlotCapMeanLL_help.html
% PlotCapMeanLL>