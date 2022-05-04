%% PlotSlicesMov
% Creates an interactive three-slice plot.
%
%% Description
% |PlotSlicesMov(underlay, infoVol, params, overlay, fn)| uses the same
% basic inputs as the function |PlotSlices| to create a video of a 4D
% volume, saved under the filename |fn|.
% 
% |PlotSlicesMov(..., infoFR)| allows the user to input an |info| structure
% or a framerate for the video to be recorded with. The default framerate
% is |1 Hz|.
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
% </table></html>
%
% Note: |applycmap| has further options for using |params| to specify
% parameters for the fusion, scaling, and colormapping process.
%
%% Dependencies
% <applycmap_help.html applycmap> | <PlotSlices_help.html PlotSlices>