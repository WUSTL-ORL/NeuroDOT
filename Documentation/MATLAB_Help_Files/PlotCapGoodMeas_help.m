%% PlotCapGoodMeas
% A Good Measurements visualization overlaid on a cap grid.
%
%% Description
% |PlotCapGoodMeas(info)| plots a visualization of the noisy channels
% determined by |FindGoodMeas| and arranges them based on the metadata in
% |info.optodes|. Good channels are depicted as green lines, bad channels
% red lines; sources and detectors are given lettering in light blue and
% red.
%
% The plot title provides tallies for all specified groupings. The next
% line of the title lists how many optodes for which only 33% of their
% measurements are good. These optodes are surrounded with white circles.
%
% |PlotCapGoodMeas(info, params)| allows the user to specify
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
% <tr><td>rlimits</td><td>(all R2D)</td><td>Limits of pair radii
% displayed.</td></tr>
% <tr><td>Nnns</td><td>(all NNs)</td><td>Number of NNs displayed.</td></tr>
% <tr><td>Nwls</td><td>(all WLs)</td><td>Number of WLs displayed.</td></tr>
% <tr><td>mode</td><td>'good'</td><td>Display mode. 'good' display channels
% above noise threshold, 'bad' below.</td></tr>
% </table></html>
%
%% Dependencies
% <PlotCapData_help.html PlotCapData> | <istablevar_help.html istablevar>
% 
%% See Also
% <FindGoodMeas_help.html FindGoodMeas> | <PlotCap_help.html PlotCap> |
% <PlotCapMeanLL_help.html PlotCapMeanLL>