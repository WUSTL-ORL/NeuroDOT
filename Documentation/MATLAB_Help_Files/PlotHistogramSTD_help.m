%% PlotHistogramSTD
% A histogram of channel standard deviations.
%
%% Description
% |PlotHistogramSTD(info)| plots a histogram of the distribution of channel
% variances for specified groupings. A red bar will show the threshold
% value.
%
% |PlotHistogramSTD(info, params)| allows the user to specify parameters
% for plot creation.
%
% |PlotHistogramSTD(info, params, bthresh)| allows the user to specify a
% threshold value. If no input is supplied, the default value of |0.075|
% will be used.
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
% <tr><td>dimension</td><td>'2D'</td><td>Specifies either a 2D or 3D plot
% rendering.</td></tr>
% <tr><td>rlimits</td><td>(all R2D)</td><td>Limits of pair radii
% displayed.</td></tr>
% <tr><td>Nnns</td><td>(all NNs)</td><td>Number of NNs displayed.</td></tr>
% <tr><td>Nwls</td><td>(all WLs)</td><td>Number of WLs displayed.</td></tr>
% <tr><td>xlimits</td><td>[0, 50]</td><td>Limits of x-axis.</td></tr>
% <tr><td>ylimits</td><td>[0, 200]</td><td>Limits of y-axis.</td></tr>
% <tr><td>bins</td><td>[0:0.5:100]</td><td>Number or array of histogram
% bins.</td></tr>
% </table></html>
%
%% Dependencies
% <PlotHistogramData_help.html PlotHistogramData> | <istablevar_help.html
% istablevar>
% 
%% See Also
% <FindGoodMeas_help.html FindGoodMeas> | <PlotCapGoodMeas_help.html
% PlotCapGoodMeas> | <matlab:doc('histogram') histogram>