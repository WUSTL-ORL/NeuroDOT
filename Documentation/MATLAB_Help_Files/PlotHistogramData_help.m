%% PlotHistogramData
% A basic histogram plotting function.
%
%% Description
% |PlotHistogramData(hist_data)| plots a histogram of the input data.
%
% |PlotHistogramData(hist_data, params)| allows the user to specify
% parameters for plot creation.
%
% |h = PlotHistogramData(...)| passes the handles of the plot line objects
% created.
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
% <tr><td>xlimits</td><td>'auto'</td><td>Limits of x-axis.</td></tr>
% <tr><td>ylimits</td><td>'auto'</td><td>Limits of y-axis.</td></tr>
% <tr><td>bins</td><td>200</td><td>Number or array of histogram
% bins.</td></tr>
% </table></html>
%
%% See Also
% <PlotHistogramSTD_help.html PlotHistogramSTD> | <FindGoodMeas_help.html
% FindGoodMeas> | <PlotCapGoodMeas_help.html PlotCapGoodMeas> |
% <matlab:doc('histogram') histogram>