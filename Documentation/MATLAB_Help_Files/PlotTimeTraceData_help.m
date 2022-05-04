%% PlotTimeTraceData
% A basic time traces plotting function.
%
%% Description
% |PlotTimeTraceData(data, info)| takes a light-level array |data| of the
% MEAS x TIME format, and plots its time traces.
% 
% |PlotTimeTraceData(data, info, params)| allows the user to specify
% parameters for plot creation.
% 
% |h = PlotTimeTraceData(...)| passes the handles of the plot line objects
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
% <tr><td>xscale</td><td>'linear'</td><td>Scaling of x-axis.</td></tr>
% <tr><td>ylimits</td><td>'auto'</td><td>Limits of y-axis.</td></tr>
% <tr><td>yscale</td><td>'linear'</td><td>Scaling of y-axis.</td></tr>
% </table></html>
%
%% See Also
% <PlotTimeTraceAllMeas_help.html PlotTimeTraceAllMeas> |
% <PlotTimeTraceMean_help.html PlotTimeTraceMean>