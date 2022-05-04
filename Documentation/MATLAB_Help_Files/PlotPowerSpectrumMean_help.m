%% PlotPowerSpectrumMean
% A mean frequency-domain time trace visualization.
%
%% Description
% |PlotPowerSpectrumMean(data, info)| takes a light-level array
% |data| of the MEAS x TIME format, and creates a plot of the absolute
% value squared of its FFT. A mean of this is taken for specified
% groupings.
% 
% |PlotPowerSpectrumMean(data, info, params)| allows the user to
% specify parameters for plot creation.
% 
% |PlotPowerSpectrumMean(data, info, params, framerate)| allows the user
% to specify the data framerate.
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
% <tr><td>useGM</td><td>1</td><td>Use Good Measurements.</td></tr>
% <tr><td>xlimits</td><td>[1e-3, 5]</td><td>Limits of x-axis.</td></tr>
% <tr><td>xscale</td><td>'log'</td><td>Scaling of x-axis.</td></tr>
% <tr><td>ylimits</td><td>[0, 3e-5]</td><td>Limits of y-axis.</td></tr>
% <tr><td>yscale</td><td>'linear'</td><td>Scaling of y-axis.</td></tr>
% </table></html>
%
%% Dependencies
% <PlotPowerSpectrumData_help.html PlotPowerSpectrumData> |
% <istablevar_help.html istablevar>
%
%% See Also
% <PlotPowerSpectrumAllMeas_help.html PlotPowerSpectrumAllMeas>