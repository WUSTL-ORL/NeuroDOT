%% PlotPowerSpectrumData
% A basic frequency-domain time traces plotting function. 
%
%% Description
% |PlotPowerSpectrumData(data, info, framerate)| takes a light-level array
% |data| of the MEAS x TIME format, and creates a plot of the absolute
% value squared of the FFT.
%
% |PlotPowerSpectrumData(data, info, framerate, params)| allows the user to
% specify parameters for plot creation.
% 
% |h = PlotPowerSpectrumData(...)| passes the handles of the plot line
% objects created.
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
% <tr><td>xscale</td><td>'log'</td><td>Scaling of x-axis.</td></tr>
% <tr><td>ylimits</td><td>'auto'</td><td>Limits of y-axis.</td></tr>
% <tr><td>yscale</td><td>'linear'</td><td>Scaling of y-axis.</td></tr>
% </table></html>
%
%% Dependencies
% <fft_tts_help.html fft_tts>
%
%% See Also
% <PlotPowerSpectrumAllMeas_help.html PlotPowerSpectrumAllMeas> |
% <PlotPowerSpectrumMean_help.html PlotPowerSpectrumMean>