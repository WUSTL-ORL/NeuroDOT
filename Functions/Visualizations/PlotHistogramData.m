function h = PlotHistogramData(hist_data, params)

% PLOTHISTOGRAMDATA A basic histogram plotting function.
%
%   PLOTHISTOGRAMDATA(hist_data) plots a histogram of the input data.
%
%   PLOTHISTOGRAMDATA(hist_data, params) allows the user to specify
%   parameters for plot creation.
%
%   h = PLOTHISTOGRAMDATA(...) passes the handles of the plot line objects
%   created.
% 
%   "params" fields that apply to this function (and their defaults):
%       fig_size    [200, 200, 560, 420]    Default figure position vector.
%       fig_handle  (none)                  Specifies a figure to target.
%                                           If empty, spawns a new figure.
%       xlimits     'auto'                  Limits of x-axis.
%       ylimits     'auto'                  Limits of y-axis.
%       bins        200                     Number or array of histogram
%                                           bins.
%
% See Also: PLOTHISTOGRAMSTD, FINDGOODMEAS, PLOTCAPGOODMEAS, HISTOGRAM.
% 
% Copyright (c) 2017 Washington University 
% Created By: Adam T. Eggebrecht
% Eggebrecht et al., 2014, Nature Photonics; Zeff et al., 2007, PNAS.
%
% Washington University hereby grants to you a non-transferable, 
% non-exclusive, royalty-free, non-commercial, research license to use 
% and copy the computer code that is provided here (the Software).  
% You agree to include this license and the above copyright notice in 
% all copies of the Software.  The Software may not be distributed, 
% shared, or transferred to any third party.  This license does not 
% grant any rights or licenses to any other patents, copyrights, or 
% other forms of intellectual property owned or controlled by Washington 
% University.
% 
% YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS 
% PROVIDED AS IS, WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR 
% IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY 
% OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY 
% THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT.  
% IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON 
% UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR 
% CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH 
% THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER 
% IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS 
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

%% Parameters and Initialization.
LineColor = 'w';
BkgdColor = 'k';
FaceAlpha = 0.5;
EdgeAlpha = 0.25;

if ~exist('params', 'var')
    params = [];
end

if ~isfield(params, 'fig_size')  ||  isempty(params.fig_size)
    params.fig_size = [200, 200, 560, 420];
end
if ~isfield(params, 'fig_handle')  ||  isempty(params.fig_handle)
    params.fig_handle = figure('Color', BkgdColor, 'Position', params.fig_size);
    new_fig = 1;
else
    switch params.fig_handle.Type
        case 'figure'
            set(groot, 'CurrentFigure', params.fig_handle);
        case 'axes'
            set(gcf, 'CurrentAxes', params.fig_handle);
    end
end
if ~isfield(params, 'xlimits')  ||  isempty(params.xlimits)
    params.xlimits = 'auto'; % STD (as a %)
end
if ~isfield(params, 'ylimits')  ||  isempty(params.ylimits)
    params.ylimits = 'auto'; % Meas per bin.
end
if ~isfield(params, 'bins')  ||  isempty(params.bins)
    params.bins = 200;
end

%% Plot Data.
hold on

% DOES NOT NEED ND INPUT SUPPORT. "histogram" always uses entire input,
% whether matrix or vector.
h = histogram(hist_data, params.bins, 'FaceAlpha', FaceAlpha,...
    'EdgeColor', LineColor, 'EdgeAlpha', EdgeAlpha);

%% Format axes.
box on

a = gca;
a.Color = BkgdColor;
a.XColor = LineColor;
a.YColor = LineColor;

%% Resize.
ylim(params.ylimits)
xlim(params.xlimits)



%
