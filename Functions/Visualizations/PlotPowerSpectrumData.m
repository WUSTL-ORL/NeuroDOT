function h = PlotPowerSpectrumData(data, info, params, framerate)

% PLOTPOWERSPECTRUMDATA A basic frequency-domain time traces plotting
% function. 
%
%   PLOTPOWERSPECTRUMDATA(data, info) takes a light-level array "data" of
%   the MEAS x TIME format, and creates a plot of the absolute value
%   squared of the FFT.
%
%   PLOTPOWERSPECTRUMDATA(data, info, params) allows the user to specify
%   parameters for plot creation.
% 
%   PLOTPOWERSPECTRUMDATA(data, info, params, framerate) allows the user to
%   specify the data framerate.
% 
%   h = PLOTPOWERSPECTRUMDATA(...) passes the handles of the plot line
%   objects created.
%
%   "params" fields that apply to this function (and their defaults):
%       fig_size    [200, 200, 560, 420]    Default figure position vector.
%       fig_handle  (none)                  Specifies a figure to target.
%                                           If empty, spawns a new figure.
%       xlimits     'auto'                  Limits of x-axis.
%       xscale      'log'                   Scaling of x-axis.
%       ylimits     'auto'                  Limits of y-axis.
%       yscale      'linear'                Scaling of y-axis.
%
% Dependencies: FFT_TTS.
%
% See Also: PLOTPOWERSPECTRUMALLMEAS, PLOTPOWERSPECTRUMMEAN.
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

dims = size(data);
Nt = dims(end);
NDtf = (ndims(data) > 2);

if ~exist('framerate', 'var')  ||  isempty(framerate)
    if (isfield(info, 'system')  &&  ~isempty(info.system))...
            &&  (isfield(info.system, 'framerate')  &&...
            ~isempty(info.system.framerate))
        framerate = info.system.framerate;
    end
end

if ~exist('params', 'var')
    params = [];
end

if ~isfield(params, 'fig_size')  ||  isempty(params.fig_size)
    params.fig_size = [200, 200, 560, 420];
end
if ~isfield(params, 'fig_handle')  ||  isempty(params.fig_handle)
    params.fig_handle = figure('Color', BkgdColor,...
        'Position', params.fig_size);
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
    params.xlimits = 'auto';
end
if ~isfield(params, 'xscale')  ||  isempty(params.xscale)
    params.xscale = 'log';
end
if ~isfield(params, 'ylimits')  ||  isempty(params.ylimits)...
        ||  (all(all(data > params.ylimits(2)))...
        &&  all(all(data < params.ylimits(1))))
    params.ylimits = 'auto';
end
if ~isfield(params, 'yscale')  ||  isempty(params.yscale)
    params.yscale = 'linear';
end

%% N-D Input.
if NDtf
    data = reshape(data, [], Nt);
end

%% Calculate & Plot Data.
[~, ftdomain, ftpower] = fft_tts(data, framerate);

h = plot(ftdomain, ftpower, 'LineWidth', 1);

%% Format Plot.
box on

a = gca;
a.Color = BkgdColor;
a.XColor = LineColor;
a.YColor = LineColor;

%% Scaling.
ylim(params.ylimits)
xlim(params.xlimits)
a.YScale = params.yscale;
a.XScale = params.xscale;


%% Labels.
% Since this is always going to be power spectrum, we label here.
ylabel('Power (A.U.)', 'Color', LineColor)
xlabel('Frequency (Hz)', 'Color' ,LineColor)



%
