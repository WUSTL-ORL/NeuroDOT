function PlotSlicesMov(underlay, infoVol, params, overlay, fn, infoFR)

% PLOTSLICESMOV Creates an interactive three-slice plot.
%
%   PLOTSLICESMOV(underlay, infoVol, params, overlay, fn) uses the same
%   basic inputs as the function PLOTSLICES to create a video of a 4D
%   volume, saved under the filename "fn".
% 
%   PLOTSLICESMOV(..., infoFR) allows the user to input an "info" structure
%   or a framerate for the video to be recorded with. The default framerate
%   is 1 Hz.
%
%   "params" fields that apply to this function (and their defaults):
%       fig_size    [20, 200, 1240, 420]    Default figure position vector.
%       fig_handle  (none)                  Specifies a figure to target.
%
%   Note: APPLYCMAP has further options for using "params" to specify
%   parameters for the fusion, scaling, and colormapping process.
%
% Dependencies: APPLYCMAP, PLOTSLICES.
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

%% Parameters and initialization.
[nVxU, nVyU, nVzU, nVtU] = size(underlay);
[nVxO, nVyO, nVzO, nVtO] = size(overlay);
nVt = [];
LineColor = 'w';
BkgdColor = 'k';

if ~exist('overlay', 'var')  ||  isempty(overlay)
    overlay = [];
    nVt = nVtU;
else
    if (nVxO ~= nVxU)  ||  (nVyO ~= nVyU)  ||  (nVzO ~= nVzU)
        error('*** Error: "underlay" spatial size does not match "overlay". ***')
    end
    if (nVtO ~= nVtU)  &&  (nVtU > 1)
        error('*** Error: "underlay" time length does not match "overlay". ***')
    end
    nVt = nVtO;
end

if exist('infoFR', 'var')
    if isnumeric(infoFR)
        FR = infoFR;
    elseif isstruct(infoFR)... % Check if it's ND2's "info" structure.
            &&  isfield(infoFR, 'system')  &&  ~isempty(infoFR.system)
        if isfield(infoFR.system, 'framerate')  &&  ~isempty(infoFR.system.framerate)
            FR = infoFR.system.framerate;
        elseif isfield(infoFR.system, 'init_framerate')  &&  ~isempty(infoFR.system.init_framerate)
            FR = infoFR.system.init_framerate;
        else
            error('*** Error: No framerate recognized in "info" structure provided. ***')
        end
    elseif ischar(infoFR)  &&  ~isnan(str2double(infoFR))
        FR = str2double(infoFR);
    else
        disp(['*** Warning: No framerate recognized in ''infoFR''',...
            ' input. Using default framerate. ***'])
        FR = 1;
    end
else
    FR = 1;
end

if ~isfield(params, 'fig_size')  ||  isempty(params.fig_size)
    params.fig_size = [20, 200, 1240, 420];
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

%% Initialize VideoWriter.
v = VideoWriter(fn);
v.FrameRate = FR; % 1 Hz per second.
open(v)

%% Create Video.
for k = 1:nVt
    PlotSlices(underlay, infoVol, params, overlay(:, :, :, k))
    annotation(gcf, 'textbox', [.1, .1, .3, .2],...
        'String', ['T = ' num2str(k) 's'], 'FitBoxToText', 'on',...
        'Color', LineColor, 'EdgeColor', LineColor, 'BackgroundColor', BkgdColor)
    frame = getframe(gcf);
    writeVideo(v, frame)
    pause(0.1)
end
close(v)




%
