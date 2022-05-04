function PlotCap(info, params)

% PLOTCAP A visualization of the cap grid.
%
%   PLOTCAP(info) plots a basic diagram of the cap based on the cap
%   metadata stored in "info.optodes". All optodes are numbered; sources
%   are colored light red, and detectors light blue.
%
%   PLOTCAP(info, params) allows the user to specify parameters for
%   plot creation.
%
%   "params" fields that apply to this function (and their defaults):
%       fig_size    [20, 200, 1240, 420]    Default figure position vector.
%       fig_handle  (none)                  Specifies a figure to target.
%                                           If empty, spawns a new figure.
%       dimension   '2D'                    Specifies either a 2D or 3D
%                                           plot rendering.
%
% Dependencies: PLOTCAPDATA.
% 
% See Also: PLOTCAPGOODMEAS, PLOTCAPMEANLL.
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
% Ns = length(unique(info.pairs.Src));
% Nd = length(unique(info.pairs.Det));
Ns = length(info.optodes.spos2);
Nd = length(info.optodes.dpos2);
LineColor = 'w';
BkgdColor = 'k';
if ~exist('params','var'),params=struct;end
if ~isfield(params,'SrcColor'),params.SrcColor=[1, 0.75, 0.75];end
if ~isfield(params,'DetColor'),params.DetColor=[0.55, 0.55, 1];end

if ~isfield(params,'rhombus'),params.rhombus=1;end

params.mode = 'textpatch';

if ~isfield(params, 'dimension')  ||  isempty(params.dimension)
    params.dimension = '2D'; % '2D' | '3D'
end
if ~isfield(params, 'fig_handle')  ||  isempty(params.fig_handle)
    if ~isfield(params, 'fig_size')  ||  isempty(params.fig_size)
        switch params.dimension
            case '2D'
                params.fig_size = [20, 200, 1240, 420];
            case '3D'
                params.fig_size = [20, 200, 560, 560];
        end
    end
    params.fig_handle = figure('Color', BkgdColor, 'Position', params.fig_size);
else
    switch params.fig_handle.Type
        case 'figure'
            set(groot, 'CurrentFigure', params.fig_handle);
        case 'axes'
            set(gcf, 'CurrentAxes', params.fig_handle);
    end
end

SrcRGB = repmat(params.SrcColor, Ns, 1);
DetRGB = repmat(params.DetColor, Nd, 1);

%% Send to PlotCapData.
PlotCapData(SrcRGB, DetRGB, info, params)


%
