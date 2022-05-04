function [mapped,map_out,params] = applycmap(overlay, underlay, params)

% APPLYCMAP Performs color mapping and fuses images with anatomical
% models.
%
%   mapped = APPLYCMAP(overlay) fuses the N-D image array "overlay" with a
%   default flat gray background and applies a number of other default
%   settings to create a scaled and colormapped N-D x 3 RGB image array
%   "mapped".
%
%   mapped = APPLYCMAP(overlay, underlay) fuses the image array "overlay"
%   with the anatomical atlas volume input "underlay" as the background.
%
%   mapped = APPLYCMAP(overlay, underlay, params) allows the user to
%   specify parameters for plot creation.
%
%   "params" fields that apply to this function (and their defaults):
%       TC          0               Direct map integer data values to
%                                   defined color map ("True Color").
%       DR          1000            Dynamic range.
%       Scale       (90% max)       Maximum value to which image is scaled.
%       PD          0               Declares that input image is positive
%                                   definite.
%       Cmap                        Color maps.
%           .P      'jet'           Colormap for positive data values.
%           .N      (none)          Colormap for negative data values.
%           .flipP  0               Logical, flips the positive colormap.
%           .flipN  0               Logical, flips the negative colormap.
%       Th                          Thresholds.
%           .P      (25% max)       Value of min threshold to display
%                                   positive data values.
%           .N      (-Th.P)         Value of max threshold to display
%                                   negative data values.
%       BG          [0.5, 0.5, 0.5] Background color, as an RGB triplet.
%       Saturation  (none)          Field the size of data with values to
%                                   set the coloring saturation. Must be 
%                                   within range [0, 1].
%
% See Also: PLOTSLICES, PLOTINTERPSURFMESH, PLOTLRMESHES, PLOTCAPMEANLL.
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
img_size = size(overlay);
overlay = +overlay(:);                  % Make sure not logical
mapped = zeros(length(overlay), 3);
overlay(~isfinite(overlay))=0;
underlay(~isfinite(underlay))=0;
if ~sum(abs(overlay(overlay~=0)))
    disp(['The Overlay has only elements equal to zero'])
    mapped=mapped+0.5;
    map_out=[0.5,0.5,0.5];
    return
end

if ~exist('underlay', 'var')  ||  isempty(underlay)
    params.underlay = 0;
else
    params.underlay = 1;
    underlay = underlay(:);
    underlay = underlay ./ max(underlay);
end

if ~exist('params', 'var')
    params = [];
end

if ~isfield(params, 'PD')  ||  isempty(params.PD)
    params.PD = 0;
end
if ~isfield(params, 'TC')  ||  isempty(params.TC)
    params.TC = 0;
end
if ~isfield(params, 'DR')  ||  isempty(params.DR)
    params.DR = 1000;
end
if ~isfield(params, 'Scale')  ||  isempty(params.Scale)
    params.Scale = 0.9 * max(overlay);
end
if ~isfield(params, 'Th')  ||  isempty(params.Th)
    params.Th.P = 0.25 * params.Scale;
    params.Th.N = -params.Th.P;
end
if ~isfield(params.Th, 'N')  ||  isempty(params.Th.N)
    params.Th.N = -params.Th.P;
end
if ~isfield(params, 'Cmap')  ||  isempty(params.Cmap) ||...
        (isfield(params.Cmap, 'P')  &&  isempty(params.Cmap.P))
    params.Cmap.P = 'jet';
else
    if ~isstruct(params.Cmap)
        temp = params.Cmap;
        params.Cmap = [];
        params.Cmap.P = temp;
    end
end
if ischar(params.Cmap.P) % Generate positive and negative Cmaps, if present.
    Cmap.P = eval([params.Cmap.P, '(', num2str(params.DR), ');']);
elseif isnumeric(params.Cmap.P)
    Cmap.P = params.Cmap.P;
end
if isfield(params.Cmap, 'flipP')  &&  ~isempty(params.Cmap.flipP)...
        &&  params.Cmap.flipP % Optional colormap flip.
    Cmap.P = flip(Cmap.P, 1);
end
if isfield(params.Cmap, 'N')  &&  ~isempty(params.Cmap.N)
    if ischar(params.Cmap.N)
        Cmap.N = eval([params.Cmap.N, '(', num2str(params.DR), ');']);
    elseif isnumeric(params.Cmap.N)
        Cmap.N = params.Cmap.N;
    end
    if isfield(params.Cmap, 'flipN')  &&  ~isempty(params.Cmap.flipN)  &&...
            params.Cmap.flipN
        Cmap.N = flip(Cmap.N, 1);
    end
    params.PD = 1;
end
if ~isfield(params, 'BG')  ||  isempty(params.BG)
    params.BG = [0.5, 0.5, 0.5];
end

%% Get Background.
bg = find(overlay == 0);

%% Truecolor
if params.TC
    % Color in RGB channels.
    Nc = size(Cmap.P, 1);
    for j = 1:Nc
        mapped(overlay == j, 1) = Cmap.P(j, 1);
        mapped(overlay == j, 2) = Cmap.P(j, 2);
        mapped(overlay == j, 3) = Cmap.P(j, 3);
    end
    
    % Color in background.
    if ~params.underlay
        mapped(bg, :) = repmat(params.BG, length(bg), 1);
    else
        mapped(bg, :) = [underlay(bg) .* params.BG(1), ...
            underlay(bg) .* params.BG(2), underlay(bg) .* params.BG(3)];
    end
else
    %% Label Data Outside Thresholds to Background.
    bg = union(bg, intersect(find(overlay <= params.Th.P), find(overlay >= params.Th.N)));
    
    %% Scaling and Dynamic Range.
    % All below threshold set to gray, all above set to colormap.
    if params.PD % Treat as pos def data?
        overlay = (overlay ./ params.Scale) .* (params.DR); % Normalize and scale to colormap
        fgP = find(overlay > 0); % Populate Pos to color
        fgN = find(overlay < 0); % Populate Neg to color
    else
        overlay = (overlay ./ params.Scale) .* (params.DR / 2); % Normalize and scale to colormap
        overlay = overlay + (params.DR / 2); % Shift data
        fgP = find(overlay ~= (params.DR / 2)); % Populate Pos to color
        overlay(overlay <= 0) = 1; % Correct for neg clip
    end
    
    overlay(overlay >= params.DR) = params.DR; % Correct for pos clip
    
    %% Color in RGB channels.
    mapped(fgP, 1) = Cmap.P(ceil(overlay(fgP)), 1);
    mapped(fgP, 2) = Cmap.P(ceil(overlay(fgP)), 2);
    mapped(fgP, 3) = Cmap.P(ceil(overlay(fgP)), 3);
    
    % Pos def treatment.
    if params.PD
        if isfield(params.Cmap, 'N')
            overlay(-overlay >= params.DR) = -params.DR; % Correct for clip.
            mapped(fgN, 1) = Cmap.N(ceil(-overlay(fgN)), 1);
            mapped(fgN, 2) = Cmap.N(ceil(-overlay(fgN)), 2);
            mapped(fgN, 3) = Cmap.N(ceil(-overlay(fgN)), 3);
        else
            bg = union(bg, fgN); % If only pos, neg values to background.
        end
    end
    
    %% Apply background coloring.
    if ~params.underlay
        mapped(bg, :) = repmat(params.BG, length(bg), 1);
    else
        mapped(bg, :) = [underlay(bg) .* params.BG(1),...
            underlay(bg) .* params.BG(2), underlay(bg) .* params.BG(3)];
    end
end

%% Apply Saturation if available
if isfield(params,'Saturation')
    nbg = setdiff(1:size(mapped, 1), bg);
     if ~params.underlay
         mapped(nbg, :) = bsxfun(@times, params.Saturation(nbg), mapped(nbg, :))...
             + bsxfun(@times, params.BG, 1 - params.Saturation(nbg));
     else
         mapped(nbg, :) = bsxfun(@times, params.Saturation(nbg), mapped(nbg, :))...
             + bsxfun(@times, params.BG, (1 - params.Saturation(nbg)) .* underlay(nbg, :));
     end
end

%% Reshape to original size + RGB.
if (img_size(2) == 1)  &&  (length(img_size) == 2)
    mapped = reshape(mapped, [img_size(1), 3]);
else
    mapped = reshape(mapped, [img_size, 3]);
end

%% Create a final color map.
if ~params.underlay
    thresh_zone_color = params.BG; % Same color as background if no underlay.
else
    thresh_zone_color = [0, 0, 0]; % Black if there is an underlay.
end

map_out = Cmap.P;
if params.TC
    map_out = params.Cmap.P;
else
    if params.PD
        thP = floor(params.Th.P / params.Scale * params.DR);
        map_out(1:thP, :) = repmat(thresh_zone_color, thP, 1);
        if isfield(Cmap, 'N')
            thN = floor(params.Th.N / params.Scale * params.DR);
            tempN = Cmap.N;
            tempN(1:thN, :) = repmat(thresh_zone_color, thN, 1);
            map_out = [flip(tempN, 1); map_out];
        end
    else
        thP = floor(params.Th.P / params.Scale * params.DR / 2 + (params.DR / 2));
        thN = ceil(params.Th.N / params.Scale * params.DR / 2 + (params.DR / 2));
        if thN<1,thN=1;end
        map_out(thN:thP, :) = repmat(thresh_zone_color, thP-thN+1, 1);
    end
end



%
