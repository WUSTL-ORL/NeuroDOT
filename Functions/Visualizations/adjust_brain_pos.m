function [Lnodes, Rnodes] = adjust_brain_pos(meshL, meshR, params)

% ADJUST_BRAIN_POS Repositions mesh orientations for display.
%
%   [Lnodes, Rnodes] = ADJUST_BRAIN_POS(meshL, meshR) takes the left and
%   right hemispheric meshes "meshL" and "meshR", respectively, and
%   repositions them to the proper perspective for display.
%
%   [Lnodes, Rnodes] = ADJUST_BRAIN_POS(meshL, meshR, params) allows the
%   user to specify parameters for plot creation.
%
%   "params" fields that apply to this function (and their defaults):
%       ctx         'std'       Defines inflation of mesh.
%       orientation 't'         Select orientation of volume. 't' for
%                               transverse, 's' for sagittal.
%       view        'lat'       Sets the view perspective.
% 
% Dependencies: ROTATE_CAP, ROTATION_MATRIX.
% 
% See Also: PLOTLRMESHES.
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

%% Parameters and Initialization
dy = -5;

if ~exist('params', 'var')  ||  isempty(params)
    params = [];
end

if ~isfield(params, 'ctx')  ||  isempty(params.ctx)
    params.ctx = 'std';
end
if ~isfield(params, 'orientation')  ||  isempty(params.orientation)
    params.orientation = 't';
end
if ~isfield(params, 'view')  ||  isempty(params.view)
    params.view = 'lat';
end

%% Choose inflation.
switch params.ctx
    case 'std'
        str = 'nodes';
    case 'inf'
        str = 'Inodes';
    case 'vinf'
        str = 'VInodes';
    case 'flat'
        str = 'Fnodes';
end
Lnodes = meshL.(str);
Rnodes = meshR.(str);

%% Standardize orientation.
% If volume is not in transverse orientation, put it there
if strcmp(params.ctx, 'std')
    switch params.orientation 
        case 's'
            Nln = size(Lnodes, 1);
            temp = cat(1, Lnodes, Rnodes);
            temp = rotate_cap(temp, [-90, 0, 90]);
            Lnodes = temp(1:Nln, :);
            Rnodes = temp((Nln + 1):end, :);
        case 'c'
            % Not yet supported.
        case 't'
            % Not yet supported.
    end
end


%% Normalize L and R nodes to their max and min, respectively.
switch params.ctx
    case 'std'
        str = 'nodes';
    case 'inf'
        str = 'Inodes';
Lnodes(:, 1) = Lnodes(:, 1) - max(Lnodes(:, 1));
Rnodes(:, 1) = Rnodes(:, 1) - min(Rnodes(:, 1));

    case 'vinf'
        str = 'VInodes';
Lnodes(:, 1) = Lnodes(:, 1) - max(Lnodes(:, 1));
Rnodes(:, 1) = Rnodes(:, 1) - min(Rnodes(:, 1));
    case 'flat'
        str = 'Fnodes';
        Lnodes(:, 1) = Lnodes(:, 1) - max(Lnodes(:, 1));
        Rnodes(:, 1) = Rnodes(:, 1) - min(Rnodes(:, 1));
        Nln = size(Lnodes, 1);
        temp = cat(1, Lnodes, Rnodes);
        temp1 = rotate_cap(temp, [-90, 0, 90]);
        temp2 = rotate_cap(temp, [-90,-0,-90]);%[-90,-10,-90]);
        Lnodes = temp1(1:Nln, :);
        Rnodes = temp2((Nln + 1):end, :);
%         Rnodes = Rnodes(:, 1,1) +100;
        
      
end


%% Rotate if necessary
if strcmp(params.view, 'lat')  ||  strcmp(params.view, 'med')

cmL = mean(Lnodes, 1); % rotate right hemi  for visualization
    cmR = mean(Rnodes, 1);
    rm = rotation_matrix('z', pi);
    
    switch params.view % Rotate
        case 'lat'
            Rnodes = (Rnodes - (repmat(cmR, size(Rnodes, 1), 1))) *...
                rm + (repmat(cmR, size(Rnodes, 1), 1));
        case 'med'
            Lnodes = (Lnodes - (repmat(cmL, size(Lnodes, 1), 1))) *...
                rm + (repmat(cmL, size(Lnodes, 1), 1));
    end
    Rnodes(:, 1) = Rnodes(:, 1) + (cmL(:, 1) - cmR(:, 1)); % Shift over to same YZ plane
    Rnodes(:, 2) = Rnodes(:, 2) - max(Rnodes(:, 2)) + min(Lnodes(:, 2)) + dy;
end


%% Translate if necessary
if strcmp(params.ctx, 'flat')
    Lnodes(:,2) = Lnodes(:,2) - 30;
    Rnodes(:,2) = Rnodes(:,2) - 40;
end


%
