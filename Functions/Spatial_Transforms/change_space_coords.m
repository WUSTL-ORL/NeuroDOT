function coord_out = change_space_coords(coord_in, space_info, output_type)

% CHANGE_SPACE_COORDS Applies a look up to change 3D coordinates into a new
% space.
% 
%   coord_out = CHANGE_SPACE_COORDS(coord_in, space_info, output_type) takes
%   a set of coordinates "coord_in" of the initial space "output_type", and
%   converts them into the new space defined by the structure "space_info",
%   which is then output as "coord_out".
% 
% See Also: AFFINE3D_IMG.
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

%% Parameters and initialization
if ~exist('output_type', 'var')
    output_type = 'coord';
end

% Define the voxel space.
nVxA = space_info.nVx;
nVyA = space_info.nVy;
nVzA = space_info.nVz;
drA = space_info.mmppix;
centerA = space_info.center; % i.e., center = coordinate of center of voxel with index [-1,-1,-1]
nV = size(coord_in, 1);

% Preallocate.
coord_out = zeros(size(coord_in));

%% Create coordinates for each voxel index.
X = (drA(1).*[nVxA:-1:1]-centerA(1))';
Y = (drA(2).*[nVyA:-1:1]-centerA(2))';
Z = (drA(3).*[nVzA:-1:1]-centerA(3))';

%% Convert coordinates to new space.
switch output_type
    case 'coord' % ATLAS/4DFP/ETC COORDINATE SPACE
        for j = 1:nV
            x = coord_in(j, 1);
            if ((floor(x) > 0)  &&  (floor(x) <= nVxA))
                coord_out(j, 1) = X(floor(x)) - drA(1) * (x - floor(x));
            elseif floor(x) < 1
                coord_out(j, 1) = X(1) - drA(1) * (x - 1);
            elseif floor(x) > nVxA
                coord_out(j, 1) = X(nVxA) - drA(1) * (x - nVxA);
            end
            
            y = coord_in(j, 2);
            if ((floor(y) > 0)  &&  (floor(y) <= nVyA))
                coord_out(j, 2) = Y(floor(y)) - drA(2) * (y - floor(y));
            elseif floor(y) < 1
                coord_out(j, 2) = Y(1) - drA(2) * (y - 1);
            elseif floor(y) > nVyA
                coord_out(j, 2) = Y(nVyA) - drA(2) * (y - nVyA);
            end
            
            z = coord_in(j, 3);
            if ((floor(z) > 0)  &&  (floor(z) <= nVzA))
                coord_out(j, 3) = Z(floor(z)) - drA(3) * (z - floor(z));
            elseif floor(z) < 1
                coord_out(j, 3) = Z(1) - drA(3) * (z - 1);
            elseif floor(z) > nVzA
                coord_out(j, 3) = Z(nVzA) - drA(3) * (z - nVzA);
            end
        end
    case 'idx' % MATLAB INDEX SPACE
        for j = 1:nV
            [~, coord_out(j, 1)] = min(abs(coord_in(j, 1) - X));
            [~, coord_out(j, 2)] = min(abs(coord_in(j, 2) - Y));
            [~, coord_out(j, 3)] = min(abs(coord_in(j, 3) - Z));
        end
    case 'idxC' % MATLAB INDEX SPACE WITH NO ROUNDING
        for j = 1:nV
            [~, foo] = min(abs(coord_in(j, 1) - X));
            coord_out(j, 1) = foo + (X(foo) - coord_in(j, 1)) / drA(1);
            [~, foo] = min(abs(coord_in(j, 2) - Y));
            coord_out(j, 2) = foo + (Y(foo) - coord_in(j, 2)) / drA(2);
            [~, foo] = min(abs(coord_in(j, 3) - Z));
            coord_out(j, 3) = foo + (Z(foo) - coord_in(j, 3)) / drA(3);
        end
end



%
