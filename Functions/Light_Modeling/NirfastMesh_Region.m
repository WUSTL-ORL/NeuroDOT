function mesh=NirfastMesh_Region(mask,meshname,param)

% This function passes the nirfast meshing code a volume segmented into a
% single region.  Then, we label the nodes of the mesh after it is created.
%  This is done to keep nirfast from meddling with region sizes etc.
%
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
% mask
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


%% Generate homogeneous mesh
Scalp=max(mask(:));
if ~exist('param','var'),param='hd';end
if ~isfield(param,'Mode'),param.Mode=0;end
mask_old=mask;
mask(mask~=0)=Scalp;
mesh=NirfastMesh(mask,meshname,param,param.Mode);


%% Label nodes based on region overlay with segmented volume
if ~param.Mode

if ~isfield(param,'Offset'),param.Offset=[0,0,0];end
if (isfield(param,'info') && norm(param.Offset-[0,0,0]))
    % if mesh is made in coord-mapped space
    Vnodes=change_space_coords(mesh.nodes,param.info,'idx'); 
else
    Vnodes=mesh.nodes;
end

% Interpolate region labels using nearest neighbor method
F = griddedInterpolant(mask_old,'nearest','nearest');
mesh.region = F(Vnodes);

% if node outside of mask, set equal to boundary region.
if isstruct(param)
    if isfield(param,'r0')
mesh.region(mesh.region==0)=param.r0; 
    else
mesh.region(mesh.region==0)=Scalp; 
    end
else
mesh.region(mesh.region==0)=Scalp; 
end

disp('<<<Saving mesh')
save_mesh(mesh,meshname);
end