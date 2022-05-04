function Mout=voxel(M,t,p,elements)

% Voxelate data (M) based on a node space (described by t,p,elements) into
% a space defined by dim.
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

%% Initialize and prepare
dims = size(M);
M2=reshape(M,[],dims(end));
Nm=size(M2,1);
Nvox=size(t,1);
Mout=zeros(Nm,Nvox); % Noptodes*Ncolors, Nvox

%% Voxellate only for voxels that actually contain nodes of the mesh
keep=~isnan(t);
Nkeep=sum(keep);
elementList=elements(t(keep),:);

% Interpolate
for n=1:Nm                                % For each measurement
    Jam = zeros(Nkeep,4);
    for j = 1:Nkeep                         % For each good voxel
        Jam(j,:) = M2(n,elementList(j,:));
    end
    Mout(n,keep) = sum(p(keep,:).*Jam,2);
end

Mout=reshape(Mout,[dims(1:(end-1)),Nvox]);