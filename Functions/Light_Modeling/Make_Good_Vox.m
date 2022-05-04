function [Gs,Gd,dc,dim]=Make_Good_Vox(Gs,Gd,dc,dim,mesh,flags)
%
% To keep filesizes under some control, remove voxels with below threshold
% sensitivity. dim.Good_Vox provides efficient bookkeeping.
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


%% Find threshold
Vthresh=flags.gthresh*0.1;
[Ncol,Nd,Nvox]=size(Gd);
Ns=size(Gs,2);
Gd=reshape(Gd,Ncol*Nd,Nvox);
Gs=reshape(Gs,Ncol*Ns,Nvox);

voxlevelT=max(squeeze(sum(abs(Gd),1))+squeeze(sum(abs(Gs),1)),[],1); % Max vox Gs&Gd


%% Find voxels below threshold
Kill=find(voxlevelT<Vthresh*max(voxlevelT));

dim.Good_Vox=setdiff([1:Nvox],Kill)';
Gd=reshape(Gd,Ncol,Nd,Nvox);
Gs=reshape(Gs,Ncol,Ns,Nvox);
Gd=Gd(:,:,dim.Good_Vox);
Gs=Gs(:,:,dim.Good_Vox);
dc=dc(:,dim.Good_Vox);

clear KillD KillS Kill

%% Include grid coordinates from mesh
dim.Ngv=length(dim.Good_Vox);
dim.tpos=mesh.source.coord;                         % mesh xyz indices
