function [vox,dim]=getvox(nodes,G,flags)

% Find voxel grid and limits.  Voxel grid is defined as a subset of the
% volume voxellated space.
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

%% Parameters and Initialization
if ~isfield(flags,'gthresh'),flags.gthresh=1e-5;end
if ~isfield(flags,'voxmm'),flags.voxmm=1;end
if ~isfield(flags.info,'nVx'),flags.info.nVx=102;end
if ~isfield(flags.info,'nVy'),flags.info.nVy=145;end
if ~isfield(flags.info,'nVz'),flags.info.nVz=58;end
if ~isfield(flags.info,'acq'),flags.info.acq='sagittal';end
if ~isfield(flags.info,'center'),flags.info.center=[50,-73,-80];end
if ~isfield(flags.info,'mmppix'),flags.info.mmppix=[1,-1,-1];end
if ~isfield(flags, 'center_shift'), flags.center_shift = [-1,1,1];end

%% Determine Threshold and nodes within the threshold range
nodelevel=max(squeeze(sum(abs(G),1)),[],1);
m=max(nodelevel);
nkeep=find(nodelevel(:)>=flags.gthresh*m);

xkeep=nodes(nkeep,1);
xmin=floor(min(xkeep));
xmax=ceil(max(xkeep));
while mod(xmax-xmin,flags.voxmm)
    xmax=xmax+1;
end
xmax=xmax-1;

ykeep=nodes(nkeep,2);
ymin=floor(min(ykeep));
ymax=ceil(max(ykeep));
while mod(ymax-ymin,flags.voxmm)
    ymax=ymax+1;
end
ymax=ymax-1;

zkeep=nodes(nkeep,3);
zmin=floor(min(zkeep));
zmax=ceil(max(zkeep));
while mod(zmax-zmin,flags.voxmm)
    zmax=zmax+1;
end
zmax=zmax-1;


%% Construct dim structure
dim=struct('xmin',xmin,'xmax',xmax,...
    'ymin',ymin,'ymax',ymax,'zmin',zmin,'zmax',zmax);

switch sign(min(nodes(:)))
    case -1 % coordinate based mesh nodes
% Have to flip x bc coord (mesh) vs indexing (vox) conventions.
dim.xv=[dim.xmax:(-flags.voxmm):dim.xmin]+0;
dim.yv=[dim.ymin:flags.voxmm:dim.ymax]+0;
dim.zv=[dim.zmin:flags.voxmm:dim.zmax]+0;

MPR111 = change_space_coords(...
    [flags.info.nVx,flags.info.nVy,flags.info.nVz], ...
    flags.info,'coord');
dim111=[dim.xv(end),dim.yv(end),dim.zv(end)];
dr=(MPR111-dim111);
dim.center=flags.info.center+dr-flags.center_shift.*flags.voxmm; % fixed 240118


    case 1 % index based mesh nodes
dim.xv=dim.xmin:flags.voxmm:dim.xmax;
dim.yv=dim.ymin:flags.voxmm:dim.ymax;
dim.zv=dim.zmin:flags.voxmm:dim.zmax;

MPRsize=[flags.info.nVx,flags.info.nVy,flags.info.nVz];
centerMPR=flags.info.center;
dr=[MPRsize(1)-(dim.xmax+1),MPRsize(2)-(dim.ymax+1),...
    MPRsize(3)-(dim.zmax+1)];
dim.center=centerMPR-(dr.*flags.info.mmppix);

end

dim.nVx=numel(dim.xv);
dim.nVy=numel(dim.yv);
dim.nVz=numel(dim.zv);
dim.nVt=dim.nVx*dim.nVy*dim.nVz;

dim.sV=flags.voxmm;
dim.acq=flags.info.acq;
dim.orientation=flags.info.acq;
dim.mmppix=[flags.voxmm,-flags.voxmm,-flags.voxmm];



%% Construct vox
vox=zeros(dim.nVx,dim.nVy,dim.nVz,3);
for x=1:dim.nVx
    for y=1:dim.nVy
        for z=1:dim.nVz
            vox(x,y,z,:)=[dim.xv(x),dim.yv(y),dim.zv(z)];
        end
    end
end