function Draw_ROIs_on_Cortex(Anat,ROI,params)
%
% This function will draw ROI spheres on a model of the cortical surface.
% Anat is a structure that describes the anatomy.
% ROI is a structure that describes the regions of interest. 
%
% Anat fields:
%   CtxL    mesh for left cortex. mesh includes 2 fields, nodes and
%           elements.
%   CtxR    mesh for right cortex
%   view    perspective of image {'lat'(default),'med','post','dorsal',)}
%          (optional)
%   ctx    cortical surface type. {'std'(default),'inf','vinf'} (optional)
% ROI fields:
%   coord   [x,y,z] coordinates for the ROIs. It is assumed that the ROIs
%           are in the same coordinate space as the anatomy.
%   radius  size of spheres to be drawn (default = 5). (optional)
%   color   color of each sphere. Typically this matches some Network
%           organization.  If this is also not present, all
%           spheres will be colored white. (optional)
%   Network Nx2 array of ROI indices (1st col) and network membership (2nd
%           column). (optional)

% Copyright (c) 2017 Washington University 
% Authors: Adam T. Eggebrecht, Zachary E. Markow
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

%% Under the hood Parameters
if ~exist('params','var'),params=struct;end
if ~isfield(params,'OL'),params.OL=1;end
if ~isfield(params,'PD'),params.PD=1;end
if ~isfield(params,'TC'),params.TC=1;end
if ~isfield(params,'Th'),params.Th.P=0.001;params.Th.N=-params.Th.P;end
if ~isfield(params,'lighting'),params.lighting='gouraud';end
if ~isfield(params,'alpha'),params.alpha=1;end 
if ~isfield(params,'ctx'),params.ctx='std';end
if ~isfield(params,'view'),params.view='lat';end

Anat.CtxL.data=zeros(size(Anat.CtxL.nodes,1),1);
Anat.CtxR.data=zeros(size(Anat.CtxR.nodes,1),1);

if ~isstruct(ROI),ROI.coord=ROI;end
Nroi=size(ROI.coord,1);
if ~isfield(ROI,'radius'), ROI.radius=repmat(5,[Nroi,1]);end
if ~isfield(ROI,'Network')
    ROI.Network=[[1:Nroi]',ones(Nroi,1)];
    Nnet=1;
else
    Nnet=max(ROI.Network(:,2));
end
if ~isfield(ROI,'color')
    if Nnet>1
        Cmap=jet(Nnet);
        ROI.color=zeros(Nroi,3);
        for j=1:Nnet
           keep=find(ROI.Network(:,2)==j);
           N=length(keep);
           ROI.color(keep,:)=repmat(Cmap(j,:),[N,1]);
        end
    else
        ROI.color=1.0.*ones(Nroi,3);
    end
end

% figCol='w'; % make fig white for easy printing. make black for good pptx


%% Draw cortex
PlotLRMeshes(Anat.CtxL,Anat.CtxR, params);
colorbar off
hold on


%% Draw ROIs
foci.lighting=params.lighting;
for j=1:Nnet
    keep=find(ROI.Network(:,2)==j);
    foci.color=ROI.color(keep,:);
    foci.radius=ROI.radius(keep,:);
    foci.location=ROI.coord(keep,:);
    foci.location=AdjustFoci(foci.location,Anat.CtxL,Anat.CtxR,params);
    hold on
    Draw_Foci(foci,10)
end
% set(gcf,'Color',figCol)