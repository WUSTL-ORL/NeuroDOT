function pos=AdjustFoci(foci,meshL,meshR,params)
%
% This function adjusts Foci locations to match the surface of the mesh
% in case it is inflated or rotated or otherwise adjusted.
% foci are the [x,y,z] locations of points to be re-drawn relative to the
% meshes.

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



%% Parameters and Initialization
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
    params.view = 'post';
end
if ~isfield(params, 'meshFoci')  ||  isempty(params.meshFoci)
    params.meshFoci = 0; % foci based on original space or mesh space
end


%% Determine home hemisphere
pos=foci;
[IDXl,dl]=knnsearch(meshL.nodes,pos);
[IDXr,dr]=knnsearch(meshR.nodes,pos);

if ~isfield(params, 'meshHome')    
    focL=dl<dr;
    focR=dr<=dl;
    
elseif strcmp(params.meshHome,'L')
    focL=ones(size(foci,1),1)==1;
    focR=focL==0;
    
elseif strcmp(params.meshHome,'R')
    focR=ones(size(foci,1),1)==1;
    focL=focR==0;
end

%% Choose inflation and adjust foci to match ctx
switch params.ctx
    case 'std'
        Lnodes=meshL.nodes;
        Rnodes=meshR.nodes;
    case 'inf'
        Lnodes=meshL.Inodes;
        Rnodes=meshR.Inodes;
    case 'vinf'
        Lnodes=meshL.VInodes;
        Rnodes=meshR.VInodes;
end

if params.meshFoci
    dxL=zeros(size(Lnodes));
    dxR=zeros(size(Lnodes));
else
    dxL=Lnodes-meshL.nodes;
    dxR=Rnodes-meshR.nodes;
end

pos(focL,:)=pos(focL,:)+dxL(IDXl(focL),:);
pos(focR,:)=pos(focR,:)+dxR(IDXr(focR),:);


%% Small adjustment to separate hemispheres
switch params.ctx
    case 'std'
        Lnodes=meshL.nodes;
        Rnodes=meshR.nodes;
    case 'inf'
        Lnodes=meshL.Inodes;
        Rnodes=meshR.Inodes;
pos(focL,1)=pos(focL,1)-max(Lnodes(:,1));
pos(focR,1)=pos(focR,1)-min(Rnodes(:,1));

Lnodes(:,1)=Lnodes(:,1)-max(Lnodes(:,1));
Rnodes(:,1)=Rnodes(:,1)-min(Rnodes(:,1));

    case 'vinf'
        Lnodes=meshL.VInodes;
        Rnodes=meshR.VInodes;
pos(focL,1)=pos(focL,1)-max(Lnodes(:,1));
pos(focR,1)=pos(focR,1)-min(Rnodes(:,1));

Lnodes(:,1)=Lnodes(:,1)-max(Lnodes(:,1));
Rnodes(:,1)=Rnodes(:,1)-min(Rnodes(:,1));
end


%% Rotate if necessary
if (strcmp(params.view ,'lat') || strcmp(params.view ,'med'))
dy=-5;
    % rotate right hemi around and move to position for visualization
    cmL=mean(Lnodes,1);
    cmR=mean(Rnodes,1);
    rm=rotation_matrix('z',pi);
    
    % Rotate
    switch params.view 
        case 'lat'
            Rnodes=(Rnodes-(repmat(cmR,size(Rnodes,1),1)))*rm +...
                (repmat(cmR,size(Rnodes,1),1));             
            pos(focR,:)=(pos(focR,:)-(repmat(cmR,sum(focR),1)))*rm +...
                (repmat(cmR,sum(focR),1));             
        case 'med'
            Lnodes=(Lnodes-(repmat(cmL,size(Lnodes,1),1)))*rm +...
                (repmat(cmL,size(Lnodes,1),1));       
            pos(focL,:)=(pos(focL,:)-(repmat(cmL,sum(focL),1)))*rm +...
                (repmat(cmL,sum(focL),1));              
    end
    pos(focR,1)=pos(focR,1)+(cmL(:,1)-cmR(:,1));    % Shift over to same YZ plane
    pos(focR,2)=pos(focR,2)-max(Rnodes(:,2))+min(Lnodes(:,2))+dy;
end