function hh=Draw_Foci_191203(foci,params)

% This function draw spheres of a specified color and radius at the
% location given by foci. Params fields include: color, faces
% radius. Only foci is strictly required as color and radius each
% have default values ([0,0,0], and 5, repectively).
%   params.color are rows of RGB values contained within [0,1].
%   params.lighting can be any of 'phong' (default),'flat'


%% Parameters and Initialization
if ~exist('params','var'), params=struct;end
if ~isfield(params,'faces'),   params.faces=10;end
if ~isfield(params,'lighting'),params.lighting='phong';end
if ~isfield(params,'AmbientStrength'),params.AmbientStrength=0.02;end
if ~isfield(params,'SpecularStrength'),params.SpecularStrength=0.9;end
if ~isfield(params,'DiffuseStrength'),params.DiffuseStrength=0.6;end
if ~isfield(params,'FaceAlpha'),params.FaceAlpha=1;end

Nfoci=size(foci,1);

if ~isfield(params,'color'),params.color=zeros(Nfoci,3);end
if size(params.color,1)~=Nfoci,params.color=repmat(params.color,[Nfoci,1]);end
if ~isfield(params,'radius'),params.radius=ones(Nfoci,1).*5;end
if size(params.radius,1)~=Nfoci,params.radius=repmat(params.radius,[Nfoci,1]);end

[x,y,z]=sphere(params.faces);

%% Draw the foci
for j=1:Nfoci
    
    rad=params.radius(j);
    
    hh=patch(surf2patch(rad*x+foci(j,1),rad*y+foci(j,2),rad*z+foci(j,3)),...
        'EdgeColor',params.color(j,:),'FaceColor',params.color(j,:),...
        'EdgeAlpha',0,'FaceAlpha',params.FaceAlpha);
    
    set(hh,'FaceLighting',params.lighting,...
        'AmbientStrength',params.AmbientStrength,...
        'SpecularStrength',params.SpecularStrength,...
        'DiffuseStrength',params.DiffuseStrength);    
end

axis image