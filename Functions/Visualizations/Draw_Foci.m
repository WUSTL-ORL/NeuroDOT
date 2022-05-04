function Draw_Foci(foci,faces)

% This function draw spheres of a specified color and radius at the
% location given by the structure foci. Fields include: location, color,
% radius. Only foci.location is strictly required as color and radius each
% have default values ([0,0,0], and 5, repectively).
%   foci.color are rows of RGB values contained within [0,1].
%   foci.lighting can be any of 'phong' (default),'flat'


%% Parameters and Initialization
if ~exist('faces','var'), faces=100;end
if ~isfield(foci,'lighting'),foci.lighting='phong';end
if ~isfield(foci,'AmbientStrength'),foci.AmbientStrength=0.02;end

[x,y,z]=sphere(faces);
Nfoci=size(foci.location,1);

if ~isfield(foci,'color'),foci.color=zeros(Nfoci,3);end
if size(foci.color,1)~=Nfoci,foci.color=repmat(foci.color,[Nfoci,1]);end

if ~isfield(foci,'radius'),foci.radius=ones(Nfoci,1).*5;end
if size(foci.radius,1)~=Nfoci,foci.radius=repmat(foci.radius,[Nfoci,1]);end


%% Draw the foci
for j=1:Nfoci
    
    rad=foci.radius(j);
    
    hh=patch(surf2patch(rad*x+foci.location(j,1),...
        rad*y+foci.location(j,2),rad*z+foci.location(j,3)),...
        'EdgeColor',foci.color(j,:),'FaceColor',foci.color(j,:),...
        'EdgeAlpha',0);
    
    set(hh,'FaceLighting',foci.lighting,'AmbientStrength',foci.AmbientStrength);
    
end

axis image