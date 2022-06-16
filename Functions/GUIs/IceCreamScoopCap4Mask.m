function maskCrop=IceCreamScoopCap4Mask(infoT1,tpos,maskOrig,maxOptodeDist)

% This function optimizes a mask for mesh generation based on the optode
% positions relative to the mask. Mask values beyond a given distance from
% any optode are set to zero for efficient light modeling.

%% Parameters and Initialization
if ~exist('maxOptodeDist','var')
    maxOptodeDist = 40; % max distance to keep in mm
end
% maxOptodeDist=maxOptodeDist./mean(abs(infoT1.mmppix));
Npos=size(tpos,1);

disp(['Setting up coordinate space'])
% Generate coordinates for voxels
nVxA = infoT1.nVx;
nVyA = infoT1.nVy;
nVzA = infoT1.nVz;
drA = infoT1.mmppix;
centerA = infoT1.center; 

% Create coordinates for each voxel index.
X = ((-centerA(1) + nVxA * drA(1)):(-drA(1)):(-centerA(1) + drA(1)))';
Y = ((-centerA(2) + nVyA * drA(2)):(-drA(2)):(-centerA(2) + drA(2)))';
Z = ((-centerA(3) + nVzA * drA(3)):(-drA(3)):(-centerA(3) + drA(3)))';


%% Find min distance to optodes for each voxel
[x,y,z] = ndgrid(X,Y,Z); 
voxXYZ = [x(:),y(:),z(:)]; 

OptVoxDist=ones(nVxA,nVyA,nVzA).*100.*maxOptodeDist;
OptVoxDist=OptVoxDist(:);
disp(['Calculating pad distances'])
for j=1:Npos
    d=pdist2(tpos(j,:),voxXYZ,'Euclidean')';
    OptVoxDist=min(d,OptVoxDist);
end
minDist = reshape(OptVoxDist,nVxA,nVyA,nVzA); 


%% Kill mask beyond max distance
disp(['Completing mask'])
distMask=(minDist < maxOptodeDist);
maskCrop = (distMask & (maskOrig>0)) .* maskOrig; 