function GVidx=Coord2GoodVox(inVox,info,dim)
%
% Inputs
%   inVox    list of coordinates of interest in original volume
%   info     structure describing volume space
%   dim      structure defining A-matrix that Good_Vox lives in
% Outputs
%           GVidx are the dim.Good_Vox indices relating to the input
%           coordinates. If the input voxel is not in dim.Good_Vox, the
%           entry should be empty.
%
x=change_space_coords(inVox,info,'coord');
xDim=change_space_coords(x,dim,'idxC');
VolIdx=sub2ind([dim.nVx,dim.nVy,dim.nVz],xDim(:,1),xDim(:,2),xDim(:,3));
[Ia,GVidx]=ismember(VolIdx,dim.Good_Vox);

