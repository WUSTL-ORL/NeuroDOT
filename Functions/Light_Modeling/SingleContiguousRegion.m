function R1=SingleContiguousRegion(data)

% This function segments out the largest continuous region of non-zero
% values from a data set.  All other non-zero voxels are set to zero.
% Output volume has 1 contiguous region of activation with voxel values
% equal to input data values.

% Binarize data set
bw_data=data;
bw_data(bw_data<0)=0;
bw_data(bw_data>0)=1;


% Separate largest contiguous region
L=bwlabeln(bw_data);
stats=regionprops(L,'Area');
A=[stats.Area];
biggest=find(A==max(A));
bw_data(L~=max(biggest))=0;


% Finalize output data
R1=data.*bw_data;