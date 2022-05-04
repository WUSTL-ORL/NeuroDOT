function [GVTD]=CalcGVTD(data)

% This function calculates the Root Mean Square across measurements (be
% they log-mean light levels or voxels) of the temporal derivative. The
% data is assumed to have measurements in the first and time in the last 
% dimension. Any selection of measurement type or voxel index must be done
% outside of this function.

%% Double check data has correct dimensions
Dsizes=size(data);
Ndim=length(Dsizes);
if Ndim>2, data=reshape(data,[],Dsizes(end));end
 
%% 1st Temporal Derivative
Ddata=data-circshift(data,[0 -1]);

%% RMS across measurements
GVTD=cat(1,0,rms(Ddata(:,1:(end-1)),1)');