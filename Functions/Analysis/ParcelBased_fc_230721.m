function [fcMaps,fcMatrix,ParcelTT]=ParcelBased_fc_230721(data,Parcels,dim,params)
%
% This function calculates Parcel Based functional connectivity between a
% volumetric data set (data) and a set of volumetrically defined parcels. 
% To increase computational efficiency, the data are assumed to be in 
% 2 dimensions with time in the 2nd dimension. If data are not in 2D, 
% data are reshaped to be in 2D. 
% Parcels are assumed to be in the same space, with both data and Parcels
% described by the dim structure.
% Parcels are assumed to be a volume of zeros (background) and integers
% defining parcel locations. The Parcel labels need not be contiguous.
% The outputs are the Fisher-z transformed maps (fcMaps), 
% the Fisher-z transformed fc matrix (fcMatrix), 
% and the ROI timecourses (seedTT). 
% A temporal mask delineating time points to keep may also be passed in as
% params.tMask.


%% Parameters and initialization
dims = size(data);
Nt = dims(end); % Assumes time is always the last dimension.
NDtf = (ndims(data) > 2);
if NDtf
    data = reshape(data, [], Nt);
end

if ~exist('params','var'), params=struct;end
if ~isfield(params,'tMask'),params.tMask=ones(Nt,1);end

% Normalize for fast fc calc
 data=normr(data(:,params.tMask==1));


uParcels=unique(Parcels(Parcels~=0));
NParcels=length(uParcels);
ParcelTT=zeros(Nt,NParcels);
fcMaps=zeros(dims(1)*dims(2)*dims(3),NParcels);


%% Calculate fc
disp(['Calculating correlations'])
for k=1:NParcels
    vol=(Parcels(:)==uParcels(k));
    Ns=nansum(vol(:));
    % Generate Seed Time Trace
    ParcelTT(:,k)=([nansum(bsxfun(@times,data,vol),1)./Ns]');
    % Generate Seed map
    fcMaps(:,k)=FisherR2Z((((ParcelTT(:,k))'*(data)')'));
end
fcMaps=reshape(fcMaps,dims(1),dims(2),dims(3),NParcels);
fcMatrix=FisherR2Z(corrcoef(ParcelTT));