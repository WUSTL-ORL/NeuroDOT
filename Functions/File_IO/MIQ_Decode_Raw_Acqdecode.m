function [M,I,Q]=MIQ_Decode_Raw_Acqdecode(data,sinB,cosB,windowB,div)
%
% % data may be multiple 'time steps' in 1st dimension.
%

%% data --> windowed data
if ~exist('div','var'),div=1;end
Nts=size(data,1);
wData=bsxfun(@times,windowB,data);

%% I, Q, M
I=dot(repmat(cosB,[Nts,1]),wData,2);
            
Q=dot(repmat(sinB,[Nts,1]),wData,2);
            
M=sqrt(I.^2+Q.^2).*div;