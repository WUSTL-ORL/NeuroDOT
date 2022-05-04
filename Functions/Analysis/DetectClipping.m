function clipping=DetectClipping(data,params)
%
% This function looks for clipped channels in raw detector time
% traces based on params fields:
%   clipMinus   minimum allowed value
%   clipPlus    maximum allowed value
%   clipMu      maximum allowed mean value
%   clipRange   maximum allowed full data range
% It is assumed that data is NxT where T is time.

%% Parameters and Initialization
if~isfield(params,'clipM') % clipping max
    params.clipM=1; 
    params.clipMu=0.5; 
end 
if~isfield(params,'clipPlus') % clippling+
    params.clipPlus = +params.clipM ;
end
if~isfield(params,'clipMinus')% clippling-
    params.clipMinus = -params.clipM ;
end
if~isfield(params,'clipRange')%clipping range
    params.clipRange = 1 ;
end


%% Detect clipping
clipping=((sum(data<=params.clipMinus,2))+...
            (sum(data>=params.clipPlus,2))+...
            (abs(mean(data,2))>=params.clipMu)+...
            ((max(data,[],2)-min(data,[],2))>params.clipRange))>0;
        
% clipping=((sum(abs(d1)>=params.clipM,2))+... % old
%     (abs(mean(d1,2))>params.clipMu)+...
%     ((max(d1,[],2)-min(d1,[],2))>1))>0;