function [M]=FFT_Decode_Raw_Acqdecode(data,windowB,Ndft,idxF,div)
%
% % data may be multiple 'time steps' in 1st dimension.
%

%% data --> windowed data
if ~exist('div','var'),div=1;end
data=bsxfun(@minus,data,mean(data,2));  % remove mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% already done in Generate_Decoding_Bases_U87
% windowB=windowB./sqrt(mean(windowB.^2));% not in AD 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wData=bsxfun(@times,windowB,data);

%% Take FFT
Y = fft(wData,Ndft,2)/Ndft;                 % fft
SSfft=sqrt(2)*abs(Y(:,1:floor(Ndft/2)+1));  % Single Sided Amplitude
if ( ~iscell(idxF) )
    % Compute magnitude at just 1 frequency
    M=sqrt(sum((SSfft(:,idxF).^2),2)).*div; 
else
    % Compute magnitude at just N frequencies using indices in idxF
    for i = 1:numel(idxF)
        M(i,:)=sqrt(sum((SSfft(:,idxF{i}).^2),2)).*div(i); 
    end
end
