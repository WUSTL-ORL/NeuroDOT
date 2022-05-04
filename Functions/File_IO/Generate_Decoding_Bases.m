function [H,S,C]=Generate_Decoding_Bases(Nsamples,sr,freq,params)
%
% This function generates the Hanning Window (H), Sin (S), and Cos (C) 
% basis vectors for frequency multiplexed decoding. The inputs are the 
% number of samples in the window, the sampling rate (in Hz), and the 
% frequency of modulation (in Hz).
% The older decoding was with a Hamming window. That window allows for
% transients at the edge to remain in the data. As of 180124, we are moving
% to a Hanning window which goes to zero at the endpoints of the window.
% 

%% Parameters and Initialization
if ~exist('params','var'),params.Wtype='Hann';end
if strcmp(params.Wtype,'Bcar')
    if ~isfield(params,'BcarPts')
        params.BcarPts=ones(max(Nsamples),1);
    end
end
    
Nsamples=ceil(Nsamples);
Nf=length(freq);
if length(Nsamples)==1
    Nsamples=ones(size(freq)).*Nsamples;
end

%% Make window and bases
for j=1:Nf
    switch params.Wtype
        case 'Hamm'
H{j}=0.54-(0.46*(cos((2*pi.*[0:(Nsamples(j)-1)])./(Nsamples(j)-1)))); %Hamm
        case 'Hann'
H{j}=0.5.*(1-cos((2*pi.*[0:(Nsamples(j)-1)])./(Nsamples(j)-1))); %Hann
        case 'Bcar'
H{j}=params.BcarPts(1:Nsamples(j)); %BoxCar
    end
H{j}=H{j}./sqrt(mean(H{j}.^2)); 
S{j}=sin((2*pi*freq(j).*[0:(Nsamples(j)-1)])/sr).*(sqrt(2)/Nsamples(j));
C{j}=cos((2*pi*freq(j).*[0:(Nsamples(j)-1)])/sr).*(sqrt(2)/Nsamples(j));
end