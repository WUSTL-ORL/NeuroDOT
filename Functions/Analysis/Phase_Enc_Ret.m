function [Mag,Phase,SigMag,freqs,y,power]=Phase_Enc_Ret(data,fs,StimF)

% This function calculates the magnitude and phase of a set of time trace
% data.  Additionally, it calculates the significant magnitudes as the
% ratio of the power at the stimulus frequency to the power in all other
% frequencies (other than harmonics of stim freq and very low frequencies.
% It is assumed that time is in the 2nd (and last) dimension.




%% Parameters and initialization
dims = size(data);
Nt = dims(end); % Assumes time is always the last dimension.
NDtf = (ndims(data) > 2);
if NDtf
    data = reshape(data, [], Nt);
end

[Nv,Nt]=size(data);             % Window length
dt=1/fs;                        % Time increment per sample
t=(0:(Nt-1)).*dt;               % Time range for data
n=Nt;                           % Length of transform
f=(0:n-1)*(fs/n);               % Frequency range
freqs=f(1:floor(n/2));          % 0-Nyquist Freq set of freqs
bkgnd=sum(abs(data),2)==0;


%% Find Stimulus and Noise components of frequency band
harms=StimF:StimF:max(freqs);
for i=1:length(harms)
    [m1,mI1(i)]=min(abs(freqs-harms(i)));
end
Sig=mI1(1);                     % Stimulus Frequency
mI1=cat(2,mI1,mI1+1,mI1-1);
mI1=unique(mI1);
mI1=[find(freqs<0.01),mI1];
Fnoise=1:length(freqs);
Fnoise(mI1)=[];                 % Set of noise frequency indices


%% Calculate FFT and get magnitude and phase
y=fft(data',n);                 % DFT
y=y(1:floor(n/2),:).';  
power= y.*conj(y)/n;            % Power of DFT
Mag=power(:,Sig);               % Magnitude at stim freq
Phase=angle(y(:,Sig));          % Phase at stim freq
Noise=sum(power(:,Fnoise),2);
SigMag=fcdf(Mag./Noise,Nt,Nt,'upper');      % these are p-values
if any(SigMag)
SigMag(SigMag==0)=min(SigMag(SigMag~=0));   % correct for very sig values
SigMag=-log10(SigMag);          % Convert to linear scale in log10(p-value)
end

%% Correct for possible background
Mag(bkgnd)=0;
Phase(bkgnd)=0;
SigMag(bkgnd)=0;

%% Re-pack data
if NDtf
    Mag = reshape(Mag, dims(1:(end-1)));
    Phase = reshape(Phase, dims(1:(end-1)));
    SigMag = reshape(SigMag, dims(1:(end-1)));
end
