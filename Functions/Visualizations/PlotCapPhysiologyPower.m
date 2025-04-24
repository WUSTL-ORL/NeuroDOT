function [Plevels, SNR_DQ]=PlotCapPhysiologyPower(data, info, params)

% PlotCapPhysiologyPower A visualization of band limited power OR
% band-referenced SNR for each optode.
%
%
%   "params" fields that apply to this function (and their defaults):
%       fig_size    [20, 200, 1240, 420]    Default figure position vector.
%       fig_handle  (none)                  Specifies a figure to target.
%                                           If empty, spawns a new figure.
%       dimension   '2D'                    Specifies either a 2D or 3D
%                                           plot rendering.
%       rlimits     (all R2D)               Limits of pair radii displayed.
%       Nnns        (all NNs)               Number of NNs displayed.
%       Nwls        (all WLs)               Number of WLs averaged and
%                                           displayed.
%       useGM       0                       Use Good Measurements.
%       Cmap.P      'hot'                   Default color mapping.
%
% Dependencies: PLOTCAPDATA, ISTABLEVAR, APPLYCMAP.
%
% See Also: PLOTCAP, PLOTCAPGOODMEAS.
% 
% Copyright (c) 2017 Washington University 
% Created By: Adam T. Eggebrecht
% Eggebrecht et al., 2014, Nature Photonics; Zeff et al., 2007, PNAS.
%
% Washington University hereby grants to you a non-transferable, 
% non-exclusive, royalty-free, non-commercial, research license to use 
% and copy the computer code that is provided here (the “Software”).  
% You agree to include this license and the above copyright notice in 
% all copies of the Software.  The Software may not be distributed, 
% shared, or transferred to any third party.  This license does not 
% grant any rights or licenses to any other patents, copyrights, or 
% other forms of intellectual property owned or controlled by Washington 
% University.
% 
% YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS 
% PROVIDED “AS IS”, WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR 
% IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY 
% OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY 
% THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT.  
% IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON 
% UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR 
% CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH 
% THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER 
% IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS 
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

%% Parameters and Initialization.
if isfield(info, 'MEAS')
    if istable(info.MEAS)
        info.MEAS = table2struct(info.MEAS, 'ToScalar', true);
    end
end
if isfield(info, 'pairs')
    if istable(info.pairs)
        info.pairs = table2struct(info.pairs, 'ToScalar', true);
    end
end

LineColor = 'w';
BkgdColor = 'k';
Nm = size(info.pairs.Src,1);
Ns = length(unique(info.pairs.Src));
Nd = length(unique(info.pairs.Det));
cs = unique(info.pairs.WL); % WLs.
llfo=zeros(Ns+Nd,1);

dims = size(data);
Nt = dims(end);
NDtf = (ndims(data) > 2);
Nm = prod(dims(1:end-1));

if isfield(info.system,'init_framerate')
    fr=info.system.init_framerate;
else
    fr=info.system.framerate;
end

% Hardcoded params for color mapping to the units scaling used here.
% scale_order = 1e-2;
% base = 1e9;
% dr = 3;
params.mode = 'patch';
params.PD = 1;
params.Th.P = 0;
params.DR = 1000;

if (~isfield(params, 'Cmap')  ||  isempty(params.Cmap))...
        ||  (~isfield(params.Cmap, 'P')  ||  isempty(params.Cmap.P))
    params.Cmap.P = 'hot';
end
if ~isfield(params, 'dimension')  ||  isempty(params.dimension)
    params.dimension = '2D';
end
if ~isfield(params, 'rlimits') 
    params.rlimits=[10,30];
end
if ~isfield(params, 'WL') % Choose wavelength    
    if isfield(params, 'lambda') % Choose wavelength from actual lambda
        params.WL=info.pairs.WL(find(info.pairs.lambda==params.lambda,1));
    else
        params.WL=2;
    end
end
if ~isfield(params, 'Nwls')  ||  isempty(params.Nwls)
    params.Nwls = cs';
end
if ~isfield(params, 'useGM')  ||  isempty(params.useGM)
    params.useGM = 0;
end
if ~params.useGM  ||  ~isfield(info, 'MEAS')  ||  (isfield(info, 'MEAS')...
        &&  ~isfield(info.MEAS, 'GI'))
    GM = ones(Nm, 1);
else
    GM = info.MEAS.GI;
end
if ~isfield(params, 'dimension')  ||  isempty(params.dimension)
    params.dimension = '2D'; % '2D' | '3D'
end
if ~isfield(params, 'OD')  
    params.OD = 0; % 0 for raw data, 1 for Optical Density'
end
if ~isfield(params, 'type')  
    params.type = 'SNR'; % 0 for raw data, 1 for Optical Density'
end
if ~isfield(params, 'fig_size')  ||  isempty(params.fig_size)
    switch params.dimension
        case '2D'
            params.fig_size = [20, 200, 1240, 420];
        case '3D'
            params.fig_size = [20, 200, 560, 560];
    end
end
if ~isfield(params, 'calc_only')
    params.calc_only = 0;
end

if ~params.calc_only
if ~isfield(params, 'fig_handle')  ||  isempty(params.fig_handle)
    params.fig_handle = figure('Color', BkgdColor, 'Position', params.fig_size);
    new_fig = 1;
else
    switch params.fig_handle.Type
        case 'figure'
            set(groot, 'CurrentFigure', params.fig_handle);
        case 'axes'
            set(gcf, 'CurrentAxes', params.fig_handle);
    end
end
end

if ~isfield(params, 'freqs')  % Freq range to find peak
    params.freqs =[0.5,2.0]; % Pulse band
elseif ischar(params.freqs)
    switch params.freqs
        case 'pulse'
            params.freqs =[0.5,2.0]; % Pulse band            
        case 'fc'
            params.freqs =[0.009,0.08]; % FC band
    end
end
if ~isfield(params, 'freqsBW')  
   params.freqsBW=[0.009,0.08]; % range from which to determine bandwidth
end




%% N-D Input.
if NDtf
    data = reshape(data, [], Nt);
end

%% Use only time in synchpts if present
if isfield(info, 'paradigm')
    if isfield(info.paradigm, 'synchpts')
        NsynchPts = length(info.paradigm.synchpts); % set timing of data
        if NsynchPts > 2
            tF = info.paradigm.synchpts(end);
            t0 = info.paradigm.synchpts(2);
        elseif NsynchPts == 2
            tF = info.paradigm.synchpts(2);
            t0 = info.paradigm.synchpts(1);
        else
            tF = size(data, 2);
            t0 = 1;
        end
    else
        tF = size(data, 2);
        t0 = 1;
    end
else    
    tF = size(data, 2);
    t0 = 1;
end


%% Calculate power
if ~params.OD, [data, info.MEAS.Phi_0]=logmean(data);end
[ftmag, ftdomain] = fft_tts(data(:,t0:tF),fr);

keep=info.pairs.r2d>=params.rlimits(1,1) & ...
        info.pairs.r2d<=params.rlimits(1,2) & ...
        info.pairs.WL==params.WL & ...
        GM;
    
[~,idxFCm]=min(abs(ftdomain-params.freqsBW(1)));% Define freq indices    
[~,idxFCM]=min(abs(ftdomain-params.freqsBW(2)));    
BWfc=round((idxFCM-idxFCm)/2);                  % Define Bandwidth 
[~,idxPm]=min(abs(ftdomain-params.freqs(1)));    
[~,idxPM]=min(abs(ftdomain-params.freqs(2)));   

[mP,idxP1]=max(ftmag(keep,idxPm:idxPM),[],2); 
idxP1=round(mean(idxP1));  % find mean peak freq from all mean
floIdx=idxPm+idxP1-BWfc;
fhiIdx=idxPm+idxP1+BWfc;

if floIdx<1, floIdx=1;end
if fhiIdx > length(ftdomain), fhiIdx = length(ftdomain);end

Pmax=sum(ftmag(:,(floIdx):(fhiIdx)).^2,2); % sum pulse power
fNoise=setdiff([(idxPm-BWfc):(idxPM+BWfc)],(floIdx):(fhiIdx));
fNoise(fNoise<1)=[];
fNoise(fNoise>max(length(ftdomain)))=[];
Control=median(ftmag(:,fNoise).^2,2).*BWfc.*2; 
% Control=2.*BWfc.*median(ftmag(:,idxPm:end),2).^2;

% Plevels=1e5.*Pmax./Control;           % SNR of signal
% Plevels=1e0.*Pmax./Control;           % SNR of signal
switch params.type
    case 'SNR'
        Plevels=10.*log10(Pmax./Control); % SNR in dB
    case 'mag'
        Plevels=Pmax./max(Pmax(:));
end


%% Populate metric for visualizations
for s=1:Ns
    Sgood=keep & info.pairs.Src==s;
    if sum(Sgood)>0
    Cvalue=mean(Plevels(Sgood)); % Average across measurements
    else, Cvalue=1;
    end
    llfo(s)=Cvalue;
end
for d=1:Nd    
    Dgood=keep & info.pairs.Det==d;
    if sum(Dgood)>0
    Cvalue=mean(Plevels(Dgood)); % Average across measurementsl
    else, Cvalue=1;
    end
    llfo(Ns+d)=Cvalue;
end

%% Scaling and colormapping
M=max(llfo(:)); %% calculate the max SNR
param.Th.P=max([min(llfo),0]);
m=param.Th.P; %min value
% m=M/2;
params.Scale=M;%/2;
% llfo=llfo-M/2;

llfo_len=size(llfo);
llfo_Re=zeros(llfo_len);
for i=1:llfo_len
    if (llfo(i)<0)
        llfo_Re(i)=0;
    else
        llfo_Re(i)=llfo(i);
    end
    
end

Med=median(llfo_Re(:));% median value
Avg=mean(llfo_Re(:));%mean value

SNR_DQ.min_SNR = m;
SNR_DQ.max_SNR = M;
SNR_DQ.med_SNR = Med;
SNR_DQ.avg_SNR = Avg;

keepNN1 = info.pairs.NN == 1 & ...
        info.pairs.WL==2 & ...
        GM;
SNR_NN1 = Plevels(keepNN1);

SNR_DQ.min_SNR_NN1 = min(SNR_NN1); 
SNR_DQ.max_SNR_NN1 = max(SNR_NN1);
SNR_DQ.med_SNR_NN1 = median(SNR_NN1);
SNR_DQ.avg_SNR_NN1 = mean(SNR_NN1);

keepRsd = info.pairs.r3d <= 35 & ...
        info.pairs.r3d >=15 & ...
        info.pairs.WL==2 & ...
        GM;
SNR_Rsd = Plevels(keepRsd);
SNR_DQ.min_SNR_Rsd = min(SNR_Rsd); 
SNR_DQ.max_SNR_Rsd = max(SNR_Rsd);
SNR_DQ.med_SNR_Rsd = median(SNR_Rsd);
SNR_DQ.avg_SNR_Rsd = mean(SNR_Rsd);



if params.calc_only
    return
else
    params.Cmap.P='hot';
    colormap(hot)
    [SDRGB, CMAP] = applycmap(llfo, [], params);
    SrcRGB = SDRGB(1:Ns, :);
    DetRGB = SDRGB(Ns+1:end, :);
    PlotCapData(SrcRGB, DetRGB, info, params);
    pos=get(gca,'pos');




    %% Add Title and colorbar.
    switch params.type
        case 'SNR'
    tcell=[{['Mean Band-limited SNR']};...
        {['r', lower(params.dimension), ' \in ','[',...
        num2str(params.rlimits(1, 1)), ', ',...
        num2str(params.rlimits(1, 2)), '] mm; ',...
        'f\in ','[',num2str(params.freqs(1, 1)), ', ',...
        num2str(params.freqs(1, 2)), '] Hz']}];    
        if length(unique(info.pairs.NN)) > 1
            CB = colorbar('YTick', [0, 0.5, 1],...
            'YTickLabel',{[num2str(m)], 'SNR_d_B',[num2str(M)]},...
            'Color', LineColor,'Location', 'southoutside','position',...
            [pos(1)+pos(3)/3,pos(2)-0.01,pos(3)/3,pos(4)/35]);
        else
            CB = colorbar('YTick', [0, 0.5, 1],...
            'YTickLabel',{[num2str(m)], 'SNR_d_B',[num2str(M)]},...
            'Color', LineColor,'Location', 'southoutside');
            disp('Only one nearest-neighbor is present in the data');
        end


    title(tcell, 'Color', LineColor)

        case 'mag'
    tcell=[{['Mean FFT Power']};...
        {['r', lower(params.dimension), ' \in ','[',...
        num2str(params.rlimits(1, 1)), ', ',...
        num2str(params.rlimits(1, 2)), '] mm; ',...
        'f\in ','[',num2str(params.freqs(1, 1)), ', ',...
        num2str(params.freqs(1, 2)), '] Hz']}];    
    
        if length(unique(info.pairs.NN)) > 1
            CB = colorbar('YTick', [0, 0.5, 1],...
                'YTickLabel',{[num2str(m)], 'Max Relative Power',[num2str(M)]},...
                'Color', LineColor,'Location', 'southoutside','position',...
                 [pos(1)+pos(3)/3,pos(2)-0.01,pos(3)/3,0.01]);
        else
            CB = colorbar('YTick', [0, 0.5, 1],...
                'YTickLabel',{[num2str(m)], 'Max Relative Power',[num2str(M)]},...
                'Color', LineColor,'Location', 'southoutside');
                title(tcell, 'Color', LineColor)
    end
    end
end
end
    %
