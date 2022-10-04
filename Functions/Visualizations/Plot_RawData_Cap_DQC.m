function info=Plot_RawData_Cap_DQC(data,info,params)
%
% This function generates plots of data quality as related to the cap
% layout including: relative average light levels for 2 sets of distances
% of source-detector measurements, the cap good measurements plot, and 
% a measure of the pulse power at each optode location.
% 
%
% Copyright (c) 2017 Washington University 
% Created By: Adam T. Eggebrecht
% Eggebrecht et al., 2014, Nature Photonics; Zeff et al., 2007, PNAS.
%
% Washington University hereby grants to you a non-transferable, 
% non-exclusive, royalty-free, non-commercial, research license to use 
% and copy the computer code that is provided here (the Software).  
% You agree to include this license and the above copyright notice in 
% all copies of the Software.  The Software may not be distributed, 
% shared, or transferred to any third party.  This license does not 
% grant any rights or licenses to any other patents, copyrights, or 
% other forms of intellectual property owned or controlled by Washington 
% University.
% 
% YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS 
% PROVIDED AS IS, WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR 
% IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY 
% OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY 
% THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT.  
% IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON 
% UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR 
% CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH 
% THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER 
% IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS 
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.


%% Parameters and Initialization
if ~exist('params','var'), params=struct;end
if ~isfield(params,'bthresh'),params.bthresh=0.075;end


% Check if there is more than one NN - for determining rlimits and suplots
if length(unique(info.pairs.NN)) < 2
    numNNs = 1;
    if ~isfield(params,'rlimits')
        params.rlimits = [1,30];
    end
else
    numNNs = 2;
    if ~isfield(params,'rlimits')
        params.rlimits=[1,20;21,35;36,43];
    elseif size(params.rlimits,1)==1
        params.rlimits=cat(1,params.rlimits,[21,35;36,43]);
    elseif size(params.rlimits,1)==2
        params.rlimits=cat(1,params.rlimits,[36,43]);
    end
end

Rlimits=params.rlimits;
if ~isfield(params, 'mode'),params.mode = 'good';end
if ~isfield(params, 'useGM'),params.useGM = 1;end
if ~isreal(data),data=abs(data);end

params.fig_handle=figure('Units','Normalized',...
    'Position',[0.1,0.1,0.75,0.4],'Color','k');

%% Check for good measurements
if  ~isfield(info,'MEAS') || ~isfield(info.MEAS,'GI')
    info = FindGoodMeas(logmean(data), info, params.bthresh,params);
end

%% Mean signal level at each optode
if numNNs == 2
    subplot(4,2,1);
    params.rlimits=Rlimits(1,:);
    info.MEAS.Phi_o=PlotCapMeanLL(data, info, params);

    subplot(4,2,2);
    params.rlimits=Rlimits(2,:);
    PlotCapMeanLL(data, info, params);
else
    subplot(1,3,1)
    params.rlimits = Rlimits(1,:);
    info.MEAS.Phi_o= PlotCapMeanLL(data,info,params);
end


%% Good (and maybe bad) measurements
if numNNs == 2
    subplot(4,2,[5:8]);
    params.rlimits=[min(Rlimits(:)),max(Rlimits(:))];
    PlotCapGoodMeas(info, params);
else
    subplot(1,3,3);
    params.rlimits=[min(Rlimits(:)),max(Rlimits(:))];
    PlotCapGoodMeas(info, params);
end

%% Cap Physiology Plot
params=rmfield(params,'mode');
if numNNs == 2
    subplot(4,2,3);                             % Close neighborhood
    params.rlimits=Rlimits(1,:);
    info.MEAS.Pulse_SNR_R1=PlotCapPhysiologyPower(data, info, params);

    subplot(4,2,4);                             % 2nd neighborhood
    params.rlimits=Rlimits(2,:);
    info.MEAS.Pulse_SNR_R2=PlotCapPhysiologyPower(data, info, params);
else
    subplot(1,3,2);                             % Close neighborhood
    params.rlimits=Rlimits(1,:);
    info.MEAS.Pulse_SNR_R1=PlotCapPhysiologyPower(data, info, params);
end




