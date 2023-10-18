function Plot_RawData_Signal_and_Noise(data,info,params)
%
% This function generates a plot of the temporal standard deviation for
% each measurement as a function of mean light level (Phi_0). Also shown is
% the temporal standard deviation of the light attenuation as given by
% y=-ln[Phi/Phi_0]. If the data are complex, as in the case of
% frequency domain data, an additional 2 plots are shown: standard
% deviation in the phase as a function of mean light level and the standard
% deviation in phase as a function of the standard deviation in
% attenuation.
% 
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

Nwl=length(unique(info.pairs.WL));
[Nm,Nt]=size(data);
if ~exist('params','var'), params=struct;end
if ~isfield(params,'bthresh'),params.bthresh=0.075;end
if ~isfield(params,'rlimits'),params.rlimits=[1,40];end
if ~isfield(params, 'yscale'),params.yscale = 'log';end
lambdas=unique(info.pairs.lambda,'stable');
WLs=unique(info.pairs.WL,'stable');



%% Calculate logmean and stds
lmdata=logmean(data);
Nlmm=size(lmdata,1);
Phi_0=mean(abs(data),2);
std_Phi_0=std(abs(data),[],2);
std_Y=std(lmdata,[],2);


%% Check for good measurements
if ~isfield(info,'MEAS') || ~isfield(info.MEAS,'GI')
    info = FindGoodMeas(lmdata, info, params.bthresh);
end


%% Make figures
switch Nm==Nlmm
    
    case 1 % CW
        params.fig_handle=figure('Units','Normalized',...
            'Position',[0.1,0.1,0.6,0.4]);
        
        subplot(1,2,1);
        loglog(Phi_0(~info.MEAS.GI),std_Phi_0(~info.MEAS.GI),'ko')
        xlabel(['\Phi_0']);ylabel(['\sigma (|\Phi|)']);hold on
        loglog(Phi_0(info.MEAS.GI),std_Phi_0(info.MEAS.GI),'b*')
        
        Ral_x=logspace(-6,2,50); % Raleigh-distributed data: snr=mu/std~1.91
        Ral_y=Ral_x./1.91;
        loglog(Ral_x,Ral_y,'k-'); 

        
        
        subplot(1,2,2);
        loglog(Phi_0(info.MEAS.GI),std_Y(info.MEAS.GI(1:Nm)),'b*');hold on;
        loglog(Phi_0(~info.MEAS.GI),std_Y(~info.MEAS.GI(1:Nm)),'ko');hold on;
        xM=ceil(max(log10(Phi_0)));
        xm=floor(min(log10(Phi_0)));
        loglog([10^xm,10^xM],[1,1].*0.075,'r')
        xlabel(['\Phi_0']);ylabel(['\sigma (|Y|)']);hold on        
        
        
    case 0 % fd
        params.fig_handle=figure('Units','Normalized',...
            'Position',[0.1,0.1,0.6,0.8]);
        
        std_Yatt=std_Y(1:Nm);
        std_Yph=std_Y((Nm+1):end);
        
        subplot(2,2,1);
        loglog(Phi_0,std_Phi_0,'k*')
        xlabel(['\Phi_0']);ylabel(['\sigma (|\Phi|)']);hold on
        loglog(Phi_0(info.MEAS.GI(1:Nm)),std_Phi_0(info.MEAS.GI(1:Nm)),'ro')
        
        subplot(2,2,2);
        loglog(Phi_0,std_Yatt,'k*');hold on;
        xM=ceil(max(log10(Phi_0)));
        xm=floor(min(log10(Phi_0)));
        loglog([10^xm,10^xM],[1,1].*0.075,'r')
        xlabel(['\Phi_0']);ylabel(['\sigma (|Y|)']);hold on
         
        subplot(2,2,3);
        loglog(Phi_0,std_Yph,'k*');
        xlabel(['\Phi_0']);ylabel(['\sigma (\theta) [Radians]']);
        
        subplot(2,2,4);
        loglog(std_Yatt,std_Yph,'k*');
        xlabel(['\sigma (|Y|)']);ylabel(['\sigma (\theta) [Radians]']);hold on
        loglog(std_Yatt(info.MEAS.GI(1:Nm)),std_Yph(info.MEAS.GI(1:Nm)),'ro')
        
        
end