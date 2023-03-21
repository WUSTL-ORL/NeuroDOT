function Plot_RawData_Time_Traces_Overview(data,info,params)
%
% This function generates a plot of raw data time traces separated by
% wavelength (columns). The top row shows time traces for all measurements
% within a source-detector distance range (defaults as 0.1 - 5.0 cm). The
% bottom row shows the same measurements but including only measurements
% passing a variance threshold (default: 0.075) as well as vertical lines
% correspondong to the stimulus paradigm.
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
Nwl=length(unique(info.pairs.WL));
[Nm,Nt]=size(data);
if ~isreal(data),data=abs(data);end
if isfield(info.system,'init_framerate')
    fr=info.system.init_framerate;
else
    fr=info.system.framerate;
end
    
dt=1/fr;                          
t=[1:1:Nt].*dt;

if ~exist('params','var'), params=struct;end
if ~isfield(params,'bthresh'),params.bthresh=0.075;end
if ~isfield(params,'rlimits'),params.rlimits=[1,40];end
if ~isfield(params, 'yscale'),params.yscale = 'log';end
lambdas=unique(info.pairs.lambda,'stable');
WLs=unique(info.pairs.WL,'stable');

params.fig_handle=figure('Units','Normalized',...
    'Position',[0.1,0.1,0.8,0.8],'Color','k');

%% Check for good measurements
if ~isfield(info,'MEAS') || ~isfield(info.MEAS,'GI')
    lmdata=logmean(data);
    info = FindGoodMeas(lmdata, info, params.bthresh,params);
end


%% Top row: all measurements broken apart by wavelength
for j=1:Nwl
    subplot(2,Nwl,j)
    keep=info.pairs.r2d>=params.rlimits(1,1) & ...
        info.pairs.r2d<=params.rlimits(1,2) & ...
        info.pairs.WL==j;
    PlotTimeTraceData(data(keep,:), t, params);
    xlim([0,max(t)+1]);
    xlabel('Time [sec]');ylabel('\Phi');
    if isfield(info.pairs,'lambda')
        title(['All ',num2str(lambdas(j)),' nm, Rsd:',...
            num2str(params.rlimits(1,1)),'-',...
            num2str(params.rlimits(1,2)),' mm'],...
            'Color','w')
    else
        title(['All WL ## ',num2str(WLs(j)),' nm, Rsd:',...
            num2str(params.rlimits(1,1)),'-',...
            num2str(params.rlimits(1,2)),' mm'],...
            'Color','w')
    end
        
end

%% Bottom row: good measurements broken apart by wavelength
for j=1:Nwl
    subplot(2,Nwl,j+Nwl)
    keep=info.pairs.r2d>=params.rlimits(1,1) & ...
        info.pairs.r2d<=params.rlimits(1,2) & ...
        info.pairs.WL==j & ...
        info.MEAS.GI;
    PlotTimeTraceData(data(keep,:), t, params);
    xlim([0,max(t)+1]);
    xlabel('Time [sec]');ylabel('\Phi')
    if isfield(info.pairs,'lambda')
        title(['Good ',num2str(lambdas(j)),' nm, Rsd:',...
            num2str(params.rlimits(1,1)),'-',...
            num2str(params.rlimits(1,2)),' mm'],...
            'Color','w')
    else
        title(['Good WL ## ',num2str(WLs(j)),' nm, Rsd:',...
            num2str(params.rlimits(1,1)),'-',...
            num2str(params.rlimits(1,2)),' mm'],...
            'Color','w')
    end    
    % Divide experimental paradigm timing by framerate if not done
if isfield(info,'paradigm') % Add in experimental paradigm timing
    if max(info.paradigm.synchpts) <= size(data,2)/info.system.framerate
       DrawColoredSynchPoints(info,0);
    else
       DrawColoredSynchPoints(info,1);

    end
end

end


