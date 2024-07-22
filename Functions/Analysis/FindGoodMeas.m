function info_out = FindGoodMeas(data, info_in, bthresh, params)

% FINDGOODMEAS Performs "Good Measurements" analysis.
%
%   info_out = FINDGOODMEAS(data, info_in) takes a light-level array "data"
%   in the MEAS x TIME format, and calculates the std of each channel
%   as its noise level. These are then thresholded by the default value of
%   0.075 to create a logical array, and both are returned as MEAS x 1
%   columns of the "info.MEAS" table. If pulse synch point information
%   exists in "info.system.synchpts", then FINDGOODMEAS will crop the data 
%   to the start and stop pulses.
%
%   info_out = FINDGOODMEAS(data, info_in, bthresh) allows the user to
%   specify a threshold value.
% 
% See Also: PLOTCAPGOODMEAS, PLOTHISTOGRAMSTD.
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

%% Parameters and Initialization.
if ~exist('info_in','var')
    info_out=struct;
else
    info_out = info_in;
end

if ~exist('params','var')
    params = struct;
end

if ~isfield(params,'GVwin')
    GVwin=600;
else
    GVwin = params.GVwin;
end

if ~isfield(info_out,'paradigm')
    info_out.paradigm=struct;
    Nt = size(data,2);
    info_out.paradigm.synchpts = [1, Nt];
    info_out.paradigm.synchtype = [1, 1];
    info_out.paradigm.Pulse_1 = [1, 2];
end
if ~exist('bthresh', 'var')
    bthresh = 0.075; % Empirically derived threshold value.
end 

omega_lp1 = 1;

if info_in.system.framerate/2 < omega_lp1 
    % Adjust Lowpass filter cutoff
    % frequency for systems with lower framerates
    omega_lp1 = (info_in.system.framerate/2)*0.90;
end



if ~isfield(params, 'r2d_dist')
    params.r2d_dist = 20;
end

if ~isfield(params, 'WL')
    params.WL = 2;
end
dims = size(data);
Nt = dims(end); % Assumes time is always the last dimension.
NDtf = (ndims(data) > 2);
if GVwin>(Nt-1), GVwin=(Nt-1);end


%% N-D Input.
if NDtf
    data = reshape(data, [], Nt);
end

%% Crop data to synchpts if necessary.
keep=(info_in.pairs.r2d<=params.r2d_dist) & (info_in.pairs.WL==params.WL);
foo=squeeze(data(keep,:)); 
foo=highpass(foo,0.02,info_in.system.framerate);% bandpass filter, omega_hp=0.02;    
foo=lowpass(foo,omega_lp1,info_in.system.framerate);% bandpass filter, omega_lp=1;    
foo=foo-circshift(foo,1,2);
foo(:,1)=0;
foob=rms(foo,1);
NtGV=Nt-GVwin;

if NtGV>1    % sliding window to grab a meaningful set for 'quiet'
    GVTD_win_means=zeros(NtGV,1);
    for j=1:NtGV
        GVTD_win_means(j)=mean(foob(j:(j+GVwin-1)));
    end    
    [~,t0]=min(GVTD_win_means);         % find min and set t0, tF
    tF=t0+GVwin-1;
    STD = std(data(:, t0:tF), [], 2);   % Calculate STD.
    
elseif isfield(info_out.paradigm, 'synchpts')
    NsynchPts = length(info_out.paradigm.synchpts); % set timing of data
    if NsynchPts > 2
        tF = info_out.paradigm.synchpts(end);
        t0 = info_out.paradigm.synchpts(2);
    elseif NsynchPts == 2
        tF = info_out.paradigm.synchpts(2);
        t0 = info_out.paradigm.synchpts(1);
    else
        tF = size(data, 2);
        t0 = 1;
    end
    STD = std(data(:, t0:tF), [], 2); % Calculate STD
else
    tF = size(data, 2);
    
    t0 = 1;
    STD = std(data, [], 2);
end

%% Populate in table of on-the-fly calculated stuff.
if (exist('t0','var') && exist('tf','var'))
    info_out.GVTDparams.t0=t0;
    info_out.GVTDparams.tF=tF;
end
if ~isfield(info_out,'MEAS')
    info_out.MEAS = struct;
    info_out.MEAS.STD = STD;
    info_out.MEAS.GI = STD <= bthresh;
else
    if istable(info_out.MEAS)
        info_out.MEAS = table2struct(info_out.MEAS, 'ToScalar', true);
    end
    info_out.MEAS.STD=STD;
    info_out.MEAS.GI=STD <= bthresh;
end
if isfield(info_out.MEAS,'Clipped')
    info_out.MEAS.GI=info_out.MEAS.GI & ~info_out.MEAS.Clipped; 
end


%
