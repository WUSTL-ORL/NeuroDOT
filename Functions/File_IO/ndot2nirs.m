function nirsData = ndot2nirs(data, info)

% Translation from NeuroDOT-compatible format to .nirs 
% This function takes in data and an info structure in NeuroDOT format
% and converts it to the .nirs format.
%
% Inputs:
%       'data': time series data arranged by [#channels x #samples]
%       'info': NeuroDOT formatted info structure
%
% Outputs
%       'nirsData': .nirs formatted data structure containing the following
%               'd': time series data arranged by [#time points x #channels]
%               'SD': structure containing the source and detector
%               information in addition to the data measurement list
%               't': time point array
%               's': time points and stimulus onsets
%               'ml': the data measurement list
%        
% Copyright (c) 2017 Washington University 
% Created By: Ari Segel and Emma Speh
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
   

%% Data
nirsData.d = data'; 


%% SD
if isfield(info, 'io')
    nirsData.SD.nDets = info.io.Nd;
    nirsData.SD.nSrcs = info.io.Ns;
else
    nirsData.SD.nDets = size(info.optodes.dpos3,1);
    nirsData.SD.nSrcs = size(info.optodes.spos3,1);
end

nirsData.SD.Lambda = unique(info.pairs.lambda)';


%% Time 't'
if isfield(info, 'misc')
    if isfield(info.misc, 'startTime')
        if info.misc.startTime == 0   
            nirsData.t = ((0:size(nirsData.d,1)-1)./info.system.framerate)';
        else
            nirsData.t = ((1:size(nirsData.d,1))./info.system.framerate)';
        end
    else
        nirsData.t = ((1:size(nirsData.d,1))./info.system.framerate)';
    end
else
    nirsData.t = ((1:size(nirsData.d,1))./info.system.framerate)';
end


%% Stimulus 's'
if isfield(info, 'paradigm')
    fields = fieldnames(info.paradigm); 
    fieldlist = [];
    for j = 1:length(fields)
        if strlength(fields(j)) >= 4
            fieldlist = [fieldlist, fields(j)];
        else
            fieldlist = [fieldlist, 'not_stim'];
        end
    end
    idxPulse = ismember(cellfun(@(x) x(1:5), fieldlist, 'UniformOutput', false), 'Pulse');
    pulses = fields(idxPulse);
    pulses = sort(pulses);
    num_synchs = size(pulses,1);
    nirsData.s = zeros(size(nirsData.t,1),num_synchs); 
    for j = 1:num_synchs
        pulse = string(pulses(j));
        pulseidx = zeros(1, length(pulse));
        for k = 1:length(info.paradigm.(pulse))  
            [~,pulseidx(k)] = min(abs(nirsData.t' - (info.paradigm.synchpts(info.paradigm.(pulse)(k))./info.system.framerate)));
        end
        nirsData.s(pulseidx,j) = 1;
    end
end


%% Optodes
if isfield(info.optodes, 'spos2')
    nirsData.SD.SrcPos2 = info.optodes.spos2;
end
if isfield(info.optodes, 'spos3')
    nirsData.SD.SrcPos3 = info.optodes.spos3;
end
if isfield(info.optodes, 'dpos2')
    nirsData.SD.DetPos2 = info.optodes.dpos2;
end
if isfield(info.optodes, 'dpos3')
    nirsData.SD.DetPos3 = info.optodes.dpos3;
end


%% MeasList
nirsData.SD.MeasList(:,1) = info.pairs.Src; 
nirsData.SD.MeasList(:,2) = info.pairs.Det; 
nirsData.SD.MeasList(:,3) = 1; 
nirsData.SD.MeasList(:,4) = info.pairs.WL;

nirsData.ml = nirsData.SD.MeasList;


%% Aux
if isfield(info, 'misc')
    if isfield(info.misc, 'aux')
        nirsData.aux = info.misc.aux.dataTimeSeries;
    end
end


end
