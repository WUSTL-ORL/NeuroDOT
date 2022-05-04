function nirs = converter_ND2_to_HOMER(data, info, aux)%, infoA)

% CONVERTER_ND2_TO_HOMER Converts workspace from ND2 to HOMER.
%
%   nirs = CONVERTER_ND2_TO_HOMER(data, info, aux) takes MEAS x TIME array
%   "data", structure "info", and optionally the "aux" structure, and
%   converts them to a "nirs" structure that can be saved into a MAT file
%   with the ".nirs" extension, and is compatible with HOMER2.
% 
% See Also: CONVERTER_HOMER_TO_ND2.
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
nirs = [];

%% Check for "info.io.nirs".
% Save all original data, then write over with any data that may have
% changed.
if isfield(info, 'io')  &&  isfield(info.io, 'nirs')
    nirs = info.io.nirs;
else % If no original data...
    
    % After looking through their source code, specifically on file IO,
    % HOMER doesn't really mind if things are missing and fills them in
    % along the way. IE, it's generally more forgiving than NIFTI. So
    % instead of filling in needless fields, we'll only fill in what we
    % have.
    
    if exist('aux', 'var')
        if isstruct(aux)
            nirs.aux = struct2array(aux);
        elseif isnumeric(aux)
            nirs.aux = aux;
        end
    end
    
    nirs.SpatialUnit = 'mm'; % We at least share this with HOMER.
    
    % Exclude this for now. May need to wrap into HOMER coordinates
    % later, like we did for NIFTI.
%     for var = {'xmin', 'ymin', 'xmax', 'ymax'}
%         if exist('infoA', 'var')  &&  isfield(infoA, var{1})
%             nirs.SD.(var{1}) = infoA.(var{1});
%         else
%             nirs.SD.(var{1}) = [];
%         end
%     end
    
end

%% Save data.
nirs.d = data';
[Nm, Nt] = size(data);

%% Save time array.
% If there was originally a time array, we'll save the last value
% because we may need it for later.
if isfield(nirs, 't')
    init_t_end = nirs.t(end);
end
nirs.t = 0:1 / info.system.framerate:(Nt - 1) / info.system.framerate';

%% Save s.
Npulse = numel(fieldnames(info.paradigm)) - 3;
nirs.s = zeros(Nt, Npulse);
for k = 1:Npulse
    nirs.s(info.paradigm.synchpts(info.paradigm.(['Pulse_' num2str(k)])), k) = k;
end
% We're not going to overwrite "s0" if it's already there.
try
    if ~isfield(nirs, 's0')
        if ~exist('init_t_end', 'var')
            init_t_end = round(Nt ./ info.system.init_framerate);
        end
        nirs.s0 = zeros(init_t_end, Npulse);
        for k = 1:Npulse
            nirs.s0(info.paradigm.init_synchpts(round(info.paradigm.(['Pulse_' num2str(k)])...
                * (info.system.init_framerate / info.system.framerate))), k) = k;
        end
    end
catch err
%     rethrow(err)
    disp('*** WARNING: "s0" interpolation failed. ***')
end

%% Save MeasList and ml.
nirs.SD.MeasList = [info.pairs.Src, info.pairs.Det, ones(Nm, 1), info.pairs.WL];
nirs.ml = nirs.SD.MeasList;

%% Save rest of SD.
nirs.SD.Lambda = unique(info.pairs.lambda);
nirs.SD.SrcPos = info.optodes.spos3;
nirs.SD.DetPos = info.optodes.dpos3;
nirs.SD.nSrcs = length(unique(info.pairs.Src));
nirs.SD.nDets = length(unique(info.pairs.Det));
if isfield(info, 'MEAS')  &&  istablevar(info.MEAS, 'GI')
    nirs.SD.MeasListAct = info.MEAS.GI;
end

%% Save CondNames, if not already there.
if ~isfield(nirs, 'CondNames')
    pulses = rmfield(info.paradigm, {'synchpts', 'synchtype', 'init_synchpts'});
    nirs.CondNames = fieldnames(pulses);
end



%
