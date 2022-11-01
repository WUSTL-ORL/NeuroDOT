function snirf2ndot(filename, output, type)
%
% snirf2ndot(filename,output,type) takes a file with the '.snirf' extension
% in the SNIRF format and converts it to NeuroDOT formatting. 

% 'Filename' is the name of the file to be converted, followed by the 
% .snirf extension.

% 'Output' is the filename (without extension) to be saved in NeuroDOT
% format. Output files are saved as '.mat.' 
% 'Type' is an optional input -
% the only currently acceptable value for this field is 'snirf.'

% Convert .snirf data format to be used in NeuroDOT
% Only the fields that are needed to run NeuroDOT are extracted

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

if ~exist('type','var')
    type = 'snirf';
end

if ~exist('output','var')
    output = filename;
end
    

snf = loadsnirf(filename);

if type == 'snirf' 
    if isfield(snf.nirs, 'data')
        data = snf.nirs.data.dataTimeSeries';% Snirf measurement list describes each channel as a column of the data
        if isfield(snf.nirs.data, 'time')
            info.system.framerate = 1/(mean(diff(snf.nirs.data.time)));
            info.system.init_framerate = info.system.framerate;
        end
    end
  
       
    if exist('snf','var')
        if ~isfield(snf.nirs, 'metaDataTags')
        else
            if ~isfield(snf.nirs.metaDataTags, 'Nd') 
            else
                info.io.Nd = snf.nirs.metaDataTags.Nd; %custom NeuroDOT field
            end
            if ~isfield(snf.nirs.metaDataTags, 'Ns')
            else
                info.io.Ns = snf.nirs.metaDataTags.Ns; %custom NeuroDOT field
            end
            if ~isfield(snf.nirs.metaDataTags, 'Nwl')
            else
                info.io.Nwl = snf.nirs.metaDataTags.Nwl; %custom NeuroDOT field
            end
            if ~isfield(snf.nirs.metaDataTags, 'comment')
            else
                info.io.comment = snf.nirs.metaDataTags.comment; %custom NeuroDOT field
            end
            if ~isfield(snf.nirs.metaDataTags, 'MeasurementDate')
%                snf.nirs.metaDataTags.MeasurementDate = 'n/a'; %required snirf field
            else
                info.io.date = join(convertCharsToStrings(snf.nirs.metaDataTags.MeasurementDate)); %required snirf field
            end
            if ~isfield(snf.nirs.metaDataTags, 'enc') 
            else
                info.io.enc = snf.nirs.metaDataTags.enc; %custom NeuroDOT field
            end
            if ~isfield(snf.nirs.metaDataTags, 'framesize')
            else
                info.io.framesize = snf.nirs.metaDataTags.framesize; %custom NeuroDOT field
            end
            if ~isfield(snf.nirs.metaDataTags, 'naux')
            else
                info.io.naux = snf.nirs.metaDataTags.naux; %custom NeuroDOT field 
            end
            if ~isfield(snf.nirs.metaDataTags, 'nblank')
            else
                info.io.nblank = snf.nirs.metaDataTags.nblank; %custom NeuroDOT field
            end
            if ~isfield(snf.nirs.metaDataTags, 'nframe')
            else
                info.io.nframe = snf.nirs.metaDataTags.nframe; %custom NeuroDOT field
            end
            if ~isfield(snf.nirs.metaDataTags, 'nmotu')
            else
                info.io.nmotu = snf.nirs.metaDataTags.nmotu; %custom NeuroDOT field
            end
            if ~isfield(snf.nirs.metaDataTags, 'nsamp')
            else
                info.io.nsamp = snf.nirs.metaDataTags.nsamp; %custom NeuroDOT field
            end
            if ~isfield(snf.nirs.metaDataTags, 'nts')
            else
                info.io.nts = snf.nirs.metaDataTags.nts; %custom NeuroDOT field
            end
            if ~isfield(snf.nirs.metaDataTags, 'pad')
            else
                info.io.pad = snf.nirs.metaDataTags.pad; %custom NeuroDOT field
            end
            if ~isfield(snf.nirs.metaDataTags, 'run')
            else
                info.io.run = snf.nirs.metaDataTags.run; %custom NeuroDOT field
            end
            if ~isfield(snf.nirs.metaDataTags, 'MeasurementTime')
            else
                info.io.time = join(convertCharsToStrings(snf.nirs.metaDataTags.MeasurementTime)); %required Snirf field
            end
            
            if ~isfield(snf.nirs.metaDataTags, 'tag')
            else
                info.io.tag = snf.nirs.metaDataTags.tag; %custom NeuroDOT field
            end
            
            if ~isfield(snf.nirs.metaDataTags, 'UnixTime')
            else
                info.io.unix_time = snf.nirs.metaDataTags.UnixTime; %optional Snirf field
            end
            if ~isfield(snf.nirs.metaDataTags, 'Nt')
            else
                info.io.Nt = snf.nirs.metaDataTags.Nt; %custom NeuroDOT field
            end
        end
        
        if ~isfield(snf.nirs,'probe')
        else
            if ~isfield(snf.nirs.metaDataTags,'CapName')
            else
                info.optodes.CapName = snf.nirs.metaDataTags.CapName; %custom NeuroDOT field
            end
            if ~isfield(snf.nirs.probe,'detectorPos2D')
            else
                info.optodes.dpos2  = snf.nirs.probe.detectorPos2D; %required Snirf/Ndot field
            end
            if ~isfield(snf.nirs.probe,'detectorPos3D')
            else
                info.optodes.dpos3 = snf.nirs.probe.detectorPos3D;%required Snirf/Ndot field
            end
            if ~isfield(snf.nirs.probe,'sourcePos2D')
            else
                info.optodes.spos2 = snf.nirs.probe.sourcePos2D; %required Snirf/Ndot field
            end
            if ~isfield(snf.nirs.probe,'sourcePos3D')
            else
                info.optodes.spos3 = snf.nirs.probe.sourcePos3D;%required Snirf/Ndot field
            end
            
        end
        
        if ~isfield(snf.nirs.data,'measurementList')
        else
            if ~isfield(snf.nirs.data.measurementList, 'sourceIndex')
            else
                info.pairs.Src = [snf.nirs.data.measurementList(:).sourceIndex]'; %required Snirf/Ndot field
            end
            if ~isfield(snf.nirs.data.measurementList, 'detectorIndex')
            else 
                info.pairs.Det = [snf.nirs.data.measurementList(:).detectorIndex]'; %required Snirf/Ndot field
            end            
            if ~isfield(snf.nirs.data.measurementList, 'wavelengthIndex')
            else  
                info.pairs.WL = [snf.nirs.data.measurementList(:).wavelengthIndex]'; %required Snirf/Ndot field
            end
            if ~isfield(snf.nirs.probe, 'wavelengths')
            else 
                wavelengths = snf.nirs.probe.wavelengths;
                info.pairs.lambda = info.pairs.WL;
                for j = 1:length(wavelengths)
                    info.pairs.lambda(info.pairs.WL == j) = wavelengths(j);
                end
            end
            if isfield(snf.nirs.data.measurementList, 'NN')
                info.pairs.NN = [snf.nirs.data.measurementList(:).NN]'; %custom Ndot field (in snirf format)
            else               
                gridTemp.spos3=snf.nirs.probe.sourcePos3D;
                gridTemp.dpos3=snf.nirs.probe.detectorPos3D;
                d = pdist2(gridTemp.spos3, gridTemp.dpos3);
                dmax = max(d(:));
                if dmax < 30
                    gridTemp.spos3=snf.nirs.probe.sourcePos3D.*10;
                    gridTemp.dpos3=snf.nirs.probe.detectorPos3D.*10;
                end
                Rad=Grid2Radius_180824(gridTemp,5);
                params.lambda= snf.nirs.probe.wavelengths;
                tempInfo=Generate_Info_from_Grid_Rad(gridTemp,Rad,params);
                [IdxA,IdxB]=ismember([tempInfo.pairs.Src,tempInfo.pairs.Det,tempInfo.pairs.WL],...
                    [info.pairs.Src,info.pairs.Det,info.pairs.WL],'rows');
                info.pairs.r3d=tempInfo.pairs.r3d(IdxA);
                info.pairs.r2d = tempInfo.pairs.r2d(IdxA);
                info.pairs.NN = tempInfo.pairs.NN(IdxA);
         
            end
                           
            if ~isfield(info.pairs, 'Mod')
            else 
                if isfield(snf.nirs.data.measurementList, 'Mod')
                    info.pairs.Mod = snf.nirs.data.measurementList(:).Mod'; %custom Ndot field (in snirf format)
                end
            end
            if ~isfield(info.pairs, 'r2d')
            else 
                if isfield(snf.nirs.data.measurementList, 'r2d')
                    info.pairs.r2d = snf.nirs.data.measurementList(:).r2d'; %custom Ndot field (in snirf format)
                end
            end
            if ~isfield(info.pairs, 'r3d')
            else 
                if isfield(snf.nirs.data.measurementList, 'r3d')
                    info.pairs.r3d = snf.nirs.measurementList(:).r3d'; %custom Ndot field (in snirf format)
                end
            end
        end
      
        if isfield(snf.nirs, 'stim')
            Total_synchs = [];
            Total_synchtypes = [];
            Npulses = length(snf.nirs.stim);
            for i = 1:length(snf.nirs.stim)
                Total_synchs = [Total_synchs; snf.nirs.stim(i).data(:,1)];
                Total_synchtypes = [Total_synchtypes; i.*ones(length(snf.nirs.stim(i).data(:,1)),1 )];
            end
            [info.paradigm.synchpts, sortedIdx] = sort(Total_synchs);
            info.paradigm.synchtype = Total_synchtypes(sortedIdx);
            for j = 1:Npulses
                info.paradigm.(['Pulse_', num2str(j)]) = find(info.paradigm.synchtype == j);            
            end
            info.paradigm.synchtimes = info.paradigm.synchpts;
            for j = 1: length(info.paradigm.synchpts)
                [~,info.paradigm.synchpts(j)] = min(abs(snf.nirs.data.time - info.paradigm.synchtimes(j)));
            end
        
        else
            fields = fieldnames(snf.nirs); 
            fieldlist = [];
        for j = 1:length(fields)
            if strlength(fields(j)) >= 4
                fieldlist = [fieldlist, fields(j)];
            else
                fieldlist = [fieldlist, 'not_stim'];
            end
        end
            idxPulse = ismember(cellfun(@(x) x(1:4), fieldlist, 'UniformOutput', false), 'stim');
            pulses = fields(idxPulse);
        if any(idxPulse) 
            Total_synchs = [];
            Total_synchtypes = [];
            Npulses = length(pulses);
            pulses = sort(pulses);
            for j = 1:Npulses
                pulse = string(pulses(j));
                Total_synchs = [Total_synchs; snf.nirs.(pulse).data(:,1)];
                Total_synchtypes = [Total_synchtypes; j.*ones(length(snf.nirs.(pulse).data(:,1)),1 )];
            end
            [info.paradigm.synchpts, sortedIdx] = sort(Total_synchs);
            info.paradigm.synchtype = Total_synchtypes(sortedIdx);
            for j = 1:Npulses
                info.paradigm.(['Pulse_', num2str(j)]) = find(info.paradigm.synchtype == j);            
            end
            info.paradigm.synchtimes = info.paradigm.synchpts;
            for j = 1: length(info.paradigm.synchpts)
                [~,info.paradigm.synchpts(j)] = min(abs(snf.nirs.data.time - info.paradigm.synchtimes(j)));
            end      
        end
        end
    end
    
        
        if ~isfield(snf.nirs.metaDataTags, 'SubjectID')
            info.misc.subject_id = 'default'; %required snirf field
        else
            info.misc.subject_id = join(convertCharsToStrings(snf.nirs.metaDataTags.SubjectID));
        end
    
% Save Output NeuroDOT File
[p,f,e]=fileparts(output);
outputfilename=fullfile(p,f);

save(outputfilename,'data','info');

end