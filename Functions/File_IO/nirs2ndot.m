function [data, info] = nirs2ndot(filename, save_file, output)

    % Translation from .nirs to NeuroDOT-compatible format
    % This function takes in a .nirs file and converts it to NeuroDOT 
    % compatible variables: data and info
        % The input 'filename' should contain the full file name including
        % the *.nirs extension.
        % data is the nirs data in the format of: N_meas x N_samples
        % info is the data structure that holds information pertaining to data aquisition
    % The flag 'save_file' can be set to zero to suppress saving out the 
    % 'data' and 'info' variables to a NeuroDOT compatible *.mat file.
        % By default, the *.mat file will be saved (save_file = 10
    % The optional input 'output' defines the name of the output *.mat file
        % By default, the output filename will match the input filename.
        
%
% Copyright (c) 2017 Washington University 
% Created By: Ari Segel
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
    if ~exist('output','var')
        output = filename;
    end
    
    if ~exist('save_file','var')
        save_file = 1;
    end
    
    
    %% Data
    nirsData = load(filename, '-mat');
    data = nirsData.d';
    
    
    %% IO and System parameters
    
    % IO
    info.io.Nd = nirsData.SD.nDets;
    info.io.Ns = nirsData.SD.nSrcs;
    info.io.Nwl = length(nirsData.SD.Lambda);
    info.io.nframe = length(nirsData.t);

    % Framerate
	info.system.framerate = 1/(mean(diff(nirsData.t)));
    if nirsData.t(1) == 0
        info.misc.startTime = 0;
    else
        info.misc.startTime = 1;
    end
    
    
    
    %% Paradigm

    % Initialize variables
    num_stim = size(nirsData.s,2);
    num_synchs = sum(nirsData.s == 1);
    field_names = cell(num_stim,1);
    
    if num_synchs == 0
        % do nothing, don't create and fill info.paradigm as it'll be empty anywyas
    else
        % Get synchs for each stim type
        synchs = struct([]);
        for j = 1:num_stim
            synchs(j).tp = find(nirsData.s(:,j) == 1);    
        end

        %combine synchs to 1 array and get synchpts
        synchTot = sort(cat(1, synchs.tp));
        info.paradigm.synchpts = synchTot;

        % Populate Pulse fields and synchtype
        % Note, synchtype will have 1's for rest, and then stimulus types will start at 2 and increase by 1 for each different stim type in the order stim was presented
        % For example, if you have a paradigm with left and right stimuli where left stim is presented first and right stim is presented second...
        % Synchtypes will be as follows: Rest = 1, Left = 2, Right = 3;
        info.paradigm.synchtype = zeros(size(info.paradigm.synchpts)); % initialize synchtype AFTER synchpts created so its the correct size
        for k = 1:num_stim
            field_names{k} = ['Pulse_', num2str(k)]; % get Pulse names
            [~,info.paradigm.(field_names{k})] = ismember(synchs(k).tp, synchTot); % set Pulse
            info.paradigm.synchtype(info.paradigm.(field_names{k})) = k; % set synchtype
        end

        % Set synchtimes
        info.paradigm.synchtimes = synchTot/info.system.framerate;
    end
    
    
    %% Optodes
    
    if (isfield(nirsData.SD, 'SrcPos') && isfield(nirsData.SD, 'DetPos')) %if non dimensional
        % check size
        dimension = size(nirsData.SD.SrcPos,2);
        
        if dimension == 3 %if 3D coords, get 3D pos, then check for 2D specific pos
            spos3 = nirsData.SD.SrcPos;
            dpos3 = nirsData.SD.DetPos;
            if (isfield(nirsData.SD, 'SrcPos2') && isfield(nirsData.SD, 'DetPos2'))
                spos2 = nirsData.SD.SrcPos2;
                dpos2 = nirsData.SD.DetPos2;
            else %if no 2D pos, tell the user they need to create their own
                disp(['No 2D optode positions in NIRS data, please create your own 2D layout.',...
                    newline, ' Place the 2D source coordinates in info.optodes.spos2',... 
                    ' and the 2D detector coordinates in info.optodes.dpos2.'])
            end
        elseif dimension == 2 % if 2D coords, check for 3D specific pos first
            if (isfield(nirsData.SD, 'SrcPos3') && isfield(nirsData.SD, 'DetPos3')) %if 3D coords, get them, then get 2D coords
                spos3 = nirsData.SD.SrcPos3;
                dpos3 = nirsData.SD.DetPos3;
                spos2 = nirsData.SD.SrcPos;
                dpos2 = nirsData.SD.DetPos;
            else
                disp(['No 3D optode coordinates in this nirs data.',...
                    newline, '    Please find the 3D SD coordinates and place the source coordinates in',...
                    ' info.optodes.spos3', newline, '    and the detector coordinates in info.optodes.dpos3.',...
                    newline, '    You may need to make a neuroDOT pad file and/or head model to do this.'])
            end
        else 
            disp('Why do the optode arrays have something other than 2 or 3 columns?????')
            disp(['No optode coordinates in this nirs data.',...
                newline, '    Please find the 3D SD coordinates and place the source coordinates in',...
                ' info.optodes.spos3', newline, '    and the detector coordinates in info.optodes.dpos3.',...
                newline, '    You may need to make a neuroDOT pad file and/or head model to do this.'])
            disp(['No 2D optode positions in NIRS data, please create your own 2D layout.',...
                newline, '    Place the 2D source coordinates in info.optodes.spos2',... 
                ' and the 2D detector coordinates in info.optodes.dpos2.'])
        end
    
    elseif (isfield(nirsData.SD, 'SrcPos3') && isfield(nirsData.SD, 'DetPos3')) %if dimension-specific
        spos3 = nirsData.SD.SrcPos3;
        dpos3 = nirsData.SD.DetPos3;
        if (isfield(nirsData.SD, 'SrcPos2') && isfield(nirsData.SD, 'DetPos2')) %check for 2D pos
            spos2 = nirsData.SD.SrcPos2;
            dpos2 = nirsData.SD.DetPos2;
        else %if no 2D pos, tell the user they need to create their own
            disp(['No 2D optode positions in NIRS data, please create your own 2D layout.',...
                newline, '    Place the 2D source coordinates in info.optodes.spos2',... 
                ' and the 2D detector coordinates in info.optodes.dpos2.'])
        end
    else 
        disp(['No optode coordinates in this nirs data.',...
                newline, '    Please find the 3D SD coordinates and place the source coordinates in',...
                ' info.optodes.spos3', newline, '    and the detector coordinates in info.optodes.dpos3.',...
                newline, '    You may need to make a neuroDOT pad file and/or head model to do this.'])
        disp(['No 2D optode positions in NIRS data, please create your own 2D layout.',...
                newline, '    Place the 2D source coordinates in info.optodes.spos2',... 
                ' and the 2D detector coordinates in info.optodes.dpos2.'])

    end
    
    % Variables for generating info.optodes and info.pairs
    lambda = nirsData.SD.Lambda;
    SD_sep = [];
    r2dArray = [];
    r3dArray = [];
    
    % Calculate SD separations
    for ii = 1:size(nirsData.SD.MeasList,1)
        SD_sep = [SD_sep; pdist2(spos3(nirsData.SD.MeasList(ii,1),:),...
            dpos3(nirsData.SD.MeasList(ii,2),:))]; %changed from SD3D to SD 12/8/22
    end
    
    % Check if SD units are in mm 
    avg_SD_sep = mean(abs(SD_sep(:))); %average SD separation
    min_SD_sep = min(abs(SD_sep(:))); %minimum SD separation
    
    % if avg between 10 and 100 units OR min is within an order of magnitude below range, coords are already in mm, so we're good
    if ((avg_SD_sep >= 10) && (avg_SD_sep < 100)) || ((min_SD_sep >= 1) && (min_SD_sep < 100))  
        mult = 1;
        
    % if between 1 and 100 OR min is within an order of magnitude below range, units are in cm  
    elseif (avg_SD_sep >= 1) && (avg_SD_sep < 10) || ((min_SD_sep >= 0.1) && (min_SD_sep < 10)) 
        mult = 10;
        
    % if between 0.1 and 1 OR min is within an order of magnitude below range, units are in dm   
    elseif (avg_SD_sep >= 0.1) && (avg_SD_sep < 1) || ((min_SD_sep >= 0.01) && (min_SD_sep < 1))
        mult = 100;
        
    % if between 0.01 and 0.1 OR min is within an order of magnitude below range, units are in m    
    elseif (avg_SD_sep >= 0.01) && (avg_SD_sep < 0.1) || ((min_SD_sep >= 0.001) && (min_SD_sep < 0.1))
        mult = 1000;
    
    % if avg outside min/max of full range (0.01, 100) AND min is not in range either, display warning
    elseif (avg_SD_sep < 0.01) || (min_SD_sep < 0.001)
        disp('Optode position units are larger than meters, data is wonky, please fix your data')
    elseif (avg_SD_sep > 100) || (min_SD_sep < 20) %using 20 here because I've seen one cap have units in mm, but 20 > min SD sep > 10
        disp('Optode position units are smaller than mm, data is wonky, please fix your data')
    end
    
    %Adjust optodes pos to be mm
    spos3 = spos3*mult;
    dpos3 = dpos3*mult;
    
    %Enforce column-wise optode arrays
    if size(spos3, 2) > size(spos3, 1)
        spos3 = spos3';
        dpos3 = dpos3';
        if exist('spos2', 'var')
            spos2 = spos2';
            dpos2 = dpos2';
        end
    end

    % Place optode locations and plot orientation in info structure
    info.optodes.spos3 = spos3;
    info.optodes.dpos3 = dpos3;
    info.optodes.plot3orientation.i='R2L';  
    info.optodes.plot3orientation.j='P2A'; 
    info.optodes.plot3orientation.k='D2V'; 
    
    
    %% Pairs
    
    % Calculate SD separations: r2d and r3d
    % Also make vectors that contain wavelengths
    for ii = 1:size(nirsData.SD.MeasList,1)
        r3dArray = [r3dArray; pdist2(spos3(nirsData.SD.MeasList(ii,1),:),...
            dpos3(nirsData.SD.MeasList(ii,2),:))]; 
        if nirsData.SD.MeasList(ii,4) == 1 
            lambdaArray(ii) = lambda(1);
        else
            lambdaArray(ii) = lambda(2);
        end
    end
    
    % Fill info.pairs fields
    info.pairs.Src = nirsData.SD.MeasList(:,1);
    info.pairs.Det = nirsData.SD.MeasList(:,2);
    info.pairs.NN = nirsData.SD.MeasList(:,3); %initialize as all 1's, we'll calculate NN below
    info.pairs.WL = nirsData.SD.MeasList(:,4);
    info.pairs.Mod = repmat(['CW'], size(nirsData.SD.MeasList,1),1);
    info.pairs.r3d = r3dArray;
    
    % Calculate NNs
    info = calc_NN(info,10); %minimum separation for binning NN = 10mm

    % Do lambda
    info.pairs.lambda = lambdaArray'; %transpose array so it is column vector
    
    
    %% Tissue
    
    % This structure contains information about the light model you will generate. 
    info.tissue.affine=eye(4); % if affine transform to atlas space is Identity
    info.tissue.affine_target='MNI'; % Atlas target from your model
    
    %% Do 2D 
    
    if exist('spos2', 'var') %Check if there are 2D pos, this variable will only be present if 2D pos are available from nirs data 
        % Correct for units (distances must be in mm)
        dpos2 = dpos2*mult;
        spos2 = spos2*mult;
        
        % Place optode pos in structure
        info.optodes.dpos2 = dpos2;
        info.optodes.spos2 = spos2;
        
        % Calculate 2D SD separations
        for ii = 1:size(nirsData.SD.MeasList,1)
            r2dArray = [r2dArray; pdist2(spos2(nirsData.SD.MeasList(ii,1),:),...
                dpos2(nirsData.SD.MeasList(ii,2),:))];
        end
        info.pairs.r2d = r2dArray;
    else % If no 2D coordinates are present in the .nirs file
        info.pairs.r2d = info.pairs.r3d;
    end
    
    % Save Output NeuroDOT File
    if save_file == 1
        [p,f,e]=fileparts(output);
        outputfilename=fullfile(p,f);

        save(outputfilename,'data','info');
    end
end
