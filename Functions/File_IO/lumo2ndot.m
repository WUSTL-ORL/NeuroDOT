function [data, info, clipping] = lumo2ndot(filename, save_file, output)
%
% lumo2ndot(filename,output,) takes a file with the '.lufr' extension
% in the LumoData format and converts it to NeuroDOT formatting. 
% 'Filename' is the name of the file to be converted, followed by the 
% .lufr extension.

% 'Save_file' is a flag which determines whether the data will be saved to
% a '.mat' file in NeuroDOT format. File will be saved be default if this
% is not set.
%
% 'Output' is the filename (without extension) to be saved in NeuroDOT
% format. Output files are saved as '.mat.' 

% Convert .lufr data format to be used in NeuroDOT
% Only the fields that are needed to run NeuroDOT are extracted

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

%% I/O
    
if ~exist('save_file','var')
    save_file = 1;
end

if ~exist('output','var')
    output = filename;
end
    
%% Load and get data/info from new LUMO data architecture

lumoData = LumoData(filename);
snirfData = lumoData.write_SNIRF([output, '.snf']); %snirf is easier to deal with than LumoData structure
events = lumoData.evts;
clipping = any(lumoData.data.chn_sat,2);
accelerometer = lumoData.data.node_acc;
gyroscope = lumoData.data.node_gyr;
if isempty(events) %handle what happens when the events structure in the LUMO data is empty
    disp([newline, 'Events structure is empty, event recording may have been turned off on the system.'])
    disp(['There will be no info.paradigm structure in the output: "info", you may need to make this yourself.', newline])
end
    
%% Run SNIRF converter
[data, info] = snirf2ndot([],1,output,snirfData);

% make sure pairs structure is doubles and not int32's
info.pairs.Src = double(info.pairs.Src);
info.pairs.Det = double(info.pairs.Det);
info.pairs.WL = double(info.pairs.WL);
info.pairs.lambda = double(info.pairs.lambda);
info.pairs.r3d = double(info.pairs.r3d);
info.pairs.r2d = double(info.pairs.r2d);
info.pairs.NN = double(info.pairs.NN);

% Add auxiliary data to info.misc
info.misc.accel = accelerometer;
info.misc.gyr = gyroscope;

%% Save Output NeuroDOT File

if save_file == 1
    [p,f,e]=fileparts(output);
    outputfilename=fullfile(p,f);

    save(outputfilename,'data','info', 'clipping');
end
    
end
