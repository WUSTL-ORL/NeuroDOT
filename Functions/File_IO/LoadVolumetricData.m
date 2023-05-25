function [volume, header] = LoadVolumetricData(filename, pn, file_type)

% LOADVOLUMETRICDATA Loads a volumetric data file.
%
%   [volume, header] = LOADVOLUMETRICDATA(filename, pn, file_type) loads a
%   file specified by "filename", path "pn", and "file_type", and returns
%   it in two parts: the raw data file "volume", and the header "header",
%   containing a number of key-value pairs in 4dfp format.
%
%   [volume, header] = LOADVOLUMETRICDATA(filename) supports a full
%   filename input, as long as the extension is included in the file name
%   and matches a supported file type.
% 
%   Supported File Types/Extensions: '.4dfp' 4dfp, 'nii' NIFTI.
% 
%   NOTE: This function uses the NIFTI_Reader toolbox available on MATLAB
%   Central. This toolbox has been included with NeuroDOT 2.
% 
% Dependencies: READ_4DFP_HEADER, READ_NIFTI_HEADER, MAKE_NATIVESPACE_4DFP.
% 
% See Also: SAVEVOLUMETRICDATA.
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
volume = [];
header = [];

if ~exist('file_type', 'var')  &&  ~exist('pn', 'var')
    [pn, filename, file_type] = fileparts(filename);
    file_type = file_type(2:end);
end


switch lower(file_type)
    case '4dfp'
        %% Read 4dfp file.
        % Read .ifh header.
        header = Read_4dfp_Header(fullfile(pn, [filename, '.', file_type, '.ifh']));
        
        % Read .img file.
        fid = fopen(fullfile(pn, [filename, '.', file_type, '.img']), 'r', header.byte);
        volume = fread(fid, header.nVx * header.nVy * header.nVz * header.nVt, header.format);
        fclose(fid);
        
        %% Put header into native space if not already.
        header = Make_NativeSpace_4dfp(header);
        
        %% Format for output.
        volume = squeeze(reshape(volume, header.nVx, header.nVy, header.nVz, header.nVt));
        
        switch header.acq
            case 'transverse'
                volume = flip(volume, 2);
            case 'coronal'
                volume = flip(volume, 2);
                volume = flip(volume, 3);
            case 'sagittal'
                volume = flip(volume, 1);
                volume = flip(volume, 2);
                volume = flip(volume, 3);  
        end
case {'nifti', 'nii' ,'nii.gz', '.nii'}   
        %% Call NIFTI_Reader function.
        %%% NOTE: When passing file types, if you have the ".nii" file
        %%% extension, you must use that as both the "ext" input AND add it
        %%% as an extension on the "filename" input.
            if strcmp(file_type, '.nii')
                file_type = 'nii';
            end
            if strcmp(file_type, 'nii.gz')
               nii.img = niftiread(fullfile(pn,[filename, '.nii.gz']));
               nii.hdr = niftiinfo(fullfile(pn,[filename, '.nii.gz']));
            else
                nii.img = niftiread(fullfile(pn,[filename, '.', file_type]));
                nii.hdr = niftiinfo(fullfile(pn,[filename, '.', file_type]));
            end
            % Implemented 2/20/2023 ES
            header = nifti_4dfp(nii.hdr, '4'); % Convert nifti format header to 4dfp/NeuroDOT format
            
            if isfield(nii, 'img')
                volume = flip(nii.img, 1); % NIFTI loads in RAS orientation. We want LAS, so we flip first dim.
                header.format = class(nii.img);
                header.original_header=rmfield(nii,'img');
            end
            
            header.version_of_keys = '3.3'; 
            header.conversion_program = 'NeuroDOT_LoadVolumetricData';
            header.filename = [filename,'.nii'];
            header.nVx = header.matrix_size(1);
            header.nVy = header.matrix_size(2);
            header.nVz = header.matrix_size(3);
            header.nVt = header.matrix_size(4);
            header.mmx = abs(header.mmppix(1));
            header.mmy = abs(header.mmppix(2));
            header.mmz = abs(header.mmppix(3));
            
            orientation = header.orientation;
            switch orientation
                case '2'
                    header.acq = 'transverse';
                case '3'
                    header.acq = 'coronal';
                case '4'
                    header.acq = 'sagittal';
            end
            if isfield(nii, 'machine')
                switch nii.machine
                    case 'ieee-le'
                        header.byte = 'l';
                    case 'ieee-be'
                        header.byte = 'b';
                end 
            end
           
           
             
if ~isempty(volume)
volume=double(volume);end
end
end
