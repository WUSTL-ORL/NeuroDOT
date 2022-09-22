function SaveVolumetricData(volume, header, filename, pn, file_type)

% SAVEVOLUMETRICDATA Saves a volumetric data file.
%
%   SAVEVOLUMETRICDATA(volume, header, filename, pn, file_type) saves
%   volumetric data defined by "volume" and "header" into a file specified
%   by "filename", path "pn", and "file_type".
%
%   SAVEVOLUMETRICDATA(volume, header, filename) supports a full
%   filename input, as long as the extension is included in the file name
%   and matches a supported file type.
%
%   Supported File Types/Extensions: '.4dfp' 4dfp, '.nii' NIFTI.
%
%   NOTE: This function uses the NIFTI_Reader toolbox available on MATLAB
%   Central. This toolbox has been included with NeuroDOT 2.
%
% Dependencies: WRITE_4DFP_HEADER, WRITE_NIFTI.
%
% See Also: LOADVOLUMETRICDATA, MAKE_NATIVESPACE_4DFP.
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
if ~exist('file_type', 'var')  &&  ~exist('pn', 'var')
    [pn, filename, file_type] = fileparts(filename);
    file_type = file_type(2:end);
end

switch lower(file_type)
    case '4dfp'
        if header.nVt~=size(volume,4)
        disp(['Warning: Stated header 4th dimension size ',...
        num2str(header.nVt),' does not equal the size of the volume ',...
        num2str(size(volume,4)),'. Updating header info.'])
        header.nVt=size(volume,4);
        end
        %% Write 4dfp header.
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
        
        Write_4dfp_Header(header, fullfile(pn, filename))
     
        %% Write 4dfp image file.
        volume = squeeze(reshape(volume, header.nVx * header.nVy * header.nVz * header.nVt, 1));
        fid = fopen([fullfile(pn, filename), '.4dfp.img'], 'w', 'l');
        fwrite(fid, volume, header.format);
        fclose(fid);
        
    case {'nifti', 'nii'}
        %% Call NIFTI_Reader functions.
        volume = flip(volume, 1); % Convert back from LAS to RAS for NIFTI.
        
        if isfield(header,'original_header') % was loaded as nii
            nii.img=volume;
            nii.hdr=header.original_header.hdr;
            nii.hdr.dime.dim=[length(size(volume)),size(volume,1),...
                                size(volume,2),size(volume,3),...
                                size(volume,4),1 1 1];
            
        else % Build nii header from information in NeuroDOT header
             % .nii headers require hdr.hk, hdr.dime, hdr.hist
            nii = struct;
            nii.img = volume;
            % Required fields
            nii.hdr.hk.sizeof_hdr = 348;  % Must be 348
            
            nii.hdr.hist.magic = 'n+1';   % 'n+1' if the image volume is in the same .nii file as the header
            nii.hdr.hist.descrip = header.filename;
            
            nii.hdr.dime.datatype = 16;   
            nii.hdr.dimebitpix = 32;      % Should correspond correctly to datatype 
            nii.hdr.dime.vox_offset = 352;% Required for an 'n+1' header
            nii.hdr.dime.dim(1) = header.nDim;
            nii.hdr.dime.dim(2) = header.nVx;
            nii.hdr.dime.dim(3) = header.nVy;
            nii.hdr.dime.dim(4) = header.nVz;
            nii.hdr.dime.dim(5) = header.nVt;
            nii.hdr.dime.dim(6) = 1;
            nii.hdr.dime.dim(7) = 1;
            nii.hdr.dime.dim(8) = 1;
            nii.hdr.dime.pixdim(1) = 1;
            nii.hdr.dime.pixdim(2) = header.mmx;
            nii.hdr.dime.pixdim(3) = header.mmy;
            nii.hdr.dime.pixdim(4) = header.mmz;
            nii.hdr.dime.pixdim(5) = 1;
            nii.hdr.dime.pixdim(6) = 0;
            nii.hdr.dime.pixdim(7) = 0;
            nii.hdr.dime.pixdim(8) = 0;
                    
            % Calculate Srow values from Center 
            nii.hdr.hist.srow_x = [header.mmx, 0, 0, -(header.center(1)- header.mmx)];
            nii.hdr.hist.srow_y = [0, header.mmy, 0, -(header.center(2)- (header.mmy*header.nVy*(-1)))];
            nii.hdr.hist.srow_z = [0, 0, header.mmz, -(header.center(3)- (header.mmz*header.nVz*(-1)))];
            
            % Optional fields must be set to zero if unused
            nii.hdr.hk.dim_info = 0;
            nii.hdr.hk.data_type = 0;
            nii.hdr.hk.db_name = 0;
            nii.hdr.hk.extents = 0;
            nii.hdr.hk.session_error = 0;
            nii.hdr.hk.regular = blanks(1);
            nii.hdr.dime.intent_p1 = 0;
            nii.hdr.dime.intent_p2 = 0;
            nii.hdr.dime.intent_p3 = 0;
            nii.hdr.dime.intent_code = 0;
            nii.hdr.dime.scl_slope = 0;
            nii.hdr.dime.scl_inter = 0;
            nii.hdr.dime.xyzt_units = 2; % Milimeters
            nii.hdr.dime.cal_max = 0;    
            nii.hdr.dime.cal_min = 0; 
            nii.hdr.dime.toffset = 0;
            nii.hdr.dime.slice_code = 0;
            nii.hdr.dime.slice_start = 0;
            nii.hdr.dime.slice_end = 0;
            nii.hdr.dime.slice_duration = 0;
            nii.hdr.hist.intent_name = '\0';
            nii.hdr.hist.qform_code = 1; % Not recommended to set to 0
            nii.hdr.hist.quatern_b = 0;
            nii.hdr.hist.quatern_c = 0;
            nii.hdr.hist.quatern_d = 0;
            nii.hdr.hist.qoffset_x = 0;
            nii.hdr.hist.qoffset_y = 0;
            nii.hdr.hist.qoffset_z = 0;
            nii.hdr.hist.sform_code = 1; % Not recommended to set to 0
            nii.hdr.hist.originator(1) = -nii.hdr.hist.srow_x(4) + 1;
            nii.hdr.hist.originator(2) = -nii.hdr.hist.srow_y(4) + 1;
            nii.hdr.hist.originator(3) = -nii.hdr.hist.srow_z(4) + 1;
            nii.hdr.hist.aux_file = '';
                           
        save_nii(nii, [fullfile(pn, filename), '.nii']);
        end
end
end



%
