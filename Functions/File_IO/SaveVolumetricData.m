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
        
    case {'nifti', 'nii','.nii'}
        % Implemented 2/20/2023 ES
        if strcmp(file_type, '.nii')
            file_type = 'nii';
        end
        %% Call NIFTI_Reader functions.
        volume = flip(volume, 1); % Convert back from LAS to RAS for NIFTI.
        
        if isfield(header,'original_header') % was loaded as nii
            nii = struct;
            nii.img=volume;
            nii.hdr = nifti_4dfp(header, 'n');
            nii.hdr.original_header=header.original_header.hdr;  
             % Required fields
            nii.hdr.Description = nii.hdr.descrip;
            if length(nii.hdr.Description) > 80
                nii.hdr.Description = 'Converted with NeuroDOT nifti_4dfp';
            end
            nii.hdr.ImageSize = nii.hdr.dim(2:4);
            if nii.hdr.datatype == 16
                nii.hdr.Datatype = 'single';
            end
            nii.hdr.BitsPerPixel = 32;
            nii.hdr.Version = 'NIfTI1';
            nii.hdr.Qfactor = nii.hdr.pixdim(1);
            nii.hdr.PixelDimensions = nii.hdr.pixdim(2:4);
            nii.hdr.SpaceUnits = 'Millimeter';
            nii.hdr.TimeUnits = 'Second';
            nii.hdr.AdditiveOffset = 0;
            nii.hdr.MultiplicativeScaling = 0;
            nii.hdr.TimeOffset = 0;
            nii.hdr.SliceCode = 'Unknown';
            nii.hdr.FrequencyDimension = 0;
            nii.hdr.PhaseDimension = 0;
            nii.hdr.SpatialDimension = 0;
            nii.hdr.DisplayIntensityRange = [0,0];
            nii.hdr.TransformName = 'Sform';
            nii.hdr.Transform.T = [nii.hdr.srow_x(1), ...
                nii.hdr.srow_y(1),nii.hdr.srow_z(1),0;...
                nii.hdr.srow_x(2), nii.hdr.srow_y(2),...
                nii.hdr.srow_z(2),0;nii.hdr.srow_x(3),...
                nii.hdr.srow_y(3),nii.hdr.srow_z(3),0;...
                nii.hdr.srow_x(4), nii.hdr.srow_y(4),...
                nii.hdr.srow_z(4),1];
            nii.hdr.Transform.Dimensionality = 4;
        else % Build nii header from information in NeuroDOT header
            nii = struct;
            nii.img = volume;
            % Required fields
            nii.hdr.raw = nifti_4dfp(header, 'n');
            nii.hdr.Description = nii.hdr.raw.descrip;
            nii.hdr.ImageSize = nii.hdr.raw.dim(2:4);
            if nii.hdr.raw.datatype == 16
                nii.hdr.Datatype = 'single';
            end
            nii.hdr.BitsPerPixel = 32;
            nii.hdr.Version = 'NIfTI1';
            nii.hdr.Qfactor = nii.hdr.raw.pixdim(1);
            nii.hdr.PixelDimensions = nii.hdr.raw.pixdim(2:4);
            nii.hdr.SpaceUnits = 'Millimeter';
            nii.hdr.TimeUnits = 'Second';
            nii.hdr.AdditiveOffset = 0;
            nii.hdr.MultiplicativeScaling = 0;
            nii.hdr.TimeOffset = 0;
            nii.hdr.SliceCode = 'Unknown';
            nii.hdr.FrequencyDimension = 0;
            nii.hdr.PhaseDimension = 0;
            nii.hdr.SpatialDimension = 0;
            nii.hdr.DisplayIntensityRange = [0,0];
            nii.hdr.TransformName = 'Sform';
            nii.hdr.Transform.T = [nii.hdr.raw.srow_x(1), ...
                nii.hdr.raw.srow_y(1),nii.hdr.raw.srow_z(1),0;...
                nii.hdr.raw.srow_x(2), nii.hdr.raw.srow_y(2),...
                nii.hdr.raw.srow_z(2),0;nii.hdr.raw.srow_x(3),...
                nii.hdr.raw.srow_y(3),nii.hdr.raw.srow_z(3),0;...
                nii.hdr.raw.srow_x(4), nii.hdr.raw.srow_y(4),...
                nii.hdr.raw.srow_z(4),1];
            nii.hdr.Transform.Dimensionality = 4;
        end
        % In order to write the header, you need a very specific set of
        % fields
        % original_header does not get saved due to it not being one of the
        % fields that niftiwrite considers
%         save_nii(nii, filename);
        niftiwrite(single(nii.img), filename, nii.hdr);

end
end


%
