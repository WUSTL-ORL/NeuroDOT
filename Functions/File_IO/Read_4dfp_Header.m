function header = Read_4dfp_Header(filename)

% READ_4DFP_HEADER Reads the .ifh header of a 4dfp file.
%
%   header = READ_4DFP_HEADER(filename) reads an .ifh text file specified
%   by "filename", containing a number of key-value pairs. The specific
%   pairs are parsed and stored as fields of the output structure "header".
%
% See Also: LOADVOLUMETRICDATA.
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

%% Parameters and initialization.
header = [];

%% Open file.
fid = fopen(filename, 'r', 'b'); % Must be big-endian format.

%% Read text.
temp = deblank(textscan(fid, '%s %s', 'Delimiter', ':='));

%% Prepare key-value pairs.
temp = temp{1}; % We only need the first set of strings; the second are just empties.

if rem(size(temp, 1), 2) ~= 0 % Removing the first line.
    temp(1) = [];
end

temp = reshape(temp, [2, size(temp, 1)/2]); % Reshape so that these can be fed into a for loop.

%% Create list of fields to load into.

%% Interpret pairs.
for pairs = temp
    
    key = pairs{1};
    value = pairs{2};
    
    switch key
        case 'version of keys'
            header.version_of_keys = value;
        case 'number format'
            header.format = value;
        case 'conversion program'
            header.conversion_program = value;
        case 'name of data file'
            header.filename = value;
        case 'number of bytes per pixel'
            header.bytes_per_pixel = str2double(value);
        case 'imagedata byte order'
            switch value
                case 'bigendian'
                    header.byte = 'b';
                case 'littleendian'
                    header.byte = 'l';
            end
        case 'orientation'
            switch value
                case '2'
                    header.acq = 'transverse';
                case '3'
                    header.acq = 'coronal';
                case '4'
                    header.acq = 'sagittal';
            end
        case 'number of dimensions'
            header.nDim = str2double(value);
        case 'matrix size [1]'
            header.nVx = str2double(value);
        case 'matrix size [2]'
            header.nVy = str2double(value);
        case 'matrix size [3]'
            header.nVz = str2double(value);
        case 'matrix size [4]'
            header.nVt = str2double(value);
        case 'scaling factor (mm/pixel) [1]'
            header.mmx = str2double(value);
        case 'scaling factor (mm/pixel) [2]'
            header.mmy = str2double(value);
        case 'scaling factor (mm/pixel) [3]'
            header.mmz = str2double(value);
        case 'patient ID'
            header.subjcode = value;
        case 'date'
            header.filedate = value;
        case 'mmppix'
            header.mmppix = str2num(value); %#ok<*ST2NM>
        case 'center'
            header.center = str2num(value);
    end
end
if ~isfield(header,'byte'),header.byte = 'b';end

%% Close file.
fclose(fid);



%
