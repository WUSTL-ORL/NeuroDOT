function Write_4dfp_Header(header, filename)

% WRITE_4DFP_HEADER Writes 4dfp header to .ifh file.
% 
%   WRITE_4DFP_HEADER(header, filename) writes the input "header" in 4dfp
%   format to an .ifh file specified by "filename".
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
fid = fopen([filename, '.4dfp.ifh'], 'wt', 'l');

%% Print input header to file.
fprintf(fid, 'INTERFILE :=\n');
fprintf(fid, ['version of keys := ', header.version_of_keys, '\n']);
fprintf(fid, 'image modality := dot\n');
fprintf(fid, 'originating system := Neuro-DOT\n');
fprintf(fid, 'conversion program := MATLABto4dfp\n');
fprintf(fid, 'original institution := Washington University\n');
fprintf(fid, ['number format := ', header.format, '\n']);
fprintf(fid, ['name of data file := ', sprintf(filename), '\n']);
fprintf(fid, ['number of bytes per pixel := ', num2str(header.bytes_per_pixel), '\n']);

switch header.byte
    case 'b'
        byte='big';
    case 'l'
        byte='little';
end
fprintf(fid, ['imagedata byte order := ' byte 'endian\n']);

switch header.acq
    case 'transverse'
        fprintf(fid, 'orientation := 2\n');
    case 'coronal'
        fprintf(fid, 'orientation := 3\n');
    case 'sagittal'
        fprintf(fid, 'orientation := 4\n');
end

fprintf(fid, ['number of dimensions := ', num2str(header.nDim), '\n']);
fprintf(fid, ['matrix size [1] := ', num2str(header.nVx), '\n']);
fprintf(fid, ['matrix size [2] := ', num2str(header.nVy), '\n']);
fprintf(fid, ['matrix size [3] := ', num2str(header.nVz), '\n']);
fprintf(fid, ['matrix size [4] := ', num2str(header.nVt), '\n']);
fprintf(fid, ['scaling factor (mm/pixel) [1] := ', num2str(header.mmx), '\n']);
fprintf(fid, ['scaling factor (mm/pixel) [2] := ', num2str(header.mmy), '\n']);
fprintf(fid, ['scaling factor (mm/pixel) [3] := ', num2str(header.mmz), '\n']);
if isfield(header, 'mmppix')
    fprintf(fid, ['mmppix := ', num2str(header.mmppix), '\n']);
end
if isfield(header, 'center')
    fprintf(fid, ['center := ', num2str(header.center), '\n']);
end

fclose(fid);



%
