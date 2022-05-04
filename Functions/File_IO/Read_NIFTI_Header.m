function header = Read_NIFTI_Header(nii)

% READ_NIFTI_HEADER Reads header of "nii" structure into 4dfp.
%
%   header = READ_NIFTI_HEADER(nii) reads the structure "nii" and converts
%   it to a 4dfp format header "header".
% 
%   NOTE: The math in this function to calculate the "center" field for
%   4dfp is not straightforward, and is based on source code shared by Tim
%   Coalson from his "nifti-4dfp" toolbox.
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

%% Parameters and Initialization
header = [];

%% Fill in missing fields.
header.version_of_keys = '3.3'; % Not entirely sure here, but it's the
% version that this math came from.
header.conversion_program = 'NeuroDOT_2_LoadVolumetricData';

%% Read out elements.
[~, header.filename, ~] = fileparts(nii.fileprefix);
header.filename = [header.filename '.nii'];
header.format = class(nii.img);
header.bytes_per_pixel = nii.hdr.dime.bitpix / 8;

switch nii.machine
    case 'ieee-le'
        header.byte = 'l';
    case 'ieee-be'
        header.byte = 'b';
end

header.acq = 'transverse';

header.nDim = nii.hdr.dime.dim(1);
header.nVx = nii.hdr.dime.dim(2);
header.nVy = nii.hdr.dime.dim(3);
header.nVz = nii.hdr.dime.dim(4);
header.nVt = nii.hdr.dime.dim(5);

header.mmx = nii.hdr.dime.pixdim(2);% updated from 1:3 191117
header.mmy = nii.hdr.dime.pixdim(3);
header.mmz = nii.hdr.dime.pixdim(4);

if ~isempty(nii.hdr.hist.flip_orient)
    if nii.hdr.hist.flip_orient(1)
        header.mmppix(1) = header.mmx;
    else
        header.mmppix(1) = -header.mmx;
    end
    if nii.hdr.hist.flip_orient(2)
        header.mmppix(2) = header.mmy;
    else
        header.mmppix(2) = -header.mmy;
    end
    if nii.hdr.hist.flip_orient(3)
        header.mmppix(3) = header.mmz;
    else
        header.mmppix(3) = -header.mmz;
    end
else
    % If we have no information, we default to the "+ - -" convention.
    header.mmppix = [header.mmx, -header.mmy, -header.mmz];
end

%% Calculate center from sform.
sform = [nii.hdr.hist.srow_x;...
    nii.hdr.hist.srow_y;...
    nii.hdr.hist.srow_z];

lengths = [header.nVx, header.nVy, header.nVz, header.nVt];

% 1. Get "order" and "orientation" by reading "sform".
used = 0;
order = [1, 2, 3, 4];
for k = 1:2
    max = -1;
    m = -1;
    for l = 1:3
        if ~bitand(used, bitshift(1, l))  &&  (abs(sform(l, k)) > max)
            max = abs(sform(l, k));
            m = l - 1;
        end
    end
    used = bitor(used, bitshift(1, m));
    order(k) = m + 1;
end

switch used
    case 3
        order(3) = 3;
    case 5
        order(3) = 2;
    case 6
        order(3) = 1;
end

orientation = 0;
for k = 1:3
    if sform(order(k), k) < 0
        orientation = bitxor(orientation, bitshift(01, k));
    end
end

% 2. Create "revorder" from "order".
for k = 1:4
    revorder(order(k)) = k;
end

% 3. Bit operations on "orientation".
orientation = bitxor(orientation, bitshift(01, revorder(1)));
orientation = bitxor(orientation, bitshift(01, revorder(2)));

% 4. Flip the Y axis data in "sform", then reorder the sform axes per
% "order".
SS = zeros(size(sform));
for k = 1:3
    if bitand(orientation, bitshift(01, k))
        for l = 1:3
            sform(l, 4) = (lengths(k) - 1) * sform(l, k) + sform(l, 4);
            sform(l, k) = -sform(l, k);
        end
    end
end

for k = 1:3
    SS(order(k), :) = sform(k, :);
end
sform = SS;

% 5. Define "spacing".
spacing = zeros(4, 1);
spacing(1:3) = abs(diag(sform(1:3, 1:3)));
spacing(1:2) = -spacing(1:2);
spacing(4) = abs(header.nVt);

% 6. Initialize t4 matrix.
t4 = zeros(4);

% 7. Normalize "sform" rows by "spacing".
sform(1:3, 1:3) = sform(1:3, 1:3) ./ repmat(spacing(1:3)', 3, 1);

% 8. Calculate DET of sform.
DET = det(sform(1:3, 1:3));

% 9. Define "t4" as ADJUGATE(sform).
a = []; b = []; c = []; d = [];
for k = 1:3
    a = mod(k, 3) + 1;
    b = mod(k + 1, 3) + 1;
    for l = 1:3
        c = mod(l, 3) + 1;
        d = mod(l + 1, 3) + 1;
        t4(l, k) = sform(a, c) * sform(b, d) - sform(a, d) * sform(b, c);
    end
end

% 10. Divide "t4" by DET.
t4 = t4 / DET;

% 11. Define "center" as product of "t4" and "sform".
center = zeros(3, 1);
for k = 1:3
    for l = 1:3
        center(k) = center(k) + sform(l, 4) * t4(k, l);
    end
end

% 12. Invert "center" and shift by "spacing".
center(1:3) = -center(1:3) + spacing(1:3);

% 13. Flip X and Z axes in "spacing" and "center".
center(1) = -(spacing(1) * (lengths(revorder(1)) + 1) - center(1));
spacing(1) = -spacing(1);

center(3) = -(spacing(3) * (lengths(revorder(3)) + 1) - center(3));
spacing(3) = -spacing(3);

% Final step.
header.center = center';

%% Save original header.
header.original_header = rmfield(nii, 'img');


