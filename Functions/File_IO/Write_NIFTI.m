function nii = Write_NIFTI(volume, header, filename)

% WRITE_NIFTI Assembles an "nii" structure from 4dfp data.
% 
%   nii = WRITE_NIFTI(volume, header, filename) converts the volumetric
%   data "volume" described by "header" from 4dfp to NIFTI format, output
%   as "nii".
% 
%   NOTE: The math in this function to calculate the "center" field for
%   4dfp is not straightforward, and is based on source code shared by Tim
%   Coalson from his "nifti-4dfp" toolbox.
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
nii = [];

if isfield(header, 'original_header')
    nii = header.original_header;
else
    % Fill in all fields that we don't care to mess with, but aren't
    % actually filled by NIFTI Reader. These are default values.
    
    % nii.filetype assigned in NIFTI Reader based on nii.hdr.hist.magic.
    % nii.fileprefix assigned below.
    % nii.machine assigned below.
    % nii.img assigned below.
    % nii.original is unique to NIFTI Reader, never saved in .nii file.
    
    nii.hdr.hk.sizeof_hdr = 348; % Default NIFTI header size.
    nii.hdr.hk.data_type = '';
    nii.hdr.hk.db_name = '';
    nii.hdr.hk.extents = 0;
    nii.hdr.hk.session_error = 0;
    nii.hdr.hk.regular = 'r';
    nii.hdr.hk.dim_info = 0;
    
    nii.hdr.dime.dim = zeros(1, 8); % Initializing for below.
    nii.hdr.dime.intent_p1 = 0;
    nii.hdr.dime.intent_p2 = 0;
    nii.hdr.dime.intent_p3 = 0;
    nii.hdr.dime.intent_code = 0;
    % nii.hdr.dime.datatype assigned below.
    % nii.hdr.dime.bitpix assigned in NIFTI Reader.
    nii.hdr.dime.slice_start = 0;
%     nii.hdr.dime.pixdim = [1, 1, 1, 1, 1, 0, 0, 0]; % This *would* be read
    % from mmppix, but NIFTI requires these values anyways.
    nii.hdr.dime.vox_offset = 0;
    nii.hdr.dime.scl_slope = 1;
    nii.hdr.dime.scl_inter = 0;
    nii.hdr.dime.slice_end = 0;
    nii.hdr.dime.slice_code = 0;
    nii.hdr.dime.xyzt_units = 10; % This is an arbitrary value based on our sample set: YMMV.
    nii.hdr.dime.cal_max = 0;
    nii.hdr.dime.cal_min = 0;
    nii.hdr.dime.slice_duration = 0;
    nii.hdr.dime.toffset = 0;
    % nii.hdr.dime.glmax assigned in NIFTI Reader.
    % nii.hdr.dime.glmin assigned in NIFTI Reader.
    
    % nii.hdr.hist.descrip assigned below.
    nii.hdr.hist.aux_file = '';
    nii.hdr.hist.qform_code = 0;
    nii.hdr.hist.sform_code = 0;
    nii.hdr.hist.quatern_b = 0; % These are arbitrarily based on sample set. Not intended to ever get used in actuality.
    nii.hdr.hist.quatern_c = 1;
    nii.hdr.hist.quatern_d = 0;
    % nii.hdr.hist.qoffset_x assigned below.
    % nii.hdr.hist.qoffset_y assigned below.
    % nii.hdr.hist.qoffset_z assigned below.
    % nii.hdr.hist.srow_x assigned below.
    % nii.hdr.hist.srow_y assigned below.
    % nii.hdr.hist.srow_z assigned below.
    nii.hdr.hist.intent_name = '';
    nii.hdr.hist.magic = 'n+1'; % Writing in latest version of NIFTI.
    % nii.hdr.hist.originator assigned below.
    % nii.hdr.hist.rot_orient is specific to NIFTI Reader, never saved.
    % nii.hdr.hist.flip_orient is specific to NIFTI Reader, never saved.
end

switch class(volume)
    case 'uint8'
        nii.hdr.dime.datatype = 2;
    case 'int16'
        nii.hdr.dime.datatype = 4;
    case 'int32'
        nii.hdr.dime.datatype = 8;
    case 'single'
        nii.hdr.dime.datatype = 16;
    case 'double'
        nii.hdr.dime.datatype = 64;
    case 'int8'
        nii.hdr.dime.datatype = 256;
    case 'uint16'
        nii.hdr.dime.datatype = 512;
    case 'uint32'
        nii.hdr.dime.datatype = 768;
    case 'uint64'
        nii.hdr.dime.datatype = 1280;
    case 'int64'
        nii.hdr.dime.datatype = 1024;
end

nii.hdr.hist.descrip = 'NeuroDOT_2_SaveVolumetricData';

%% Do math to convert "center" into NIFTI format.

% 1. Initialize "order", "lengths", and "sform".
lengths = [header.nVx, header.nVy, header.nVz];
if isfield(header, 'nVt')
    lengths(4) = header.nVt;
else
    lengths(4) = 1;
end
order = [1, 2, 3, 4];

sform = zeros(3, 4);

% 2. Swap "order" values based on orientation.
switch header.acq
    case 'sagittal' % 4
        temp = order(2);
        order(2) = order(1);
        order(1) = temp;
    case 'coronal' % 3
        temp = order(3);
        order(3) = order(2);
        order(2) = temp;
    case 'transverse' % 2
        % Do nothing.
end

% 3. Fill in "sform".
for k = 1:3
    sform(order(k), k) = 1;
end

% 4. Initialize "center" and "spacing".
spacing = zeros(4,1);
spacing(1) = header.mmppix(1);
spacing(2) = header.mmppix(2);
spacing(3) = header.mmppix(3);

center = zeros(4,1);
center(1) = header.center(1);
center(2) = header.center(2);
center(3) = header.center(3);

% 5. Flip "center" and "spacing" based on orientation.
switch header.acq
    case 'sagittal' % 4
        center(3) = -center(3);
        spacing(3) = spacing(3);
        center(3) = spacing(3) * (lengths(3) + 1) - center(3); % file_content.length(3) IS nVz
    case 'coronal' % 3
        center(1) = -center(1);
        spacing(1) = -spacing(1);
        center(1) = spacing(1) * (lengths(1) + 1) - center(1);
    case 'transverse' % 2
end

% 6. Adjust "center" to account for 1-indexing.
center(1:3) = -(center(1:3) - spacing(1:3));

% 7. Apply "spacing" and "center" to "sform".
t4 = sform;
sform(1:3, 1:3) = t4(1:3, 1:3) .* repmat(spacing(1:3)', 3, 1);
sform(1:3, 4) = sform(1:3, 4) + sum(t4(1:3, 1:3) .* repmat(center(1:3)', 3, 1), 2);

% 8. Get "order" and "orientation" from "sform", and create "revorder".
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

for k = 1:4
    revorder(order(k)) = k;
end

% 9.  Flip the Y axis data in "sform", then reorder the sform axes per
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

% We can assign stuff now!
nii.hdr.hist.qoffset_x = center(revorder(1));
nii.hdr.hist.qoffset_y = center(revorder(2));
nii.hdr.hist.qoffset_z = center(revorder(3));
nii.hdr.hist.srow_x = sform(1, :);
nii.hdr.hist.srow_y = sform(2, :);
nii.hdr.hist.srow_z = sform(3, :);
nii.hdr.hist.sform_code = 2;
nii.hdr.hist.qform_code = 0;

nii.hdr.dime.dim(1) = header.nDim;
nii.hdr.dime.dim(2) = lengths(revorder(1));
nii.hdr.dime.dim(3) = lengths(revorder(2));
nii.hdr.dime.dim(4) = lengths(revorder(3));
nii.hdr.dime.dim(5) = lengths(revorder(4));

nii.hdr.dime.pixdim=[1,header.mmx,header.mmy,header.mmz,1,0,0,0];

%% Other stuff.
nii.fileprefix = filename;
nii.machine = ['ieee-', header.byte, 'e'];
nii.img = volume;



%
