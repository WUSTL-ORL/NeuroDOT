function [header_out, img_out] = nifti_4dfp(header_in, img_in, mode)
% Function to convert between nifti and 4dfp header formats.
% Inputs:
%   header_in: 'struct' of either a nifti or 4dfp header
%   mode: flag to indicate direction of conversion between nifti and 4dfp
%   formats 
%       - 'n': 4dfp to Nifti
%       - '4': Nifti to 4dfp
%
% Output:
%   header_out: header in Nifti or 4dfp format
%
%
%%
switch mode
    % Convert from 4dfp to Nifti style header
    case 'n'
        % Initialize default parameters
        AXES_4DFP =  [0,1,2,3];
        axes = AXES_4DFP;
        order = axes;
        f_offset = 0;
        isFloat = 1;
        isSigned = 1;
        bperpix = 4;
        number_of_dimensions = 4;
        timestep = 1;
        haveCenter = 1;
        scale = 1.0;
        offset = 0.0;
        
        % Check the orientation 
        if isfield(header_in, 'acq')
            if strcmp(header_in.acq, 'transverse')
                header_in.orientation = 2;
            elseif strcmp(header_in.acq, 'sagittal')
                header_in.orientation = 4;
            end
        end
        
        % Initialize sform matrix to zeros (3x4)
        for i = 0:2
            for j = 0:3
                sform(i+1,j+1) = 0.0; 
            end
        end
        
        % Assign diagonal of sform to ones
        for i = 0:2
            sform(order(i+1)+1,i+1) = 1.0;
        end
        
        
        % Adjust order of dimensions depending on the orientation
        switch header_in.orientation
            case 4
                tempv = order(2);
                order(2) = order(1);
                order(1) = tempv;
            case 3
                tempv = order(3);
                order(3) = order(2);
                order(2) = tempv;
            % case 2 is default case, so do nothing here
        end

        dims = number_of_dimensions;
        if isfield(header_in, 'nVx')
            header_in.matrix_size = [header_in.nVx, header_in.nVy, header_in.nVz, header_in.nVt];
        end
        dim(1) = header_in.matrix_size(1);
        dim(2) = header_in.matrix_size(2);
        dim(3) = header_in.matrix_size(3);
        dim(4) = header_in.matrix_size(4);

        spacing(1) = header_in.mmppix(1);
        spacing(2) = header_in.mmppix(2);
        spacing(3) = header_in.mmppix(3);

        center(1) = header_in.center(1);
        center(2) = header_in.center(2);
        center(3) = header_in.center(3);

%         disp('input: ')
%         disp(center)
%         disp('spacing: ')
%         disp(spacing)

        switch header_in.orientation
            case 2
            case 4
                center(3) = -center(3);
                spacing(3) = -spacing(3);
                center(3) = spacing(3) *(dim(3) +1)-center(3);
            case 3
                center(1) = -center(1);
                spacing(1) = -spacing(1);
                center(1) = spacing(1)*(dim(1) +1) - center(1);
        end

        % Adjust for fortran 1-indexing
        for i = 0:2
            center(i+1) = center(i+1) -spacing(i+1);
            center(i+1) = -center(i+1);
        end

        % Add sform to t4trans
        for i = 0:2
            for j = 0:3
                t4trans(i+1,j+1) = sform(i+1,j+1);
            end
        end

        for i = 0:2
            %apply spacing
            for j = 0:2
                sform(i+1,j+1) = t4trans(i+1,j+1)*spacing(j+1);
            end
            for j = 0:2
                sform(i+1,4) = sform(i+1,4) + center(j+1)* t4trans(i+1,j+1);
            end
        end

%         disp('sform output: ')
%         disp(sform)
        % Save Nifti

        % to_lpi
        used = 0;
        for i = 0:1
            max = -1.0;
            k = -1;
            for j = 0:2
                if ((bitcmp(bitand(used, (bitshift(1,j))))) & abs(sform(j+1,i+1)) > max)
                    max = abs(sform(j+1,i+1));
                    k = j;
                end
            end
            used = bitor(used, bitshift(1,k));
            order(i+1) = k;
        end
        switch used
            case 3
                order(3) = 2;
            case 5
                order(3) = 1;
            case 6
                order(3) = 0;
        end
%         order = [0,1,2,3];

        orientation = 0;
        for i = 0:2
            if sform(order(i+1)+1,i+1) < 0.0
                orientation = bitxor(orientation, bitshift(1, i));
            end
        end

        % auto_orient_header 
        for i = 0:2
            % Flip axes
            if bitand(orientation, (bitshift(1,i)))
                for j = 0:2
                    sform(j+1,4) = (dim(i+1)-1)*sform(j+1,i+1) + sform(j+1,4);
                    sform(j+1, i+1) = -sform(j+1,i+1);
                end
            end        
        end

        % Re order axes to x, y, z, t
        for i = 0:2
            for j =0:2
                sform(i+1,order(j+1)+1) = sform(i+1,j+1);
            end
        end

        % Load it back into the sform
        for i = 0:2
            for j = 0:2
                sform(i+1,j+1) = sform(i+1,j+1);
            end
        end

        % Define revorder
        revorder = zeros(1,4); % Line 425
        for i = 0:3            % Line 426
            revorder(order(i+1)+1) = i;
        end

%         disp('sform: ')
%         disp(sform)
        spacing = [0.0,0.0,0.0]; % Initialize spacing
        for i = 0:2
            for j = 0:2
                spacing(i+1) = spacing(i+1) + sform(j+1, i+1)*sform(j+1, i+1);
            end
            spacing(i+1) = sqrt(spacing(i+1));
        end
        
        % Create output Nifti-style header
        header_out = struct;
        header_out.dim = [...
            4, ...
            dim(revorder(1)+1)...
            dim(revorder(2)+1), ...
            dim(revorder(3) +1),...
            dim(revorder(4)+1),0,0,0];
        header_out.pixdim = [1.0,...
            spacing(1),...
            spacing(2), ...
            spacing(3), ...
            timestep, ...
            0,0,0];
        header_out.srow_x = [sform(1,1), sform(1,2),...
            sform(1,3), sform(1,4)];
        header_out.srow_y = [sform(2,1), sform(2,2),...
            sform(2,3), sform(2,4)];
        header_out.srow_z = [sform(3,1), sform(3,2),...
            sform(3,3), sform(3,4)];
        header_out.sform_code = 3;
        header_out.sizeof_hdr = 348;
        header_out.aux_file = '';
        header_out.descrip = [header_in.filename,'.4dfp.ifh converted with nifti_4dfp'];
        header_out.vox_offset = 352;
        NIFTI_UNITS_MM = 2;
        NIFTI_UNITS_SEC = 8;
        header_out.xyzt_units = NIFTI_UNITS_MM + NIFTI_UNITS_SEC;
        header_out.bitpix = 32;
        header_out.datatype = 16;
        header_out.dim_info = '';
        header_out.qform_code = 0;
        header_out.intent_p1 = 0;
        header_out.intent_p2 = 0;
        header_out.intent_p3 = 0;
        header_out.intent_code = 0;
        header_out.slice_start = 0;
        header_out.scl_slope = 0;
        header_out.scl_inter = 0;
        header_out.slice_end = 0;
        header_out.slice_code = 0;
        header_out.cal_max = 0;
        header_out.cal_min = 0;
        header_out.slice_duration = 0;
        header_out.toffset = 0;
        header_out.quatern_b = 0;
        header_out.quatern_c = 0;
        header_out.quatern_d = 0;
        header_out.qoffset_x = 0;
        header_out.qoffset_y = 0;
        header_out.qoffset_z = 0;
        header_out.intent_name = '';
        c = newline;
        header_out.magic = ['n+1', c];
%         disp('srow: ')
%         disp(header_out.srow_x)
%         disp(header_out.srow_y)
%         disp(header_out.srow_z)
        
        %% New dev 6/12/23
        outmem = zeros(size(img_in));
        orig_sform = sform;
        %% auto_orient
        nan_found = 0; i = 0; val_flip = zeros(1,4);
        in_val = zeros(4,1);
        out_val = zeros(4,1);
        target_length = zeros(4,1);
        in_length = header_out.dim(1:4);
        voxels = img_in;
        rData = voxels;
        for i = 0:3
            target_length(order(i+1)+1) = in_length(i+1);
            val_flip(i+1) = bitand(orientation, bitshift(1, i));
        end
        
        % Flip 
        if header_in.orientation == 2
            val_flip = zeros(1,4);
        end
        [~, idx_flip] = find(val_flip > 1);
        new_order = 1:ndims(img_in);
        if any(idx_flip) > 0
            idx_new = flip(idx_flip);
            new_order(idx_new) = flip(new_order(idx_new));           
            img_xfm = permute(img_in, new_order);
        else
            img_xfm = img_in;
        end
        if header_out.dim(5) <= 1
            val_flip = val_flip(1:3);
        end
        for k = 1:length(val_flip)
            if any(orig_sform(1:3,k) < 0)
                img_xfm = flip(img_xfm, k);
            end
        end    
         img_out = img_xfm;

    case '4'
        % Convert from Nifti to 4dfp style header

        % Look for the raw structure- if the raw is passed in, do as
        % normal; if the whole header, then look for
        % header_in.raw.
        %% Parse_Nifti
        % Initialize parameters
        AXES_NII = [0,1,2,3]; % From common-format.h 
        axes = AXES_NII; % nifti_format.c:342
        
        if isfield(header_in, 'raw')
            header_in = header_in.raw;
        end
        header_in.dim(1) = [];  % dim(1) contains the number of dimensions, which is not used in the computations

        % Parse_nifti
        order = axes; % nifti_format.c:344

        for i = 1:4  % nifti_format.c:371
            for j = 1:4 % nifti_format.c:372
                sform(i,j) = 0.0; % nifti_format.c:373
            end
        end
        if header_in.sform_code > 0  % nifti_format.c:377
            for i = 1:4 % nifti_format.c:379
                sform(1,i) = header_in.srow_x(i); % nifti_format.c:381
                sform(2,i) = header_in.srow_y(i); % nifti_format.c:382
                sform(3,i) = header_in.srow_z(i); % nifti_format.c:383
            end

        else
            if header_in.qform_code > 0 % nifti_format.c:406
                b = header_in.quatern_b;
                c = header_in.quatern_c;
                d = header_in.quatern_d;
                a = sqrt(1.0 - (b^2 + c^2 + d^2));
                % generate rotation matrix (Sform)
                sform(1,1) = a^2 + b^2 - c^2 - d^2;
                sform(1,2) = 2*b*c - 2*a*d;
                sform(1,3) = 2*b*d + 2*a*c;
                sform(2,1) = 2*b*c + 2*a*d;
                sform(2,2) = a^2 + c^2 - b^2 - d^2;
                sform(2,3) = 2*c*d - 2*a*b;
                sform(3,1) = 2*b*d - 2*a*c;
                sform(3,2) = 2*c*d + 2*a*b;
                sform(3,3) = a^2 + d^2 - c^2 - b^2;

                if header_in.pixdim(1) < 0.0
                    header_in.pixdim(4) = -header_in.pixdim(4);%/* read nifti1.h:1005 for yourself, i can't make this stuff up */
                end


                for j = 1:4
                    sform(j,1) = sform(j,1)*header_in.pixdim(2);
                    sform(j,2) = sform(j,2)*header_in.pixdim(3);
                    sform(j,3) = sform(j,3)*header_in.pixdim(4);
                end
                sform(1,4) = header_in.qoffset_x;
                sform(2,4) = header_in.qoffset_y;
                sform(3,4) = header_in.qoffset_z;   
            else %  do it the originless way
                sform(1,1) = header_in.pixdim(2);
                sform(2,2) = header_in.pixdim(3);
                sform(3,3) = header_in.pixdim(4);
            end
        end
%         disp(sform)

%% Save_4dfp (Line 412 in 4dfp_format.c) 
        % to_lpi (in transform.c) 
        used = 0;
        for i = 0:1
            max = -1.0;
            k = -1;
            for j = 0:2
                if ((bitcmp(bitand(used, (bitshift(1,j))))) & abs(sform(j+1,i+1)) > max)
                    max = abs(sform(j+1,i+1));
                    k = j;
                end
            end
            used = bitor(used, bitshift(1,k));
            order(i+1) = k;
        end
        switch used
            case 3
                order(3) = 2;
            case 5
                order(3) = 1;
            case 6
                order(3) = 0;
        end

        orientation = 0;
        for i = 0:2
            if sform(order(i+1)+1,i+1) < 0.0
                orientation = bitxor(orientation, bitshift(1, i));
            end
        end

        revorder = zeros(1,4); % 4dfp_format.c:425
        for i = 0:3 % 4dfp_format.c:426
            revorder(order(i+1)+1) = i;
        end
        orientation = bitxor(orientation, bitshift(1, revorder(1)));
        orientation = bitxor(orientation, bitshift(1, revorder(2)));
        %% Auto_orient_header
        temp_sform = zeros(3,4);
        orig_sform = sform;

        for i = 0:2
            % Flip axes
            if bitand(orientation, (bitshift(1,i)))
                for j = 0:2
                    sform(j+1,4) = (header_in.dim(i+1)-1)*sform(j+1,i+1) + sform(j+1,4);
                    sform(j+1, i+1) = -sform(j+1,i+1);
                end
            end        
        end
        
        % Re-order axes to x, y, z, t
        for i = 0:2
            for j =0:2
                temp_sform(i+1,order(j+1)+1) = sform(i+1,j+1);
            end
        end
        
        % Load it back into the sform
        for i = 0:2
            for j = 0:2
                sform(i+1,j+1) = temp_sform(i+1,j+1);
            end
        end
        % Initialize spacing
        spacing = [0,0,0,1];
        for i = 1:3
            for j = 1:3
                spacing(i) = spacing(i) + sform(j,i)^2;
            end
            spacing(i) = sqrt(spacing(i));
        end
        spacing(1) = -spacing(1);% keep the +, -, - convention in the .ifh, x and z are flipped later */
        spacing(2) = -spacing(2);% we do this here to specify we want the flips to take place before the t4 transform is applied */

        % Initialize t4trans
        t4trans = repmat(0.0, 4,4); % 4dfp_format.c:453
        t4trans(4,4) = 1.0; % 4dfp_format.c:456
        
        % First, invert sform to get t4
        for i = 1:3 % 4dfp_format.c:460
            for j = 1:3 % 4dfp_format.c:462
                sform (i,j) = sform(i,j)/spacing(j); % 4dfp_format.c:464
            end
        end
        
        
        %% Calculate Determinant of 3x3 rotation matrix
        determinant = 0;
        for i = 0:2 % 4dfp_format.c:467-477
            % Determinant
            temp = 1.0;
            temp2 = 1.0;
            for j = 0:2
                temp = temp*sform(j+1,mod((i+j),3)+1);
                temp2 = temp2*sform(j+1,mod(i-j+3,3)+1);
            end
            determinant = determinant + temp-temp2;
            temp = 1.0;
        end
        
        
        %% Adjugate
        % Since mod() is performing math, given the i,j in C, use the same i and j,
        % then add 1 at the end since a,b,c,d are indices
        t4trans = repmat(0.0, 4,4); % Line 4dfp_format.c:453
        t4trans(4,4) = 1.0; % 4dfp_format.c:456
        for i = 0:2 % 4dfp_format.c:480-488
            a = mod((i + 1),3)+1; 
            b = mod((i + 2),3)+1;
            for j = 0:2
                c = mod((j+1),3)+1;
                d = mod((j+2),3)+1;
                t4trans(j+1,i+1) =(sform(a,c)*sform(b,d)) - (sform(a,d) *  sform(b,c));
            end
        end
        
        
        %% Divide t4trans by determinant
        t4trans(1:3,1:3) = t4trans(1:3,1:3)./determinant;  

        
        %% Calculate center
        % Initialize center and assign values from sform multiplied by
        % t4trans (4dfp_format.c:497-505)
        center = [0,0,0];
        for i = 0:2
            for j = 0:2
                center(i+1) = center(i+1) + (sform(j+1,4)*t4trans(i+1,j+1));
            end
        end

        % center lines 4dfp_format.c:513-518
        for i = 0:2
            center(i+1) = -center(i+1);
            center(i+1) = center(i+1) + spacing(i+1);
        end

        % center 4dfp_format.c:522-527
        center(1) = (spacing(1) * (header_in.dim(revorder(1)+1)+1))-center(1);
        center(1) = -center(1);
        spacing(1) = -spacing(1);

        center(3) = (spacing(3) * (header_in.dim(revorder(3)+1)+1))-center(3);
        center(3) = -center(3);
        spacing(3) = -spacing(3);


        header_out = struct;
        header_out.matrix_size = [header_in.dim(revorder(1)+1),...
            header_in.dim(revorder(2)+1),...
            header_in.dim(revorder(3) + 1),...
            header_in.dim(revorder(4) + 1)];
        header_out.acq = 'transverse';
        header_out.nDim = 4;
        header_out.orientation = 2;
        header_out.scaling_factor = abs(spacing);
        header_out.mmppix = [spacing(1), spacing(2), spacing(3)];
        header_out.center = [center(1), center(2), center(3)];

        
        %% New dev 6/12/23

        outmem = zeros(size(img_in)); 
        %% auto_orient
        nan_found = 0; i = 0; val_flip = zeros(1,4);
        in_val = zeros(4,1);
        out_val = zeros(4,1);
        target_length = zeros(4,1);
        in_length = header_in.dim(1:4);
        voxels = img_in;
        rData = voxels;

        for i = 0:3
            target_length(order(i+1)+1) = in_length(i+1);
            val_flip(i+1) = bitand(orientation, bitshift(1, i));
        end
        
        % Flip 
        if orientation == 2
            val_flip = zeros(4,1);
        end
        [~, idx_flip] = find(val_flip > 1);
        new_order = 1:ndims(img_in);
        if any(idx_flip) > 0
            idx_new = flip(idx_flip);
            new_order(idx_new) = flip(new_order(idx_new));           
            img_xfm = permute(img_in, new_order);
        else
            img_xfm = img_in;
        end
        
        for k = 1:length(val_flip)
            if any(orig_sform(k, 1:3) < 0) && orig_sform(k, 4) > 0
                img_xfm = flip(img_xfm, k);
            end
        end  
         img_out = img_xfm;
end

end



