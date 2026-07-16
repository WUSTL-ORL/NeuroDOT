% Example code from NeuroDOT_Tutorials_Spatial_Transformations.m
%
% This file includes the code from the Appendix examples in the PowerPoint
% Tutorial demonstrating use of affine3d_img. 
%
% The three examples used in the main section of the tutorial are taken 
% from larger processing pipelines. The relevant code may be found in the 
% following mfiles in the NeuroDOT documentation:
% 
% Example 1 -- NeuroDOT visualization pipeline
% .../Documentation/Scripts/NeuroDOT_Visualization_Script.m
%
% Example 2 -- Visualizing the Flat Field Reconstruction
% .../Documentation/Scripts/Generating_a_Flat_Field_Reconstruction.m
%
% Example 3 -- Mapping to subject-specific anatomy
% .../Documentation/Scripts/Script_for_Participant_Specific_Head_Modeling.m
%
% The code snippets below are shown in Appendix A2 and A3.
%
% Each example is self-contained. To run, either highlight code, then 
% option-click and select "Evaluate Selection in Command Window" or click 
% on snippet %%title line and use the "Run Section" button in the toolbar.

%% Example 1: Translate the x-origin

[imgA,infoA] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
PlotSlices(imgA);

infoB = infoA; 					                % initialize infoB
infoB.center(1) = infoB.center(1) + 50; 		% change the relevant field
imgB = affine3d_img(imgA,infoA,infoB); 		    % generate the new image data
PlotSlices(imgB);					            % plot the transformed template

%% Example 2: Translate the y-origin

[imgA,infoA] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
PlotSlices(imgA);

infoB = infoA; 					                % initialize infoB
infoB.center(2) = infoB.center(2) + 50; 		% change the relevant field
imgB = affine3d_img(imgA,infoA,infoB); 		    % generate the new image data
PlotSlices(imgB);					            % plot the transformed template

%% Example 3: Translate the z-origin

[imgA,infoA] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
PlotSlices(imgA);

infoB = infoA; 					                % initialize infoB
infoB.center(3) = infoB.center(3) + 50; 		% change the relevant field
imgB = affine3d_img(imgA,infoA,infoB); 		    % generate the new image data
PlotSlices(imgB);					            % plot the transformed template

%% Example 4: Change the image width

[imgA,infoA] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
PlotSlices(imgA);

infoB = infoA; 					                % initialize infoB
infoB.nVx = infoA.nVx / 2;				        % change the relevant field
imgB = affine3d_img(imgA,infoA,infoB); 		    % generate the new image data
PlotSlices(imgB);					            % plot the transformed template

%% Example 5: Change the image depth

[imgA,infoA] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
PlotSlices(imgA);

infoB = infoA; 					                % initialize infoB
infoB.nVy = infoA.nVy / 2;				        % change the relevant field
imgB = affine3d_img(imgA,infoA,infoB); 		    % generate the new image data
PlotSlices(imgB);					            % plot the transformed template

%% Example 6: Change the image height

[imgA,infoA] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
PlotSlices(imgA);

infoB = infoA; 					                % initialize infoB
infoB.nVz = infoA.nVz / 2;				        % change the relevant field
imgB = affine3d_img(imgA,infoA,infoB); 		    % generate the new image data
PlotSlices(imgB);					            % plot the transformed template

%% Example 7: Change voxel resolution

[imgA,infoA] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
PlotSlices(imgA);

infoB = infoA; 					                % initialize infoB
infoB.nVx=infoA.nVx/10;infoB.nVy=infoA.nVy/10;infoB.nVz=infoA.nVz/10;
infoB.mmppix=10*infoA.mmppix; 			        % change the relevant field
imgB = affine3d_img(imgA,infoA,infoB); 		    % generate the new image data
PlotSlices(imgB);					            % plot the transformed template

%% Example 8: Scale by 1/2

[imgA,infoA] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
PlotSlices(imgA);

infoB = infoA; 					                % initialize infoB
t4 = [2 0 0 0; 0 2 0 0; 0 0 2 0; 0 0 0 1];	    % scale x 0.51
imgB = affine3d_img(imgA,infoA,infoB,t4); 	    % generate the new image data
PlotSlices(imgB);					            % plot the transformed template

%% Example 9: Scale by 2

[imgA,infoA] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
PlotSlices(imgA);

infoB = infoA; 					                % initialize infoB
t4 = [.5 0 0 0; 0 .5 0 0; 0 0 .5 0; 0 0 0 1];	% scale x 2
imgB = affine3d_img(imgA,infoA,infoB,t4); 	    % generate the new image data
PlotSlices(imgB);					            % plot the transformed template

%% Example 10: Rotate about x-axis

[imgA,infoA] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
PlotSlices(imgA);

infoB = infoA; 					                % initialize infoB
theta = pi/4; c = cos(theta); s = sin(theta); 	% rotate by 45 degrees 
t4 = [ 1 0 0 0 ; 0 c s 0 ; 0 -s c 0; 0 0 0 1];	
imgB = affine3d_img(imgA,infoA,infoB,t4); 	    % generate the new image data
PlotSlices(imgB);					            % plot the transformed template

%% Example 11: Rotate about y-axis

[imgA,infoA] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
PlotSlices(imgA);

infoB = infoA; 					                % initialize infoB
theta = pi/4; c = cos(theta); s = sin(theta); 	% rotate by 45 degrees 
t4 = [ c 0 -s 0; 0 1 0 0; s 0 c 0; 0 0 0 1 ];	
imgB = affine3d_img(imgA,infoA,infoB,t4); 	    % generate the new image data
PlotSlices(imgB);					            % plot the transformed template

%% Example 12: Rotate about z-axis

[imgA,infoA] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
PlotSlices(imgA);

infoB = infoA; 					                % initialize infoB
theta = pi/4; c = cos(theta); s = sin(theta); 	% rotate by 45 degrees 
t4 = [ c -s 0 0; s c 0 0; 0 0 1 0; 0 0 0 1 ];	
imgB = affine3d_img(imgA,infoA,infoB,t4); 	    % generate the new image data
PlotSlices(imgB);					            % plot the transformed template

%% Example 13: Shear along x

[imgA,infoA] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
PlotSlices(imgA);

infoB = infoA; 					% initialize infoB
t4 = [1 .5 0 0 ; 0 1 0 0; 0 0 1 0; 0 0 0 1];	% Sxy = 0.5
imgB = affine3d_img(imgA,infoA,infoB,t4); 	% generate the new image data
PlotSlices(imgB);					% plot the transformed template

%% Example 14: Shear along y

[imgA,infoA] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
PlotSlices(imgA);

infoB = infoA; 					                % initialize infoB
t4 = [1 0 0 0 ; .5 1 0 0; 0 0 1 0; 0 0 0 1];	% Syx = 0.5
imgB = affine3d_img(imgA,infoA,infoB,t4); 	    % generate the new image data
PlotSlices(imgB);					            % plot the transformed template

%% Example 15: Translate the x origin

[imgA,infoA] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
PlotSlices(imgA);

infoB = infoA; 					                % initialize infoB
t4 = [1 0 0 -50; 0 1 0 0; 0 0 1 0; 0 0 0 1];	% tx = 50
imgB = affine3d_img(imgA,infoA,infoB,t4); 	    % generate the new image data
PlotSlices(imgB);					            % plot the transformed template

%% Example 16: Scale and translate

[imgA,infoA] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
PlotSlices(imgA);

infoB = infoA; 					                % initialize infoB
t4 = [2 0 0 -50; 0 2 0 0; 0 0 2 0; 0 0 0 1];	% scale *and* translate
imgB = affine3d_img(imgA,infoA,infoB,t4); 	    % generate the new image data
PlotSlices(imgB);					            % plot the transformed template

%% Example 17: Demonstrating interpolation effects

[~,info111] = LoadVolumetricData('Segmented_MNI152nl_on_MNI111_nifti.nii');
[img333,info333] = LoadVolumetricData('Segmented_MNI152nl_on_333_nifti.nii');

imgA = affine3d_img(img333,info333,info111,[],'nearest');
PlotSlices(imgA);

imgB = affine3d_img(img333,info333,info111,[],'linear');
PlotSlices(imgB);
