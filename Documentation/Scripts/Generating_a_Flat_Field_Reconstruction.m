% Example code to generate and plot a Flat Field Reconstruction (FFR) using
% a 24x28 pad. The required A matrix can be downloaded from NITRC:
%
%   https://www.nitrc.org/frs/download.php/12843/A_AdultV24x28.mat
%
% Explanation of the FFR can be found as part of the light modeling tutorial:
% 
%  NeuroDOT_Tutorial_Generating_a_Light_Model_Pad_24x28_With_AlignMe.pptx
%
% The FFR code from that tutorial is reproduced here for convenience.

% generate the FFR

load('A_AdultV24x28.mat','info','A');
keep = (info.pairs.WL==2) & (info.pairs.r3d<40);
A = squeeze(A(keep,:));
iA = Tikhonov_invert_Amat(A,0.01,0.1);
iA = smooth_Amat(iA,info.tissue.dim,5); 
FFR = makeFlatFieldRecon(A,iA);
FFR_as_vol = Good_Vox2vol(FFR,info.tissue.dim);

% load an atlas to plot the FFR on

template_fname = 'Segmented_MNI152nl_on_MNI111_nifti.nii';
[template,template_info] = LoadVolumetricData(template_fname);

% transform the atlas to the data space and plot

template_clipped = affine3d_img(template, template_info, info.tissue.dim, [], 'nearest'); 
PlotSlices(template_clipped,info.tissue.dim,[],FFR_as_vol)

% transform the data to the atlas space and plot

FFR_in_wholebrainspace = affine3d_img(FFR_as_vol, info.tissue.dim, template_info);
PlotSlices(template,template_info,[],FFR_in_wholebrainspace)

