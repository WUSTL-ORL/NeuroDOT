function [stats,PSFidx]=PSF_Metrics_200630(PSF,dim,GVidx,nSig)
%
% This function calculates the full-width at half maximum, localization
% error, and effective resolution for a point spread function,
% PSF=ia*a(:,vox) at a voxel within dim.Good_Vox.
% It is assumed that the dim structure contains coordinates for
% each Good_Vox index in dim.GVdCoord. GVidx is the dim.Good_Vox index at
% which the perturbation was modeled.
% This function calculates these metrics for any number of input
% perturbations: size(GVidx,1)==size(PSF,2).
% Inputs:
%   PSF     point-spread-functions (Nvox by Npsf)
%   dim     meta data of space for PSFs
%   GVidx   indices of dim.Good_Vox corresponding to PSF location
%   nSig    Noise std image volume (number of Good_Vox by 1)
%
% Outputs:
%   stats   structure containing metrics: fwhm, LocError, EffRes, SNR
%   PSFidx  cell array of indices of non-zero PSF values

%
% Copyright (c) 2017 Washington University
% Created By: Adam T. Eggebrecht & Jason Trobaugh
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
th=0.5;
Npsf=size(GVidx,1);
fvhm=zeros(Npsf,1,'single');
fwhm=zeros(Npsf,1,'single');
LocError=zeros(Npsf,1,'single');
EffRes=zeros(Npsf,1,'single');
MeanAct=zeros(Npsf,1,'single');
PSFidx=cell(Npsf,1);

%% Calculate % section indented 6/30/20 JWT, and added some comments below
parfor j=1:Npsf
    pt=dim.GVdCoord(GVidx(j),:);
    ptIdx=dim.Good_Vox(GVidx(j));
    PSFj=squeeze(PSF(:,j));                        % Grab PSF
    tempMax=max(PSFj);                          % max amplitude in PSF
    PSFj=PSFj./tempMax;                          % Normalize
    keep=PSFj>=th;                                 % Threshold
    keepIdx=find(keep);                            % indices of 'activation'
    foo=zeros(dim.nVx,dim.nVy,dim.nVz,'single');
    foo(dim.Good_Vox(keep))=1;                        % Expand to full vol
    CC = bwconncomp(foo,26);                    % test connectivity
    test=cellfun(@(x) any(x==ptIdx),CC.PixelIdxList);  % grab island with pt
    if any(test)
        [~,idx]=max(test);
    else
        [~,idx]=max(cellfun(@numel,CC.PixelIdxList));  % grab largest
    end
    keep1=ismember(dim.Good_Vox(keep),CC.PixelIdxList{idx});
    PSFj=PSFj(keepIdx(keep1)); % reduce PSF to voxels above th
    MeanAct(j)=mean(PSFj.*(tempMax)); % mean of PSF above threshold, after removing normalization 
    SigTotal(j)=sum(PSFj.*(tempMax)); % integrated signal over above-threshold region
    r=dim.GVdCoord(keepIdx(keep1),:);      % Coords of points above th
    PSFidx{j}=keepIdx(keep1);
    x=r(:,1);y=r(:,2);z=r(:,3);
    fvhm(j)=numel(x)*abs(prod(dim.mmppix));      % simple fvhm via counting % adapted 6/30/20 for units in mm, JWT
    fwhm(j)=max(max(pdist2(r,r)));                 % fwhm
    Centroid=[mean((x-mean(x)).*PSFj)+mean(x),...  % Centroid
        mean((y-mean(y)).*PSFj)+mean(y),...
        mean((z-mean(z)).*PSFj)+mean(z)];
    LocError(j)=norm(Centroid-pt);                 % Localization Error
    EffRes(j)=2*max(pdist2(r,pt));                 % Effective Resolution
end

%% Outputs
stats.fvhm=fvhm;
clear fvhm
stats.fwhm=fwhm;
clear fwhm
stats.LocError=LocError;
clear LocError
stats.EffRes=EffRes;
clear EffRes
stats.MeanAct=MeanAct;
stats.SigTotal = SigTotal; 

%% SNR if nSig is passed in
if exist('nSig','var')
    if numel(nSig)>0
        parfor j=1:Npsf
            SNR(j)=MeanAct(j)/mean(nSig(PSFidx{j}));
            SNRstd(j)=MeanAct(j)/std(nSig(PSFidx{j}));
            NoiseMean(j) = mean(nSig(PSFidx{j})); 
        end
    end
    stats.SNR=SNR;
    stats.SNRstd=SNRstd;
    stats.NoiseMean = NoiseMean; 
end



