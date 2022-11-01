function [mask]=threshold5R_new(T1w,T2w)

% This segmentation procedure works by an iteration of masking and
% thresholding with region growing.  
%
% Step 0: T1 and T2 volumes are put into same resolution, FOV, and perspective
% Step 1: Volumes are individually normalized
% Step 2: The whole head is separated from the background using T1
% Step 3: The scalp/skull region is separated from the Brain/CSF region
%           using a combination of T1 and T2
% Step 4: The scalp is separated from the skull using T1
% Step 5: CSF is extracted using T1 and T2
% Step 6: Separate the WM and GM using Otsu's method in the brain-masked
%           T1.  Force GM and WM to be composed of a single contiguous region each.
% Step 7: Unlabeled voxels are identified and labeled based on neighborhood
%           population

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

%% set parameters
CSF_T2_thresh=0.175;                             % CSF T2 thresh
CSF_T1_thresh=0.175;                              % CSF T1 thresh
Skull_th=0.075;                                  % Skull T1/T2 thresh

Rbk=2;SE_bk_shrink = strel('ball', Rbk, 0);    % Head shrink 1 H/BG
Rhc=2;SE_head_close=strel('ball',Rhc,0);       % Head close

Rs=2;SE_h_shrink = strel('ball', Rs, 0);       % Head shrink 1 SMALL HEAD
Rs2=3;SE_h_shrink2 = strel('ball', Rs2, 0);     % head shrink 2 SHRINK
Rs3=3;SE_h_shrink3 = strel('ball', Rs3, 0);     % head dilate 3 DILATE
Rb=3;SE_b = strel('ball', Rb, 0);               % bknd dilate 
Rg=2;SE_g = strel('ball', Rg, 0);               % GM dilate 
Rbk=3;SE_bk_shrink2 = strel('ball', Rbk, 0);    % Final T1 Head shrink 1
killY=[1:20];                                   % Remove extra neck, etc


%% set up T1 and T2
T1a=T1w;
T1a=T1a./max(T1a(:)); % normalize

T2a=T2w;
T2a=T2a./max(T2a(:));


%% separate head from background
% Replace ball-dilates with this!!: OutVol=Points2Spheres(InVol,radius)
head_T1_bw=ExtractHead(T1a);        % Extract T1 head shape
head_T1_bw=imfill(head_T1_bw,'holes');
head_T2a_data=head_T1_bw.*T2a;      % mask T2 head shape
head=+(head_T2a_data>0);               % mask FOV based on T2 FOV
head=imfill(head,'holes');
bknd=1-head;                        % define background
bknd=+imdilate(logical(bknd),SE_bk_shrink); % Dilate
bknd=+(imclose(bknd,SE_head_close));
head=1-bknd;                       
head=smooth3(head,'gaussian',[7,7,7],2);    % Smooth head
head=+(head>0.5);
bknd=1-head;


head_T1a_data=T1a.*head;            % masked T1
head_T2a_data=T2a.*head;            % masked T2


%% separate scalp and skull from rest of head
Brain_CSF_temp=head_T1a_data.*head_T2a_data;            % T1 & T2
Brain_CSF_temp=Brain_CSF_temp./max(Brain_CSF_temp(:));  % Normalize
Brain_CSF_thresh=graythresh(Brain_CSF_temp);            % Graythresh
Brain_CSF_bw=+(Brain_CSF_temp>=Brain_CSF_thresh);       % B/W
head_small=+imerode(logical(head),SE_h_shrink);         % Small head
Brain_CSF_bw=Brain_CSF_bw.*head_small;                  % Brain_CSF
Brain_CSF_bw=imfill(Brain_CSF_bw,'holes');              % Fill
Brain_CSF_bw=SingleContiguousRegion(Brain_CSF_bw);      % SCR

Brain_CSF_bw=+imerode(logical(Brain_CSF_bw),SE_h_shrink2);  % Erode 2
Brain_CSF_bw=SingleContiguousRegion(Brain_CSF_bw);          % SCR
Brain_CSF_bw=imclose(Brain_CSF_bw,SE_head_close);       % Fill
Brain_CSF_bw=imfill(Brain_CSF_bw,'holes');              % Fill
Brain_CSF_bw=+imdilate(logical(Brain_CSF_bw),SE_h_shrink3); % Dilate
%Brain_CSF_bw=Points2Spheres(Brain_CSF_bw,Rs3);
Brain_CSF_bw=smooth3(Brain_CSF_bw,'gaussian',[7,7,7],2);    % Smooth brain
Brain_CSF_bw=+(Brain_CSF_bw>0.5);
Brain_CSF_bw=SingleContiguousRegion(Brain_CSF_bw);          % SCR

Brain_CSF_T1_data=Brain_CSF_bw.*T1a;          % T1 data
Brain_CSF_T2_data=Brain_CSF_bw.*T2a;          % T2 data

Scalp_skull_bw=head-Brain_CSF_bw;                       % Scalp/skull BW
Scalp_skull_T1_data=Scalp_skull_bw.*head_T1a_data;      % Scalp-skull T1


%% separate the scalp from the skull using the T1 volume
ScSk_thresh=graythresh(Scalp_skull_T1_data(Scalp_skull_T1_data~=0));
% Scalp_bw=+(Scalp_skull_T1_data>=Skull_th);        
Scalp_bw=+(Scalp_skull_T1_data>=ScSk_thresh);              % Scalp BW
Skull_bw=Scalp_skull_bw-Scalp_bw;                       % Skull BW
Skull_bw=SingleContiguousRegion(Skull_bw);              % Skull SCR
bkgnd_perim=imdilate(bknd,SE_b);                        % Dilate bknd
Skull_bw=Skull_bw.*(1-bkgnd_perim);                     % Remove outer skull
Skull_bw=SingleContiguousRegion(Skull_bw);              % Skull SCR
Scalp_bw=Scalp_skull_bw-Skull_bw;                       % New Scalp


%% separate CSF from head using T1 & T2 volumes
CSF_bw=zeros(size(Brain_CSF_T1_data));                  % Initialize
keep=intersect(find(Brain_CSF_T2_data>=CSF_T2_thresh),...
    find(Brain_CSF_T1_data<=CSF_T1_thresh));            % Threshold T1 and T2
CSF_bw(keep)=1;                                         % Set


%% separate brain into gray matter (GM) and white matter (WM)
Brain_bw=Brain_CSF_bw-CSF_bw;
Brain_T1_data=Brain_bw.*head_T1a_data;
Brain_T2_data=Brain_bw.*head_T2a_data;

% WM
WM_threshA=graythresh(Brain_T1_data(Brain_T1_data>0));
WM_threshM=0.3;
WM_thresh=mean([WM_threshA,WM_threshM]);
WM_bw=+(Brain_T1_data>=WM_thresh);
WM_bw=SingleContiguousRegion(WM_bw);              % WM SCR

% GM
GM_bw=Brain_bw-WM_bw;
GM_bw=imclose(GM_bw,SE_g);                       % Close GM
GM_bw=imfill(GM_bw,'holes');                     % Fill
WM_bw=Brain_bw-GM_bw;                            % Fix WM


%% Put it all together

% Set un-segmented T1 to scalp.
head_T1_bwC=imclose(head_T1_bw,strel('ball',4,0));
head_T1_bwC=imfill(head_T1_bwC,'holes');
bknd_T1=1-head_T1_bwC;                        % define background
bknd_T1=+imdilate(logical(bknd_T1),SE_bk_shrink2); % Dilate
head_T1=1-bknd_T1;                        % Fix Head
head_T1=smooth3(head_T1,'gaussian',[7,7,7],2);    % Smooth head
head_T1=+(head_T1>0.5);

% Set values to region numbers: [1,2,3,4,5]=[csf,wm,gm,sk,sc]
mask=head_T1.*5;        %0.6

mask(Scalp_bw==1)=5;  %0.6
mask(Skull_bw==1)=4;  %0.7
mask(GM_bw==1)=3;     %0.8
mask(WM_bw==1)=2;     %0.9
mask(CSF_bw==1)=1;    %1.0

mask(:,killY,:)=0;

%{
% create gray values for each tissue type
MaxT1=max(T1w(:));
mask=zeros(size(T1a));
mask(bknd>0)=MaxT1 * 0.5;
mask(Scalp_bw>0)=MaxT1 * 1.2;
mask(Skull_bw>0)=MaxT1 * 1.4;
mask(CSF_bw>0)=MaxT1 * 1.8;
mask(WM_bw>0)=MaxT1 * 2.0;
mask(GM_bw>0)=MaxT1 * 1.6;
View4dfp(mask,'sagittal');


empty=find(mask==0);
disp(['Labeling ',num2str(length(empty)),' unlabeled voxels'])
tic
if numel(empty)>0
    for i=1:length(empty)
        temp=zeros(size(mask));
        temp(empty(i))=1;
           
        temp_neighbor=imdilate(temp,ones(3,3,3)); %compare with 26 nns
        neighbor=mask(temp_neighbor>0);
        neighbor(neighbor==0)=[];
                
        [uniques,numUnique] = count_unique(neighbor);
        [~,Imax]=max(numUnique);
        mask(empty(i))=uniques(Imax);    
    end
    empty2=mask==0;
    mask(empty2)=MaxT1*1.8;
end
toc
mask(bkgnd>0)=MaxT1 * 0.0;

Stats.scalp_fract=length(find(mask==1.2*MaxT1))/length(find(mask~=0));
Stats.skull_fract=length(find(mask==1.4*MaxT1))/length(find(mask~=0));
Stats.csf_fract=length(find(mask==1.8*MaxT1))/length(find(mask~=0));
Stats.gm_fract=length(find(mask==1.6*MaxT1))/length(find(mask~=0));
Stats.wm_fract=length(find(mask==2.0*MaxT1))/length(find(mask~=0));
figure;
bar([Stats.scalp_fract,Stats.skull_fract,Stats.csf_fract,Stats.gm_fract,Stats.wm_fract]);
set(gca,'XTickLabel',{'Scalp','Skull','CSF','GM','WM'});
ylabel('Fraction Head Volume')
clear T1a
%}