function IM_d1f=ExtractHead(T1,th,R)

% Extract adult head from background of a T1-weighted image


%% Parameters
if ~exist('th','var')
Foreground_th=0.25; % Foreground_th=0.25; % MNI!
else
    Foreground_th=th;
end
if ~exist('R','var'),R=4;end

center=R+1;
NSsize=(2*R)+1;
seed=zeros(NSsize,NSsize,NSsize);
seed(center,center,center)=1;
kernel=Points2Spheres(seed,R);

head_T1_bw=zeros(size(T1)); % Initialize

%%
T1=T1./max(T1(:));                                   % Normalize
head_T1_bw(T1>Foreground_th)=1;                      % Threshold
IM_d1 = imdilate(head_T1_bw,kernel);
IM_d1=imfill(IM_d1,26,'holes');           % Fill  holes
% head_T1_bw=+(convn(logical(head_T1_bw),kernel,'same')>=1); % Dilate

L=bwlabeln(IM_d1);                      % Only 1 contiguous head region
stats=regionprops(L,'Area');
A=[stats.Area];
biggest=find(A==max(A)); 
IM_d1(L~=max(biggest))=0;

IM_d1=imfill(IM_d1,26,'holes');           % Fill remaining holes
IM_d1=imfill(IM_d1,26,'holes');           % Fill remaining holes
IM_d1 = imerode(IM_d1,kernel);
IM_d1f = imclose(IM_d1,kernel);

L=bwlabeln(IM_d1f);                      % Only 1 contiguous head region
stats=regionprops(L,'Area');
A=[stats.Area];
IM_d1f(L~=find(A==max(A(:))))=0;