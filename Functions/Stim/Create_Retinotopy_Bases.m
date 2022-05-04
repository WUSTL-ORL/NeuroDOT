function stimbase=Create_Retinotopy_Bases(stimsize,NumChecksRad,NumChecksTH)

% This function generates a set of basis arrays from which to generate
% traditional retinotopic stimuli.
%
% stimbase.x=pixel x-position
% stimbase.y=pixel y-position
% stimbase.r=pixel radius
% stimbase.th=pixel angle
% stimbase.On=image BW without crosshair
% stimbase.Off=image BW without crosshair
% stimbase.blank=blank gray screen without crosshair

% stimbase(:,:,8)=image A with crosshair
% stimbase(:,:,9)=image B with crosshair
% stimbase(:,:,10)=blank gray screen with crosshair
% Keep in mind the the angle goes from 0 to <360 with 0 along the -y axis
% and that the radius of the crosshair is 30 pixels.
%
% For PTB3, images must be nXn, uint8, btwn 0-255, and uint8.


%% Parameters and Initialization
if ~exist('stimsize','var'),stimsize=1023;end         % Granularity of stim images
if ~exist('NumChecksRad','var'),NumChecksRad=16; end  % Number of checks in radial direction; hcp=12
if ~exist('NumChecksTH','var'),NumChecksTH=24; end    % Number of check in the angular direction; hcp=24

maxrad=floor(stimsize/2)+1;
Dtheta=360/NumChecksTH;     % Angular width of checks (degrees)
rad(100)=maxrad;            % Maximum radius of stimulation

stimbase=struct('x',zeros(stimsize,stimsize),'y',zeros(stimsize,stimsize),...
    'r',zeros(stimsize,stimsize),'th',zeros(stimsize,stimsize),...
    'On',zeros(stimsize,stimsize),'Off',zeros(stimsize,stimsize),...
    'Blank',ones(stimsize,stimsize).*0.5);


%% Generate bases
disp(['Establishing geometry'])
aVec=[1:stimsize]';

stimbase.y=-repmat(aVec-maxrad,1,stimsize);
stimbase.x=repmat(aVec'-maxrad,stimsize,1);
stimbase.r=sqrt(stimbase.x.^2+stimbase.y.^2);
stimbase.th=atan2d(stimbase.y,stimbase.x)+90;
stimbase.th(stimbase.th<0)=stimbase.th(stimbase.th<0)+360;
stimbase.th=-stimbase.th+360;

% make a list of radii, to keep the aspect ratio for each annulus
% correct
for i=1:99
    j=100-i;
    rad(j)=floor(rad(j+1)-pi*rad(j+1)/NumChecksRad);
end
clear i j

% Fix the angles to be 0<=theta<360
theta=Dtheta.*[1:NumChecksTH];


%% Create B/W stimuli
disp(['Generating checkerboards'])

% Create angle windows
wTh=bsxfun(@ge,stimbase.th(:),theta(1:(end-1))).*...
    bsxfun(@lt,stimbase.th(:),theta(2:end));
wTh=cat(2,wTh,stimbase.th(:)<theta(1));

% Create radius windows
wRo=bsxfun(@ge,stimbase.r(:),rad(1:2:end)).*...
        bsxfun(@lt,stimbase.r(:),rad(2:2:end));
    
wRe=bsxfun(@ge,stimbase.r(:),rad(2:2:(end-1))).*...
        bsxfun(@lt,stimbase.r(:),rad(3:2:(end-1)));

stimA=zeros(stimsize*stimsize,1);
for th=1:NumChecksTH
    stimA=stimA+(sum(bsxfun(@times,wRo,wTh(:,th)),2)-...
        sum(bsxfun(@times,wRe,wTh(:,th)),2)).*((-1).^th); 
end
stimA=reshape(stimA,stimsize,stimsize);

stimA=(stimA+1).*0.5; % make values [0, 0.5, 1].

stimbase.On=stimA; 
stimbase.Off=-stimA+1; 


%% View Stimuli base
disp(['Visualizing'])
fields=fieldnames(stimbase);
figure;set(gcf,'Color','w')
for i=1:6
    subplot(2,3,i);   
    if any(strcmp(fields{i},{'On','Off'}))
        imshow(stimbase.([fields{i}]))
    else
    imagesc(stimbase.([fields{i}]))
    end
    axis image
    title([fields{i}])
end

%% Organize output with parameters
stimbase.rad=rad;
stimbase.theta=theta;
stimbase.stimsize=stimsize;
stimbase.NumChecksRad=NumChecksRad;
stimbase.NumChecksTH=NumChecksTH;