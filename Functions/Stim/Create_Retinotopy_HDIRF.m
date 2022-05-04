function stim=Create_Retinotopy_HDIRF(stimbase)
%
% This function generates a set of retinotopy annuli from the stimbase
% structure.


%% Parameters and Initialization
% NcheckRadout=12;
% NcheckRad=2;
rad=stimbase.rad;           
radCH=84;
radCenter=87;
radOuter=99;

stim=struct;
stim.blank=stimbase.Blank;

%% Make crosshair
crossv=(abs(stimbase.x)<=1).*(abs(stimbase.y)<=rad(radCH));
crossh=(abs(stimbase.y)<=1).*(abs(stimbase.x)<=rad(radCH));
cross=(crossv+crossh)~=0;

stim.blank=uint8((stim.blank+cross).*255);

% Initialize
On=stimbase.On;
Off=stimbase.Off;

% Clear Outer Ring
outer=stimbase.r>rad(radOuter);
On(outer)=0.5;  % Gray out the non-nonstim pixels
Off(outer)=0.5;

% Clear out the center
center=(stimbase.r<rad(radCenter));
On(center)=0.5;  % Gray out the non-nonstim pixels
Off(center)=0.5;

% Add cross
On=On+cross;
Off=Off+cross;

% Put in Grid
stim.((['stim',num2str(1),'a']))=uint8(On.*255);
stim.((['stim',num2str(1),'b']))=uint8(Off.*255);