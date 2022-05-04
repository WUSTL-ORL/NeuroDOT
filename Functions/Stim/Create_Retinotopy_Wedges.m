function stim=Create_Retinotopy_Wedges(stimbase,NumChecks)

% This function generates checkerboard wedges using the stimbase structure
% and the desired number of checks as the width of the wedges.
%
%
%
%% Parameters and Initialization
% NumChecks=3;
Nframes=stimbase.NumChecksTH;
rad=stimbase.rad;           
theta=stimbase.theta;
radCH=84;
radOuter=99;
radCenter=87;

stim=struct;
stim.blank=stimbase.Blank;


%% Make crosshair
crossv=(abs(stimbase.x)<=1).*(abs(stimbase.y)<=rad(radCH));
crossh=(abs(stimbase.y)<=1).*(abs(stimbase.x)<=rad(radCH));
cross=(crossv+crossh)~=0;


%% Generate wedges
stim.blank=uint8((stim.blank+cross).*255);

for j=1:Nframes    
    minangle=theta(j)-theta(1);% Make Wedge
    maxangle=(j+NumChecks-1)*theta(1);    
    
    On=stimbase.On;        % initialize
    Off=stimbase.Off;
    
    ii=(stimbase.th>=minangle & stimbase.th<maxangle);
    if maxangle>=360
        ii=(ii+(stimbase.th<(maxangle-360)))~=0;
    end
        
    On(~ii)=0.5;  % Gray out the non-nonstim pixels
    Off(~ii)=0.5;
    
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
    stim.((['stim',num2str(j),'a']))=uint8(On.*255);
    stim.((['stim',num2str(j),'b']))=uint8(Off.*255);
end
