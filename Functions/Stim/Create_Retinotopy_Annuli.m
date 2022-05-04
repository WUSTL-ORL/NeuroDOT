function stim=Create_Retinotopy_Annuli(stimbase,NcheckRadout,NcheckRad)
%
% This function generates a set of retinotopy annuli from the stimbase
% structure.


%% Parameters and Initialization
% NcheckRadout=12;
% NminradcheckRad=2;
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


%% Generate annuli
for j=1:NcheckRadout
    
    On=stimbase.On;
    Off=stimbase.Off;
    
    jr=max([radCenter+j-(NcheckRad-1),radCenter]);
    k=min([radCenter+j+1,radOuter]);
    minrad=rad(jr);
    maxrad=rad(k);
    
    ii=(stimbase.r>minrad).*(stimbase.r<=maxrad);
    
    On(~ii)=0.5;  % Gray out the non-nonstim pixels
    Off(~ii)=0.5;    
    
    % Add cross
    On=On+cross;
    Off=Off+cross;
    
    % Put in Grid
    stim.((['stim',num2str(j),'a']))=uint8(On.*255);
    stim.((['stim',num2str(j),'b']))=uint8(Off.*255);
end


%% Fix wrap-around
% if NcheckRad>1
%     for j=1:NcheckRad-1
% idx=find(stim.((['stim',num2str(1),'a']))~=...
%     stim.((['stim',num2str(j+NcheckRad-1),'a'])));
% 
% stim.(['stim',num2str(NcheckRadout-j+1),'a'])(idx)=...
%     stim.((['stim',num2str(1),'a']))(idx);
% stim.(['stim',num2str(NcheckRadout-j+1),'b'])(idx)=...
%     stim.((['stim',num2str(1),'b']))(idx);
%     end
% end