function PlotSD(spos,dpos,type,F,textColor,rad)

% This function plots the source and detector positions for a cap.
% If the user wants the positions placed onto a current figure, then
% include the handle F.
%
% spos is the set of coordinates for the sources.
% dpos is the set of coordinates for the detectors.
% type selects what kind of markers to use for SD.  'norm' shows stars and numbers.  'render' creates balls.
% F is the figure handle you want to place these SD markers on.
% textColor is an option if you don't want to use black text on SD labels.
% 
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

if ~exist('type','var'), type='norm';end
if exist('F','var'),set(0,'CurrentFigure',F);hold on;else figure;end
if ~exist('textColor','var'), textColor='k';end
if ~exist('rad','var'), rad=2;end

if strcmp(textColor,'w')==1
    set(gcf,'Color','k');set(gca,'Color','k')
end

switch type
    case 'norm'
        
if size(spos,2)==2
    plot(spos(:,1),spos(:,2),'r*');
    hold on
    for i=1:size(spos,1)
        text(spos(i,1),spos(i,2),num2str(i),'Color',textColor)
    end
    plot(dpos(:,1),dpos(:,2),'bo');
    for i=1:size(dpos,1)
        text(dpos(i,1),dpos(i,2),num2str(i),'Color',textColor)
    end
    axis image
    xlabel('X')
    ylabel('Y')
else
    plot3(spos(:,1),spos(:,2),spos(:,3),'r*');
    hold on
    for i=1:size(spos,1)
        text(spos(i,1),spos(i,2),spos(i,3),num2str(i),'Color',textColor,...
            'FontSize',20)
    end
    plot3(dpos(:,1),dpos(:,2),dpos(:,3),'bo');
    for i=1:size(dpos,1)
        text(dpos(i,1),dpos(i,2),dpos(i,3),num2str(i),'Color',textColor,...
            'FontSize',20)
    end
    axis image
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
end


    case 'discs'
MarkerSize=10;
plot(spos(:,1),spos(:,2),'or','MarkerEdgeColor','k','MarkerFaceColor',[1 0.5 0.5],...
    'MarkerSize',MarkerSize);
    hold on
plot(dpos(:,1),dpos(:,2),'ob','MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 1],...
    'MarkerSize',MarkerSize);


    case 'render'

% If present, plot sources and detectors and anchors
[x,y,z]=sphere(100);
% rad=2;
for i=1:size(spos,1) 
    hh=patch(surf2patch(rad*x+spos(i,1),rad*y+spos(i,2),rad*z+spos(i,3)),...
       'EdgeColor','red','FaceColor','red','EdgeAlpha',0);
%     set(hh,'FaceLighting','phong','AmbientStrength',0.02);
    set(hh,'FaceLighting','flat','AmbientStrength',0.02);
%     text(Scoord(i,1),Scoord(i,2),Scoord(i,3),num2str(i),'Color','k',...
%         'FontSize',20);
end

for i=1:size(dpos,1) 
    hh=patch(surf2patch(rad*x+dpos(i,1),rad*y+dpos(i,2),rad*z+dpos(i,3)),...
       'EdgeColor','blue','FaceColor','blue','EdgeAlpha',0);
%     set(hh,'FaceLighting','phong','AmbientStrength',0.02);
    set(hh,'FaceLighting','flat','AmbientStrength',0.02);
%     text(Dcoord(i,1),Dcoord(i,2),Dcoord(i,3),num2str(i),'Color','k',...
%         'FontSize',20);
end

    case 'text'
        TextSize=10;
if size(spos,2)==2
    plot(spos(:,1),spos(:,2),'k.');
    set(gcf,'Color','k')
    hold on
    for i=1:size(spos,1)
        text(spos(i,1),spos(i,2),num2str(i),'Color',[1 .75 .75],...
            'HorizontalAlignment','center','FontSize',TextSize,'FontWeight','bold')
    end
    for i=1:size(dpos,1)
        text(dpos(i,1),dpos(i,2),num2str(i),'Color',[.55 .55 1],...
            'HorizontalAlignment','center','FontSize',TextSize,'FontWeight','bold')
    end
    axis image
    xlabel('X')
    ylabel('Y')
else
    plot3(spos(:,1),spos(:,2),spos(:,3),'r*');
    hold on
    for i=1:size(spos,1)
        text(spos(i,1),spos(i,2),spos(i,3),num2str(i),'Color',textColor)
    end
    plot3(dpos(:,1),dpos(:,2),dpos(:,3),'bo');
    for i=1:size(dpos,1)
        text(dpos(i,1),dpos(i,2),dpos(i,3),num2str(i),'Color',textColor)
    end
    axis image
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
end
end