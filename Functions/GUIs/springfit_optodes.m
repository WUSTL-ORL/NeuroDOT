function tpos = springfit_optodes(mesh, rad, spos3, dpos3, anchorpts, anchor)

% SPRINGFIT_OPTODES Performs the springfit algorithm on a cap grid.
% 
%   tpos = SPRINGFIT_OPTODES(mesh, rad, spos3, dpos3, anchorpts, anchor)
%   applies an ORL-developed algorithm that fits a cap grid onto a mesh by
%   minimizing a spring energy calculation.
% 
% See Also: CAP_FITTER.
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

%% Parameters and Initialization.
k.h = 1; % Head Spring constant
k.a = 1.5; % Anchor Spring constant
k.nn = 1.5; %1; % Nearest Neighbor Spring constant
psize0 = 10 ^- 1; % Step size (for head)
loopmax = 10; % Max tries to half step size before giving up on that optode.
ethresh = 1; % Energy minimum threshold
dthresh = 0.1; % Difference threshold
compnum = 2; % Back-iterations to compare

info.nn1 = sort(rad.nn1);   % Set of 1st nn in meas
info.nn2 = sort(rad.nn2);   % Set of 2nd nn in meas
info.nn12 = sort(rad.nn12);   % Set of 2nd nn in meas
info.meas = rad.meas; % S-D pairs
info.r = rad.r;
rad.srcnum=max(rad.meas(:,1));
rad.detnum=max(rad.meas(:,2));

% load(['grid_',padname],'srcnum','detnum')
info.meas(:,2)=info.meas(:,2)+rad.srcnum; % relative to optode not S,D
opnum=rad.srcnum+rad.detnum;
info.srcnum=rad.srcnum;
info.detnum=rad.detnum;
tpos=[spos3;dpos3];
% tpos=[[tpos(:,1)],[tpos(:,2)],[tpos(:,1)]];
% Move CM back to anchrs + 25
% MuAnch=mean(anchor,1);
% MuCap1=mean([tpos(anchorpts',:)],1);
% deltaAnch=MuAnch-MuCap1;
% for i=1:3
%     tpos(:,i)=tpos(:,i)+deltaAnch(i);
% end
% tpos(:,1)=tpos(:,1)-35; % Move behind head
% % tpos(:,2)=tpos(:,2)+5; % Move up on head
% tpos=1.05*tpos;

%% Initialize spring fit.
% Prep Mesh
TR=TriRep(mesh.elements,mesh.nodes);% Make MATLAB triangulation class
[elemb nodesb]=freeBoundary(TR);    % Get Boundary
TRb=TriRep(elemb,nodesb);           % triangulation of boundary
faceb=incenters(TRb);               % centers of faces
[~,MfInd]=max(faceb(:,2));    % furthest face in Y
tempNorm=faceNormals(TRb);          % Get norms then test to chose direction
if tempNorm(MfInd,2)<0, fn=-tempNorm;
else fn=tempNorm;
end

% Initial Energy
[nnb,surfsep]=knnsearch(faceb,tpos); % initial separation of all optodes from surface

f1=figure;
plottpos(nodesb,elemb,tpos,anchor,anchorpts,info,1,f1); % plot initial position

% Calculate energy of each optode individually
E0=0;
for i=1:opnum
    % Find nn optodes, not measurements for given optode
    MeasOp=[find(info.meas(:,1)==i);find(info.meas(:,2)==i)];
    OnSpring=setxor(unique(info.meas(intersect(MeasOp,info.nn12),:)),i);
    
    E(i)=springerror(tpos,i,OnSpring,anchor,anchorpts,...
        info,k,surfsep(i)); % initial energy
    E0=E0+E(i);
end

loopnum=0;
looptry=0;
onhead=[];
Eiter=E0;

% Check to see how close we are to desired config.
targetcheck(info,tpos,loopnum,looptry,loopmax,E0,onhead,surfsep);

% pause

%% Minimize energy of spring configuration
% matlabpool open 8
while 1  % do til converge
    loopnum=loopnum+1; % count iterations
    E0=Eiter(end);
    
    % Move each optode in parallel in comparison with previous config.
    % (1) Find surface position and normal
    % (2) Find gradient to lower energy
    % (3) Move optode distance towards head (ideally)
    % (4) Calc change in E to keep or discard movement
    
    tic
    parfor op=1:opnum
        % (1)
        %         psizeOP(op)=psize0;
        surfpos=faceb(nnb(op),:);
        %         surfN=fn(nnb(op),:);
        
        % (2)
        MeasOp=[find(info.meas(:,1)==op);find(info.meas(:,2)==op)]; %#ok<*PFBNS>
        OnSpring=setxor(unique(info.meas(intersect(MeasOp,info.nn12),:)),op);
        
        pgrad=Egrad(tpos,op,OnSpring,info,anchorpts,anchor,k,surfpos);
        
        % (3) & (4)
        looptry(op)=0;
        tpos2=tpos;
        while 1
            looptry(op)=looptry(op)+1;
            
            if looptry(op)>loopmax % give up
                break
            end
            
            % (3) Here is the new optode position
            tpos2(op,:)=tpos(op,:)+psize0*pgrad;
            
            % Head issues pending... if optode is inside head
            
            % (4)Energy of new position
            [~,surfsep]=knnsearch(faceb,tpos2);
            Etry=springerror(tpos2,op,OnSpring,anchor,anchorpts,info,k,surfsep(op));
            
            % Check to see if the iteration is good
            if Etry<E(op) % if we reduced the energy
                E(op)=Etry; % This is the new current energy
                tpos1(op,:)=tpos2(op,:);
                break
            else % try again with smaller step
                %                 psizeOP(op)=0.5*psizeOP(op);
                tpos1(op,:)=tpos(op,:);
                continue
            end
        end
        
    end
    
    tpos=tpos1;
    
    Eiter(loopnum)=sum(E);
    [nnb,surfsep]=knnsearch(faceb,tpos); % Find the closest points on the surface for each optode.
    OnHead(loopnum)=length(find(surfsep<1));
    [MaxTries(loopnum),NN1ave(loopnum),NN2ave(loopnum),ADFH(loopnum)]=...
        targetcheck(info,tpos,loopnum,looptry,loopmax,Eiter(loopnum),OnHead(loopnum),surfsep);
    t(loopnum)=toc;
    disp([' Elapsed time for loop = ',num2str(t(loopnum)),' seconds.'])
    
    
    % Plot optode positions every 2 iterations
    if ~mod(loopnum,2)
        plottpos(nodesb,elemb,tpos,anchor,anchorpts,info,0,f1);
    end
    
    %%% Reasons for Terminating Loop %%%
    
    % (1) Energy is below threshold (at minimum).  (This won't happen.)
    if E0<ethresh
        disp('Energy below threshold: terminating loop')
        break
    end
    % (2)Difference between loops is small. (This is what will happen, either
    % for good reasons, or because it will reject so many changes that it stalls).
    if loopnum>compnum
        if abs(E0-Eiter(loopnum-compnum))<dthresh
            disp('Energy difference below threshold: terminating loop')
            break
        end
    end
    % (3) If all optodes failed movement, then just stop instead of waiting for
    % Delta-E to be low.
    if numel(find(looptry>loopmax))==opnum
        disp('No optodes could move: terminating loop')
        break
    end
    
    % Plot values for loops
    if loopnum==1, f2=figure;end
    set(0,'CurrentFigure',f2)
    subplot(2,2,1); hold off;   % Optodes with max tries and on head
    plot(MaxTries,'-ro','LineWidth',2);hold on;
    plot(OnHead,'-bo','LineWidth',2)
    title('Max Tries (r) and on head (b)')
    ylabel('Number of Optodes')
    xlabel('Loop number')
    xlim([0 loopnum])
    subplot(2,2,2);             % Energy
    semilogy(Eiter./1e6,'-ko','LineWidth',2);
    title('Energy of configuration')
    ylabel('Joules')
    xlabel('Loop number')
    xlim([0 loopnum])
    subplot(2,2,3); hold off;   % Average NN distances, and 10* Average distance from head
    semilogy(NN1ave,'-ro','LineWidth',2);hold on;
    semilogy(NN2ave,'-bo','LineWidth',2);
    semilogy(10.*ADFH,'-go','LineWidth',2)
    title('NN1ave (r), NN2ave (b), 10* ADFH (g)')
    ylabel('Millimeters')
    xlabel('Loop number')
    xlim([0 loopnum])
    subplot(2,2,4);             % Time per loop
    plot(t,'-ko','LineWidth',2);
    title('Time per loop')
    ylabel('Seconds')
    xlabel('Loop number')
    xlim([0 loopnum])
    
    pause(0.1)
    
end

% matlabpool close


disp(['Total time = ',num2str(sum(t))])
end


%% springerror(): energy function to be minimized
function E=springerror(tpos,op,OnSpring,anchor,anchorpts,info,k,surfsep)

% Etot=Eanchor+Enn1+Esurf

% Energy from anchors
Eanchor=0;
if intersect(op,anchorpts)
    ind=find(anchorpts==op);
    Eanchor=Eanchor+(1/2)*k.a*norm(tpos(anchorpts(ind),:)-anchor(ind,:))^2;
end

% Energy from nn
Enn=0;
for m=1:length(OnSpring) % loop through measurements
    sep=norm(tpos(OnSpring(m),:)-tpos(op,:));
    if op>info.srcnum
        MeasM=intersect(find(info.meas(:,1)==OnSpring(m)),find(info.meas(:,2)==op));
    else
        MeasM=intersect(find(info.meas(:,2)==OnSpring(m)),find(info.meas(:,1)==op));
    end
    Enn=Enn+(1/2)*k.nn*(sep-info.r(MeasM))^2;
end

% Energy from surface
Esurface=(1/2)*k.h*surfsep^2;

% Add all contributions into one energy
E=Eanchor+Enn+Esurface;

end


%% Egrad(): gradient of energy function
function [gradout]=Egrad(tpos,op,OnSpring,info,anchorpts,anchor,k,surfpos)

% Initialize
gradout=[0 0 0];

% Grid Spring Forces
for m=1:length(OnSpring)
    dr=tpos(op,:)-tpos(OnSpring(m),:);
    sep=norm(dr);
    if op>info.srcnum
        MeasM=intersect(find(info.meas(:,1)==OnSpring(m)),find(info.meas(:,2)==op));
    else
        MeasM=intersect(find(info.meas(:,2)==OnSpring(m)),find(info.meas(:,1)==op));
    end
    gradout=gradout-k.nn*(1-info.r(MeasM)/sep)*dr;
end

% Anchor Springs
if intersect(op,anchorpts)
    ind=find(anchorpts==op);
    gradout=gradout-k.a*(tpos(anchorpts(ind),:)-anchor(ind,:));
end

% Surface Spring
gradout=gradout-k.h*(tpos(op,:)-surfpos);

end


%% targetcheck()
function [MaxTries,NN1ave,NN2ave,ADFH]=...
    targetcheck(info,tpos,loopnum,looptry,loopmax,E0,onhead,surfsep)

% How far are we from desired optode separations?
m1=0; % reset errors
m2=0;
for m=1:info.srcnum*info.detnum
    s=info.meas(m,1);
    d=info.meas(m,2);
    sep=norm(tpos(s,:)-tpos(d,:));
    
    if ismember(m,info.nn1)
        m1=m1+sep;
    end
    
    if ismember(m,info.nn2)
        m2=m2+sep;
    end
end
m1=m1/numel(info.nn1); % Average
m2=m2/numel(info.nn2);


MaxTries=numel(find(looptry>loopmax));
NN1ave=m1;
NN2ave=m2;
ADFH=mean(surfsep);


disp(['Loop: ',num2str(loopnum)])
disp([' Optodes w/ Max Tries: ',num2str(MaxTries)])
disp([' Optodes on Head: ',num2str(onhead)])
disp([' E: ',num2str(E0)])
disp([' Ave NN1 Distances (mm): ',num2str(m1)])
disp([' Ave NN2 Distances (mm): ',num2str(m2)])
disp([' Ave Distance from head (mm): ',num2str(ADFH)])

end


%% plotpos()
function plottpos(nodesb,elemb,tpos,anchor,anchorpts,info,m,f1)

set(0,'CurrentFigure',f1)
if m==1
    
    % Mesh
    % plot3(nodes(:,1),nodes(:,2),nodes(:,3),'c.')
    grayLevel=0.5;
    set(gcf,'Color',[1 1 1]);
    h=patch('Faces',elemb,'Vertices',nodesb,'EdgeColor','black',...
        'FaceColor',[grayLevel grayLevel grayLevel]);
    hold on;
    
    % Lower lighting
    light('Position',[-140,90,-100],'Style','local ');
    light('Position',[-140,-350,-100],'Style','local ');
    light('Position',[300,90,-100],'Style','local ');
    light('Position',[300,-350,-100],'Style','local ');
    % light('Position',[80,0,-500],'Style','local ');
    
    % Higher lighting
    light('Position',[-140,90,360],'Style','local ');
    light('Position',[-140,-350,360],'Style','local ');
    light('Position',[300,90,360],'Style','local ');
    light('Position',[300,-350,360],'Style','local ');
    
    set(h,'FaceLighting','flat','AmbientStrength',0.25);
    
    % Anchors
    [x,y,z]=sphere(100);
    rad=4;
    for i=1:size(anchor,1)
        hh=patch(surf2patch(rad*x+anchor(i,1),rad*y+anchor(i,2),rad*z+anchor(i,3)),...
            'EdgeColor','green','FaceColor','green','EdgeAlpha',0);
        set(hh,'FaceLighting','phong','AmbientStrength',0.02);
    end
    
    axis image
    view(-153,-38)
    xlabel('x')
    ylabel('y')
    zlabel('z')
end



hold on;

% Sources
for i=1:info.srcnum
    text(tpos(i,1),tpos(i,2),tpos(i,3),num2str(i),'FontSize',10,'Color','r')
end

% Detectors
for i=info.srcnum+1:(info.srcnum+info.detnum)
    text(tpos(i,1),tpos(i,2),tpos(i,3),num2str(i-info.srcnum),'FontSize',10,'Color','b')
end

% Anchorpoints
for i=1:length(anchorpts)
    text(tpos(anchorpts(i),1),tpos(anchorpts(i),2),tpos(anchorpts(i),3),...
        num2str(anchorpts(i)),'FontSize',10,'Color',[1 0 1])
end

% Rotate viewing angle by 30 degrees each iteration
shift=30;
RotMat=[cosd(shift),-sind(shift),0;sind(shift),cosd(shift),0;0,0,1];
view([get(gca,'CameraPosition')*RotMat])

% pause to allow visualization
pause(0.01)

end