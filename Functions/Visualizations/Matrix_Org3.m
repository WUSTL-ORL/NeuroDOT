function Matrix_Org3(Min,key,buffer,lims,Cmap,fig,CmapData)

% This function organizes a sorted matrix using lines to add visual
% separation between modules and block colors along the left and bottom.
% This 3rd version of Matrix_Org is designed to support rectangular input
% matrices. It will also work with square matrices. If the matrix is
% square, the code essentially works the same as Matrix_Org2. If the
% matrix is rectangular, then key and Cmap must be 2x1 cell arrays, where
% each cell contains the normal key/colormap. Cell 1 is the vertical and
% Cell 2 is the horizontal, which is assumed to be the larger dimension.
% Key is assumed to be organized with col1 as indices of Min and col2 as
% module assignment.
% Rectangles are drawn to denote each module, with colormap for modules
% passed in Cmap.
% Buffer is width of rectangles relative to matrix size
% lims sets limits for imagesc matrix display (current colormap is defined
% as jet below in 'Draw Matrix' section.
% Cmap is the colormap (passed, defined explicitly) for the rectangles.
% fig is a flag that creates a new figure if set to 1.
% CmapData is a colomap for the data itself, else default is jet.


%% Set Parameters
if ~exist('fig','var'), fig=1; end
if ~exist('Cmap','var'), Cmap=colorcube(max(key(:,2))); end
if fig==1, figure('Color','k'); end
if ~exist('CmapData','var'),CmapData=jet(1000); end
Mx=size(Min,2);
My=size(Min,1);
if Mx==My, type='Square';else type='Rect';end

switch type
    case 'Square'
        Nets=unique(key(:,2));
        Nrect=length(Nets);
    case 'Rect'
        NetsV=unique(key{1}(:,2));
        NrectV=length(NetsV);
        NetsH=unique(key{2}(:,2));
        NrectH=length(NetsH);
        if ~iscell(Cmap)
            foo=Cmap;
            Cmap=cell(2,1);
            Cmap{1}=foo;Cmap{2}=foo;
        end
end

%% Draw matrix
imagesc(Min,[lims(1) lims(2)]);
colormap(CmapData);hold on
axis([-buffer+0.5,Mx+0.5,0.5,My+buffer+0.5])


%% Draw lines
switch type
    case 'Square'
        for j=1:(Nrect-1)
            if any(find(key(:,2)==Nets(j)))
                % Horizontal
                plot([0.5,Mx+0.5],...
                    [find(key(:,2)==Nets(j),1,'last')+0.5,...
                    find(key(:,2)==Nets(j),1,'last')+0.5],'k');
                % Vertical
                plot([find(key(:,2)==Nets(j),1,'last')+0.5,...
                    find(key(:,2)==Nets(j),1,'last')+0.5],...
                    [0.5,Mx+0.5],'k');
            end
        end
        
    case 'Rect'        
        % Horizontal
        for j=1:(NrectV-1)
            if any(find(key{1}(:,2)==NetsV(j)))
                plot([0.5,Mx+0.5],...
                    [find(key{1}(:,2)==NetsV(j),1,'last')+0.5,...
                    find(key{1}(:,2)==NetsV(j),1,'last')+0.5],'k');
            end
        end
            % Vertical
        for j=1:(NrectH-1)
            if any(find(key{2}(:,2)==NetsH(j)))
                plot([find(key{2}(:,2)==NetsH(j),1,'last')+0.5,...
                    find(key{2}(:,2)==NetsH(j),1,'last')+0.5],...
                    [0.5,Mx+0.5],'k');
            end
        end
end


%% Draw Rectangles
switch type
    case 'Square'
        
for j=1:Nrect
    if any(find(key(:,2)==Nets(j)))
        % Y-axis rectangles
        x=-buffer+0.5;
        y=find(key(:,2)==Nets(j),1,'first')-0.5;
        w=buffer;
        h=find(key(:,2)==Nets(j),1,'last')-y+0.5;
        rectangle('Position',[x,y,w,h],'FaceColor',Cmap(Nets(j),:))
        % X-axis rectangles
        x=find(key(:,2)==Nets(j),1,'first')-0.5;
        y=My+0.5;
        w=find(key(:,2)==Nets(j),1,'last')-x+0.5;
        h=buffer;
        rectangle('Position',[x,y,w,h],'FaceColor',Cmap(Nets(j),:))
    end
end

    case 'Rect'
        % Y-axis rectangles
    for j=1:NrectV
        if any(find(key{1}(:,2)==NetsV(j)))
            x=-buffer+0.5;
            y=find(key{1}(:,2)==NetsV(j),1,'first')-0.5;
            w=buffer;
            h=find(key{1}(:,2)==NetsV(j),1,'last')-y+0.5;
            rectangle('Position',[x,y,w,h],'FaceColor',Cmap{1}(NetsV(j),:))
        end
    end
        
        % X-axis rectangles
    for j=1:NrectH
        if any(find(key{2}(:,2)==NetsH(j)))
            x=find(key{2}(:,2)==NetsH(j),1,'first')-0.5;
            y=My+0.5;
            w=find(key{2}(:,2)==NetsH(j),1,'last')-x+0.5;
            h=buffer;
            rectangle('Position',[x,y,w,h],'FaceColor',Cmap{2}(NetsH(j),:))
        end
    end
end

%% If data is lower-triangular, add white
if (length(unique(triu(Min,1)))==1)
% Add white patch with white edges
vX=reshape(repmat([(1+0.5):(Mx+0.5),(Mx+0.5)],[2,1]),[],1);
vX(end-1:end)=[];
vY=reshape(repmat([(0+0.5):(My-0.5),(0+0.5)],[2,1]),[],1);
vY(1)=[];vY(end)=[];
v=[vX,vY];
f=[1:length(v),1];
% patch('Faces',f,'Vertices',v,'FaceColor','k','EdgeColor','k');
patch('Faces',f,'Vertices',v,'FaceColor','w','EdgeColor','w');

% Fix black edges on perimeter
vX=reshape(repmat([(1+0.5):(Mx+0.5)],[2,1]),[],1);
vX(end)=[];
vY=reshape(repmat([(0+0.5):(My-0.5)],[2,1]),[],1);
vY(1)=[];
v=[vX,vY;vX(end),(vY(end)+buffer+1);...
    (vX(1)-buffer-1),(vY(end)+buffer+1);(vX(1)-buffer-1),(vY(1))];
f=[1:length(v),1];
patch('Faces',f,'Vertices',v,'FaceColor','none','EdgeColor','k');
set(gca,'Box','off','XTick',[],'YTick',[],'DataAspectRatio',[1,1,1],...
    'TickDir','out')
else
v=[(Mx+0.5),(My+0.5+buffer);(0.5-buffer),(My+0.5+buffer);...
    (0.5-buffer),0.5;(Mx+0.5),0.5];
f=[1:length(v),1];  
patch('Faces',f,'Vertices',v,'FaceColor','none','EdgeColor','k');    
    
set(gca,'Box','off','XTick',[],'YTick',[],'DataAspectRatio',[1,1,1],...
    'TickDir','out')
end

%% Wrap up
set(gca,'Box','off','XTick',[],'YTick',[],'DataAspectRatio',[1,1,1],...
    'TickDir','out')
freezeColors;