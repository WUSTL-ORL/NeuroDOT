function [Area]=Active_Area_Surface_Mesh(nodes,elements,data)
%
% This function calculates the area of some active region in a surface
% triangular mesh. nodes: 3D positions of nodes. elements:  connectivity of
% mesh. data contains array of nodes to include. Threshold is currently set
% to zero. 
% Triangles of interest have points a, b, c


keep=find(data>0);
tf=ismember(elements,keep);
tf=sum(tf,2)==3; 
triags=elements(tf,:);

a=nodes(triags(:,1),:);         % vertices
b=nodes(triags(:,2),:);
c=nodes(triags(:,3),:);

A=sqrt(sum(((b-c).^2),2));      % lengths
B=sqrt(sum(((a-c).^2),2));
C=sqrt(sum(((b-a).^2),2));

P=0.5.*(A+B+C);                 % semiperimeter

Area=sum(sqrt(P.*(P-A).*(P-B).*(P-C))); % Total Area
