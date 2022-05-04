function [grid,Rad]=Make_Sphere_Cap(Ndy,Ndx,NN1sep,Srad,type)
% 
%
% This function generates a flat cap with Ndy detectors in vertical and Ndx
% detectors in the horizonal direction. Sources columns lie both between
% and outside of detector columns (as with classic vis cap). NN1sep is the
% 1st nearest neighbor separation. CapName is given as a moniker in the
% grid and radius filenames.
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



%% Parameters and Initialization
if ~exist('type','var'), type='hem';end
Nsy=Ndy-1;
Nsx=Ndx+1;
ddx=NN1sep*sqrt(2); % Like-column spacing


%% Make 2D grid topology
dx=zeros(Ndx*Ndy,1);
dy=zeros(Ndx*Ndy,1);
sx=zeros(Nsx*Nsy,1);
sy=zeros(Nsx*Nsy,1);
n=0;
for j=1:Ndx
    for k=1:Ndy
        n=n+1;
        dy(n)=k;
        dx(n)=j;
    end
end
n=0;
for j=1:Nsx
    for k=1:Nsy
        n=n+1;
        sy(n)=k;
        sx(n)=j;
    end
end
dpos=[dx,dy];
spos=[sx,sy];
detnum=length(dpos);
srcnum=length(spos);
spos=bsxfun(@minus,spos,mean(spos));
dpos=bsxfun(@minus,dpos,mean(dpos));


%% Scale grid
dpos=dpos.*ddx; % Dsd in mm
spos=spos.*ddx;

grid.detnum=detnum;
grid.srcnum=srcnum;
grid.dpos=dpos;
grid.spos=spos;
tpos=cat(1,grid.spos,grid.dpos);

Rad=Grid2Radius_180824(grid,1.5);
Rad.rU2d=unique(Rad.r);
Rad.r2d=Rad.r;


%% Add spherical curvature
tpos3=Bend_Cap_Flat_to_Curved(tpos,Srad,type);
grid.spos3=tpos3(1:grid.srcnum,:);
grid.dpos3=tpos3((grid.srcnum+1):end,:);
n=0;
for j=1:grid.detnum
    for k=1:grid.srcnum
        n=n+1;
        Rad.r(n)=norm(grid.spos3(k,:)-grid.dpos3(j,:));
    end
end
grid.Srad=Srad;
grid.spos2=grid.spos;
grid.dpos2=grid.dpos;