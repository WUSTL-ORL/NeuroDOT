function tpos3=Bend_Cap_Flat_to_Curved(tpos,hrad,type,aX)
%
% This function bends a 2D cap (tpos) around a given single radius 
% of curvature (hrad), assuming either spherical symmetry (type=='hem') or
% cylindrical symmetry (type=='cyl'). If 'cyl', must also include
% the preferred axis (aX) for the cylinder ('x',or,'y')
% The origin is placed at the center of mass of the tpos points.
%
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

%%
if ~exist('aX','var'),aX='y';end
Npos=size(tpos,1);
CofM=mean(tpos,1);
tpos3=zeros(Npos,3);


switch type
    case 'hem'
        for n=1:Npos
            x=tpos(n,1)-CofM(1);
            y=tpos(n,2)-CofM(2);
            r0=sqrt(x^2+y^2); % Distance from Origin
            theta=r0/hrad;
            phi=atan2(y,x);
            tpos3(n,:)=hrad*[sin(theta)*cos(phi),...
                sin(theta)*sin(phi),cos(theta)]+[CofM,0];
        end
        
    case 'cyl'
        switch aX
            case 'y'
        for n=1:Npos
            x=tpos(n,1)-CofM(1);
            y=tpos(n,2)-CofM(2);
            dTheta=x/hrad;
            X=x*cos(dTheta/2);
            Z=x*sin(dTheta/2);
            tpos3(n,:)=[X,y,-Z]+[CofM,0];
        end
            case 'x'
        for n=1:Npos
            y=tpos(n,1)-CofM(1);
            x=tpos(n,2)-CofM(2);
            dTheta=x/hrad;
            X=x*cos(dTheta/2);
            Z=x*sin(dTheta/2);
            tpos3(n,:)=[X,y,-Z]+[CofM,0];
        end
        end
end