function Plot_Grid_Design(grid,Rad)
%
% This function generates a plot of a 2D and 3D representation of a grid
% design.
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



figure('Color','w','Units','Normalized','Position',[.1,.1,.8,.7]);
subplot(2,2,1);
PlotSD(grid.spos,grid.dpos,'norm',gcf)
title([{['2D grid: Ns = ',num2str(grid.srcnum),', Nd = ',...
    num2str(grid.detnum)]};{[' ']}])
subplot(2,2,2);
histogram(Rad.r2d,[1:0.1:(ceil(max(Rad.r2d)+1))])
xlabel('r_s_d [mm]');ylabel('N measurements')
title(['Closest r_s_d: ',num2str(min(Rad.r2d)),' mm'])
subplot(2,2,3)
PlotSD(grid.spos3,grid.dpos3,'norm',gcf);view([64,46]);
% title([{['3D representation of grid']};...
%     {['Radius of curvature: ',num2str(grid.Srad),' mm']}])
subplot(2,2,4);
histogram(Rad.r3d,[1:0.1:(ceil(max(Rad.r3d)+1))])
xlabel('r_s_d [mm]');ylabel('N measurements')
