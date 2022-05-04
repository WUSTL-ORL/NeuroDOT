function [f,v] = vol2isomesh(volume, dim, isovalue)

% VOL2ISOMESH Interpolates volumetric data onto a isosurface mesh.
% 
%   [f,v] = VOL2ISOMESH(volume, dim, params) creates an isosurface given 
%   the isovalue for data within the volume and described by the volumetric
%   metadata in dim. If no isovalue is provided, a value of 1 will be used.
%   This function corrects for often confusing spatial arrangemet of
%   volumetric data.
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
nVx = dim.nVx;
nVy = dim.nVy;
nVz = dim.nVz;
dr = dim.mmppix;
center = dim.center;

if ~exist('isovalue', 'var')  
    isovalue = 1;
end

%% Define coordinate space of volumetric data
x = (dr(1).*[nVx:-1:1]-center(1))';
y = (dr(2).*[nVy:-1:1]-center(2))';
z = (dr(3).*[nVz:-1:1]-center(3))';

[X,Y,Z]=meshgrid(y,x,z);


%% Calculate isosurface and correct for space coordinates
[f,v]=isosurface(X,Y,Z,volume,isovalue);
v=v(:,[2,1,3]);
