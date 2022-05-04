function mesh=PrepareMeshForNIRFAST(mesh,meshname,tpos3)

% This function writes to disc the necessary *.meas, *.source, and *.link
% files needed for NIRFAST to solve the fwd light problem.
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


%% Check NIRFAST-style files are present, else make them
if ~exist([meshname,'.node'],'file')
%     mesh.source.fixed = 1;
    save_mesh(mesh,meshname);
end


%% Make .meas file
% This is just a dummy measurement, since we don't use NIRFAST's detectors
fid=fopen([meshname,'.meas'],'w');
fprintf(fid,'%f %f %f\n',mean(mesh.nodes(:,1)),...
    mean(mesh.nodes(:,2)),mean(mesh.nodes(:,3)));
fclose(fid);


%% Make .link file
% This is also just a placeholder since everything links to the one "detector"
fid=fopen([meshname,'.link'],'w');
for i=1:size(tpos3,1)
    fprintf(fid,'%f\n',1);
end
fclose(fid);


%% Make Source/Detector File
% Here in the .source file, we list all sources and detectors.
% This way, we always use NIRFAST's source Green's function, which use a
% better treatment of the boundary conditions.
% Later we will unwrap source and detector Green's functions.
% Structure: [xpos ypos zpos]
fid=fopen([meshname,'.source'],'w');
% fprintf(fid,'%s\n','fixed');
for j=1:size(tpos3,1)
    fprintf(fid,'%f %f %f\n',tpos3(j,1),tpos3(j,2),tpos3(j,3));
end
fclose(fid);


%% Load in full mesh
mesh=load_mesh([meshname]);
mesh.source.fixed = 1;