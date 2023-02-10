function [A,dim,Gsd,Gs,Gd,dc]=makeAnirfaster(mesh,flags)

% This function (1) Calculates the green's functions, and 
%               (2) generates the A-matrix.  
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


%% (1) Parameters and Initialization
if ~isfield(flags,'makeA'),flags.makeA=1;end
if ~isfield(flags,'Hz'),flags.Hz=0;end
opnum=size(mesh.source.coord,1);
labels=flags.labels;
op=flags.op;
numc=length(op.lambda);

f=fieldnames(labels);
numf=numel(f);
region=mesh.region;
regtypes=unique(region);
numr=numel(regtypes);

if numr~=numf
    error(['** The number of regions in the mesh does not ',...
        'mesh the number of described regions **'])
end

% initialize arrays
a=zeros(size(region));
kappa=zeros(size(region));
ri=zeros(size(region));
numNodes=size(mesh.nodes,1);
G=zeros(numc,opnum,numNodes); % run only once 
dc=zeros(numc,numNodes);



%% (2) Choose wavelengths and get Gfuncts then save them
for lambda=1:2
    TimeWL=tic;
    for n=regtypes'
        roi=find(region==n);
        
        floop=['r',num2str(n)];
        
        if ~isfield(labels,floop)
            error(['** There is not description matching region number ',...
                num2str(n),' **'])
        elseif ~isfield(op,['mua_',labels.(floop)]) || ...
                ~isfield(op,['musp_',labels.(floop)]) || ...
                ~isfield(op,['n_',labels.(floop)])
            error(['** Not all optical properties are specified for region ',...
                labels.(floop),' **'])
        end
        
        a(roi)=op.(['mua_',labels.(floop)])(lambda);                % mua
        kappa(roi)=1/(3*(op.(['mua_',labels.(floop)])(lambda)+...   % D.C.
            op.(['musp_',labels.(floop)])(lambda)));
        ri(roi)=op.(['n_',labels.(floop)])(lambda);                 % Ind. of Ref.
    end
    
    mesh.mua=a;
    mesh.kappa=kappa;
    mesh.ri=ri;
    
    disp(['>Making Greens Functions, Wavelength ',num2str(lambda)])
    tic;[data]=femdata_stnd_FD(mesh,flags.Hz);toc;
    
    G(lambda,:,:)=data.phi';% G:wl-src-nodes
    dc(lambda,:)=mesh.kappa;
    toc(TimeWL)
end

if flags.Hz==0
    G(G<0)=0; % Fix negative values for CW case
end
Gs=single(G(:,1:flags.srcnum,:));
Gd=single(G(:,(flags.srcnum+1):end,:));

disp(['<< Saving Green''s Functions'])
save(['GFunc_',flags.tag,'.mat'],'G','Gs','Gd','dc','flags','-v7.3')

clear G


%% (3) Create A matrix from Green's Funcions, mesh, and optical props.
% if flags.makeA
[A,dim,Gsd,Gs,Gd,dc]= GtoAmat(Gs,Gd,mesh,dc,flags);
% else
%     A=[];
%     dim=[];
%     Gsd=[];
% end