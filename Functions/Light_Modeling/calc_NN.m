function info_out = calc_NN(info_in, dr)
%
% This function takes in a neuroDOT info structure (info_in), calculates the
% Nearest Neigbor value for all measurement pairs, populates the output
% structure info_out.pairs.NN field with those values.

% The input dr provides a minimum separation for sources and detectors to be 
% grouped into different neighbors. dr defaults to 10. All distances are in
% millimeters


%% Parameters and Initialization

info_out = info_in;

if ~exist('dr', 'var')
    dr = 10; %default = 10mm minimum separation for SD to be grouped into different neighbors
end

Nm = length(info_in.pairs.r3d);

NN=zeros(Nm,1); %initialize NN vector


%% Calculate NN's
RadMaxR=ceil(max(info_in.pairs.r3d)); %maximum SD separation across all SD pairs


d=0; % s-d distance
c=0; % which nn are we on?
while any(info_in.pairs.r3d>d) % as long as there are still s-d pairs left to group
    if ((d==0) && (dr>9)),...
        nnkeep=find(info_in.pairs.r3d>=d & info_in.pairs.r3d<(d+(2*dr)));d=d+dr;
    else nnkeep=find(info_in.pairs.r3d>=d & info_in.pairs.r3d<(d+dr));
    end 
     % find pairs within 1 mm range
    if isempty(nnkeep) % if there are none
    else % if we find pairs
        c=c+1; % increment nn count
            NN(nnkeep)=c;
        if c>RadMaxR; break; end % stop at nn9
    end
    d=d+dr; % incremement search distance
end

info_out.pairs.NN=NN;
    
    
end