function DrawColoredSynchPoints(info,SPfr)
%
% This function draws vertical lines over a plot/imagesc type axis to
% delineate synchronization points on time traces.
% This function assumes the standard NeuroDOT info structure and paradigm
% format: 
%   info.paradigm       Contains all paradigm timing information
%   *.synchpts          Contains sample points of interest
%   *.synchtype         Optional field that contains marker that can
%                           distinguish between synchronization types,
%                           traditionally as sound boops of differing
%                           frequencies, i.e., 25, 30, 35 (Hz).
%   *.Pulse_1           Conventionally used to index synchpts that
%                           denote 'bookends' of data collection as
%                           well as resting (or OFF) periods of a given
%                           stimulus.
%   *.Pulse_2           Conventionally used to index synchpts that
%                           correspond to the start of a given stimulus
%                           epoch (e.g., the start of a flickering
%                           checkerboard)
%   *.Pulse_3           Index of synchpts of a 2nd stimulus type
%   *.Pulse_4           Index of synchpts of a 3rd stimulus type
% The 2nd input, SPfr, selects whether or not to adjust synchpoint timing
% based on framerate of data. (Default=0).


%% Parameters and Initialization
if ~isfield(info,'paradigm'),return;end
if ~isfield(info.paradigm,'synchpts'),return;end
h=gca;
xLim=get(h,'XLim');
yLim=get(h,'YLim');
hold on
if ~exist('SPfr','var'),SPfr=0;end
if SPfr
    fr=info.system.framerate;
    synchs=info.paradigm.synchpts./fr;
else
    synchs=info.paradigm.synchpts;
end
for j=1:4
    if ~isfield(info.paradigm,['Pulse_',num2str(j)])
        info.paradigm.(['Pulse_',num2str(j)])=[];
    end
end

%% Draw lines
for j=1:length(synchs)    % Draw synch pt bars
    if ismember(j,info.paradigm.Pulse_1)
        plot([1,1].*synchs(j),yLim,'r','LineWidth',1)
    elseif ismember(j,info.paradigm.Pulse_2)
        plot([1,1].*synchs(j),yLim,'g','LineWidth',1)
    elseif ismember(j,info.paradigm.Pulse_3)
        plot([1,1].*synchs(j),yLim,'b','LineWidth',1)
    elseif ismember(j,info.paradigm.Pulse_4)
        plot([1,1].*synchs(j),yLim,'m','LineWidth',1)
    else
        if any(h.Color)
            plot([1,1].*synchs(j),yLim,'k','LineWidth',1)   
        else
            plot([1,1].*synchs(j),yLim,'w','LineWidth',1)
        end
    end
end