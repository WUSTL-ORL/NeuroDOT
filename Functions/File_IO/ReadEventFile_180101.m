function Etimes=ReadEventFile_180101(fn)
%
% This function loads in an Event file, generated in PsychToolBox3 stimulus
% presentations, and outputs timing information into the structure Etimes.
% These files are used to synchronize MRI acquisition with the stimuli
% used. The first Event (Initialize) and the Initial_Fixation both occur
% before starting the BOLD sequence. The first BOLD frame should align with
% TR_wait_num_1 (and the time for this frame is set to zero).
% Etimes.synchpts are times (in seconds) of various events during the run.
% Etimes.TRidx are the indices of Events corresponding to BOLD frames.
% Etimes.Pulse_1 is the first BOLD frame (with time set to zero)
% Etimes.Pulse_2 is the first stimulus-type occurance
% Etimes.Pulse_3 is the second stim type, etc (may be fixation)
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

%% Load the Event file
[Event_num,Comp_time,Elapsed_time,dt,Label]=textread(fn,...
    '%d %f %f %f %s',-1,'headerlines',1);

%% Adjust times
TR1=find(strcmp(Label,'TR_wait_num_1'),1);

Elapsed_time=Elapsed_time(TR1:end);
Elapsed_time=Elapsed_time-Elapsed_time(1);
Label=Label(TR1:end);
Etimes.synchpts=Elapsed_time;


%% Label MRI BOLD frame times
TRlog=(strncmp(Label,'TR_wait',7)|strcmp(Label,'t'));
Etimes.TRidx=find(TRlog);

%% Label stimulation happenings
Event_types=unique(Label(~TRlog),'stable');
Event_types=Event_types(~strncmp(Event_types,'Initial',7));
Etimes.Pulse_1=Etimes.TRidx(1);
Etimes.Pulse_1_Label='Initial_Frame';

% List Fixation Last
if any(strncmp(Event_types,'Fixation',8))
    toMove=Event_types(strncmp(Event_types,'Fixation',8));
    Event_types(strncmp(Event_types,'Fixation',8))=[];
    Event_types=cat(1,Event_types,toMove);
end

% Label and list event types and their Pulse numbers
Nevents=size(Event_types,1);
for j=1:Nevents
    Etimes.(['Pulse_',num2str(j+1),'_Label'])=Event_types{j};
    Etimes.(['Pulse_',num2str(j+1)])=find(strcmp(Label,Event_types(j)));
end