function [ s,w ] = run_freesurfer_command(FScmd, FSsetup)
% RUN_FREESURFER_COMMAND sets up the environment for calling FreeSurfer
% recon_all.
%
%
%
% Copyright (c) 2024 Washington University 
% Created By: Michael S. Jones
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameters and Initialization
% full path to FreeSurfer top level install dir
freesurferdir = FSsetup.freesurferdir;

% shell used to run FreeSurfer
freesurfershell = FSsetup.freesurfershell;

% full path to FreeSurfer setup script, run before any FS command
freesurfersetup = FSsetup.freesurfersetup;

% full path to FreeSurfer environmental setup script, run before any FS command
freesurferenvironment = FSsetup.freesurferenvironment;

% full path to FreeSurfer subjects directory
freesurfersubjectsdir = FSsetup.freesurfersubjectsdir;

% freesurfer output format ('nii' or 'nii.gz')
freesurferoutputformat = FSsetup.freesurferoutputformat;


% Check whether ${FSDIR}/bin is already in path

pth = getenv('PATH');

FSbin = fullfile(freesurferdir,'bin');

% Add colons to beginning and end of path in case FSbin is 
% at beginning or end and not bracketed by them

for path2check = {FSbin}
    sf = strfind([ ':' pth ':' ],[ ':' path2check{1} ':' ]);
    if (isempty(sf))
        pth = [ pth ':' path2check{1} ];
    end
end


%% Set Up Environment
ENV = {...
    'PATH',pth;...
    'SUBJECTS_DIR',freesurfersubjectsdir;...
    'FREESURFER_HOME',freesurferdir;...
    'FSF_OUTPUT_FORMAT',freesurferoutputformat;...
    };


FSsetup='';

for setupscript = {deblank(freesurfersetup) deblank(freesurferenvironment)}
    if ~isempty(setupscript{1})
        if setupscript{1}(end)~=';', setupscript{1}=[setupscript{1} ';']; end
        FSsetup=[ FSsetup 'source ' setupscript{1} ];
    end
end


switch freesurfershell

    case 'none'
        for e = 1:size(ENV,1)
            setenv(ENV{e,1},ENV{e,2});
        end
        cmd = [FSsetup FScmd];

    case 'bash'
        cmd=[freesurfershell ' -c '''];
        for e = 1:size(ENV,1)
            cmd = [cmd sprintf('export %s=%s;',ENV{e,1},ENV{e,2})];
        end
        cmd = [cmd FSsetup FScmd ''''];

    case {'csh', 'tcsh'}
        cmd=[freesurfershell ' -c '''];
        for e = 1:size(ENV,1)
            cmd = [cmd sprintf('setenv %s %s;',ENV{e,1},ENV{e,2})];
        end
        cmd = [cmd FSsetup FScmd ''''];
end

fprintf('running command: %s', cmd);

[s,w] = system(cmd);


%% Return Errors
if (s)
    [ ~,wenv ] = system('/usr/bin/env');
    fprintf('*** LINUX ERROR FROM SHELL %s\n\n *** WHILE RUNNING COMMAND\n%s',w,cmd);
    fprintf('\n\n*** WITH ENVIRONMENT VARIABLES\n%s',wenv);
    error(' ');
end



end