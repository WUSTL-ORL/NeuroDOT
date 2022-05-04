function [t4,scale,header]=loadAviT4(filepath,filename)

% this function loads in the t4 file which transforms an mpr to the atlas
% space.
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
s = filesep;  % Update by ZEM 2018/10/16 to circumvent inconsistent behavior of getenv('OS').
if ~isempty(filepath)
    if ~strcmp(filepath(end),s)  % Update by ZEM 2018/10/16 to make general across operating systems.
        filepath=[filepath,s];
    end
end

%%
fid=fopen([filepath,filename]);
tline = fgetl(fid);
header={};
while ischar(tline)
    if tline(1) ~='t' && tline(1)~='s' && tline(1)~=' ' 
        header=cat(1,header,tline(1:end));
    end
    if tline(1:2)=='t4'
        for i=1:4
            tline = fgetl(fid);
            t4(i,:)=str2num(tline);
        end
    end
    if tline(1:5)=='scale'
        scale=str2num(tline(7:end));
    else scale=[];
    end
    tline = fgetl(fid);
end
fclose(fid);