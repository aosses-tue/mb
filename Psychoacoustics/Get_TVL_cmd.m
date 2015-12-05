function Get_TVL_cmd(dir,exp2filter)
% function Get_TVL_cmd(dir,exp2filter)
%
% 1. Description:
%       It extracts the instantaneous, short-term and long-term loudness
%       using the Time-Varying Loudness model of Moore and Glasberg.
% 
% 2. Stand-alone example:
%       dir = [Get_TUe_data_paths('piano') '04-PAPA' delim '02-Exported-as-segments' delim 'Csh5' delim];
%       exp2filter = '*.wav';
%       Get_TVL_cmd(dir,exp2filter);
% 
%       dir = [Get_TUe_data_paths('piano') '04-PAPA' delim '02-Exported-as-segments' delim 'A4-mod' delim];
%       exp2filter = '*.wav';
%       Get_TVL_cmd(dir,exp2filter);
% 
%       dir = [Get_TUe_data_paths('piano') '04-PAPA' delim '02-Exported-as-segments' delim 'A4-mod' delim];
%       exp2filter = 'NS19*.wav';
%       Get_TVL_cmd(dir,exp2filter);
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 27/11/2015
% Last update on: 27/11/2015 
% Last use on   : 27/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

if nargin < 2
    exp2filter = '*.wav';
end

if nargin < 1
    % dir = [Get_TUe_data_paths('piano') '04-PAPA' delim '02-Exported-as-segments' delim 'C2' delim];
    dir = [Get_TUe_data_paths('piano') '04-PAPA' delim '02-Exported-as-segments' delim 'C2-mod' delim];
end

dirtvl = Get_TUe_data_paths('tvl');
filenames = Get_filenames(dir,exp2filter);

for i = 1:length(filenames)
    % command4system='tvl -i 1k200r.wav -c 100 -p > 1k200r.out';
    fn      = [dir filenames{i}];
    fnnoext = Delete_extension(fn,'wav');
    
    command4system=sprintf('%sresample_2_32k -i %s -o x32.wav',dirtvl,fn);
    [s r] = system(command4system);
    
    % command4system=sprintf('%stvl -i x32.wav -c 100 -s -10 -p > %s.out',dirtvl,fnnoext); % considers a gain of -10 dB
    command4system=sprintf('%stvl -i x32.wav -c 100 -p > %s.out',dirtvl,fnnoext);
    [s r] = system(command4system);
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
