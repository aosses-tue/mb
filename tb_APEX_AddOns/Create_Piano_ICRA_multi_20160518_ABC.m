function Create_Piano_ICRA_multi_20160518_ABC
% function Create_Piano_ICRA_multi_20160518_ABC
%
% 1. Description:
%       Generate APE experiments to compare triads.
% 
% 2. Stand-alone example:
%   Create_Piano_ICRA_multi_20160518_ABC
% 
% 3. Additional info:
%       Tested cross-platform: No
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 28/04/2016
% Last update on: 28/04/2016 
% Last use on   : 28/04/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_main   = 'C:\Users\aosses\Dropbox\LISTENING-TESTS\SOUNDS-AC\06-Exp-TUe-1-similarity\'; % [Get_TUe_data_paths('ex_Experiments') '2015-APEX-my-experiments' delim 'Piano' delim 'Pilot-ICRA-v2' delim];
dir_stimuli= [dir_main '01-after-loudness-balancing'    delim];
% dir_dst    = [dir_main '02-to-be-used-in-AB-comparison' delim];
dir_dst    = [Get_TUe_paths('outputs')  '02-to-be-used-in-AB-comparison' delim];
Mkdir(dir_dst);

%           0        1       2          3            4         5       6   
pianos = {'GH05','GRAF28','JBS36','JBS51-4486','JBS51-4544','JBS73','NS19'};

rampin_len = 50;   % ms
rampout_len = 150; % ms
dur_piano_samples = 1.3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

files = Get_filenames(dir_stimuli,'*.wav');

for i = 1:length(files)

    fprefix = strsplit(files{i},'.');
    fprefix = fprefix{1};
    
    fname = [dir_stimuli fprefix '.wav'];
    [insig fs] = Wavread(fname);
    outsig = insig(1:round(dur_piano_samples*fs));
    outsig = Do_cos_ramp(outsig,fs,rampin_len,rampout_len);
    
    fnameout = sprintf('%s%s-dur-%.0f-ms.wav',dir_dst,fprefix,dur_piano_samples*1000);
    
    Wavwrite(outsig,fs,fnameout);
    
end

disp('')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
