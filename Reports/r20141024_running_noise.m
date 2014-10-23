function [m s statsT] = r20141024_running_noise(nStimuli, method)
% function [m s statsT] = r20141024_running_noise(nStimuli, method)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       nStimuli = 1; % determine thresholds for 10 ms test stimuli
%       nStimuli = 2; % determine thresholds for  5 ms test stimuli
%       r20141024_running_noise(nStimuli);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 23/10/2014
% Last update on: 23/10/2014 % Update this date manually
% Last use on   : 23/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

if nargin == 0
    close all
    method = 'dau1996'; % without overshoot limiting
end

%% Reading and adjusting onset of test stimuli

lvl_ref = 85; 
idx = 13; 
    
pathaudio = [Get_TUe_paths('outputs') 'white-run' delim];
Mkdir(pathaudio);

if nargin == 0
    nStimuli = 4;
end

insigN = [];
RM = [];
statsMean   = [];

statsStd_intra = []; % intra audio files, in a segment

for i = 1:nStimuli
    filename = [pathaudio method '-white-BP-run-' num2str(i)];
    
    try 
        [x fs] = Wavread([filename '.wav']);
    catch
        options.bSaveNoise      = 1;
        options.fs              = 48000;
        options.dB_SPL_noise    = 77;
        nTag                    = 3;
        [x fs] = Create_noise_dau1996(nTag,filename,options);
    end
    
    hopsize = 10e-3*fs;
    bufsize = 20e-3*fs;
        
    insigN = [insigN x]; % just to store noises
    
    [RMtmp fc t] = Get_internal_representations(x,fs,method);
    RM = [RM RMtmp(:,idx)];
    
    [statsN statsT]   = Get_stats_from_audio_excerpt(RM(:,end) , fs, bufsize, hopsize);
    
    statsT = statsT(:);
    statsMean       = [statsMean statsN.mean'];
    statsStd_intra  = [statsStd_intra statsN.std'];

end

[m s] = barweb_prepare_data(statsMean);

if nargout == 0
    figure;
    barweb(m,s);

    figure;
    plot(statsT,s,'o-')
    xlabel('time [s]')
end

statsStd = statsMean';

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
