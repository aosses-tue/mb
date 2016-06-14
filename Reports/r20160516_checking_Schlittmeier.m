function r20160516_checking_Schlittmeier
% function r20160516_checking_Schlittmeier
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 16/05/2016
% Last update on: 16/05/2016 
% Last use on   : 16/05/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bPart4_real = 1;
dir_stim = 'D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2016-09-05-ICA\Stimuli\';
N = 88200;

if bPart4_real
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    files = {'Spr_1_1-mono.wav', 1.11, 24; ... % Ni = 1058400 (24 s), lvl = 56.1 dB, 67.2
             'Spr_1_2-mono.wav', 1.21,  9; ... % Ni =  396900 ( 9 s), lvl = 60.0 dB, 69.4
             'Spr_2-mono.wav'  , 0.38,  0; ... % Ni =       1       , lvl = 63.6 dB, 67.8
             'Tier1-mono.wav'  , 1.77,0.5; ... % Ni =   22050 (0.5 s), lvl = 64.5 dB, 73.4 dB (peak)
             'RR2-mono.wav'    , 0.02,  0; ... % Ni = 1, lvl = 60.1 dB
             'M5-mono.wav'     , 0.56,  1; ... % Ni =   44100 (1 s), lvl = 58.2 dB
             'M6_2-mono'       , 0.21, 26};    % Ni =  1146600 (26 s) lvl = 62.1 dB
    
	% idx = size(FS,1);
    for i = 1:length(files)
        [insig fs] = Wavread([dir_stim files{i,1}]);
        
        Ni = round(files{i,3}*fs)+1;
        t = (1:length(insig))/fs; t = t(:);
        % t_b = buffer(t,N,round(5*N/10),'nodelay');
        % insig_b = buffer(insig,N,round(5*N/10),'nodelay');
        
        insig = From_dB(-6)*insig(Ni:Ni+N-1,1);
        SPL(i) = rmsdb(insig)+100;
        SPLpeak(i) = max(To_dB(abs(insig)/2e-5));
        
        % [FS(idx+i,1) fi outs] = FluctuationStrength_Garcia(insig, fs, N);
        % FS_theo(idx+i,1) = files{i,2};
        % 
        % figure;
        % plot(1:47,outs.mdept,1:47,outs.kp); grid on
        % legend('m','k')
        % title(num2str(idx+i))
    end

    disp('')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
