function h = tf_electrodogam
% function h = tf_electrodogam
% 
% 1. Description:
%       Fig1: 8-Ch Electrodogram of the syllable 'bus' using score
%       Fig2: Spectrogram of 'bus'
%       Fig3: Spectrogram of CI-simulated 'bus'. The simulation was done using
%             a noise vocoder of 8Ch
% 
%       Tested: yes
% 
% Created by Tom Francart
% Edited by Alejandro Osses
% Last use: 10/03/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

fs=16000;
NumCh   = 8; % Num channels vocoder

dirM    = Get_TUe_paths( 'MATLAB' );
dir     = Get_TUe_paths( 'db_speechmaterials' );

addpath( [dirM 'tb_Misc' delim 'score'       delim] ) % where score_process_elec.m is
addpath( [dirM 'tb_Misc' delim 'audis-cisim' delim] ) % where vocoder.m is

% dir = '/mnt/x/ExpORL/speechmaterials/'; % Tom's path in Linux computer

[bus, rate] = wavread([dir 'dutch' delim 'nva' delim '1' delim 'bus.wav']);
bus=resample(bus, fs, rate);
sound(bus,fs);

p.fs=fs;
p.rate = 900;
p.electrodes = [1:8];

[seq,p_nmt]=score_process_elec(bus, p)

bus_vocoder=vocoder(bus,fs,NumCh);
sound(bus_vocoder,fs);

h=figure;
showsequence(seq);
xlim([0 0.9]);
ylim([1 9]);
title('Electrodogram: Bus')
% savefig(f, 'got-ci-electrodogram-bus');

h(end+1)=figure;
spectrogram_l(bus);
title('Bus')
% savefig(f, 'got-ci-spectrogram-bus');


h(end+1)=figure;
spectrogram_l(bus_vocoder);
title('Bus - vocoder')
% savefig(f, 'got-ci-spectrogram-bus-vocoder');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename])

function spectrogram_l(d)
fs=16000;

cmap='jet'; 
Slim=[-70,-1];
w=@hamming; T=[18,1]; nfft=1024; 


s = filter([1 -0.95],1,d);   % apply a typical preemphasis filter

[S,F,T,P] = spectrogram(s, T(1)/1000*fs, (T(1)-T(2))/1000*fs, 1024, fs, 'yaxis');
% set(gca,'yscale', 'log')
% set(gca,'YTick', [250 500 1000 2000 4000]);
% ylim([200 4500]);

S=abs(S+eps);
S = S/max(S(:));
S = 20*log10(S);
S(S<Slim(1)) = Slim(1);

% surf(T,F, S, 'EdgeColor','none');
imagesc(T,F,S(end:-1:1,:), Slim);
ax=gca;
yticks = [0:1000:8000];
set(ax,'YTick', yticks);
set(ax,'YTickLabel', yticks(end:-1:1));

% axis xy;
% axis tight;
% view(0,90);
% colormap(cmap);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
