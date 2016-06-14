function r20160507_analysis_piano
% function r20160507_analysis_piano
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 07/05/2016
% Last update on: 07/05/2016 
% Last use on   : 07/05/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

% dir = 'D:\Databases\dir01-Instruments\Piano\01-Chabassier\SONS\F3\';
% file = [dir 'pressionsimu.wav'];

% file = 'D:\Databases\dir01-Instruments\Piano\04-PAPA\03-Exported-as-segments\C2\GH05-C2_1.wav';
file = 'D:\Databases\dir01-Instruments\Piano\04-PAPA\03-Exported-as-segments\F3\JBS50-F3_2.wav';
[x fs] = Wavread(file);

N = (2^12);
Ni = Get_max_waveform(x,fs,10e-3);
Nf  = Ni+N-1; 

% file = [Get_TUe_data_paths('piano') '01-Chabassier' delim 'SONS' delim 'Cd5' delim 'pressionexpe.wav'];
% [insig fs] = Wavread(file);

insig = transpose(x);

[AdB F] = searchp_impulseV2(insig,fs,Ni,Nf);
figure;
plot(F,AdB,'r*','LineWidth',2); hold on, grid on
ha = gca;
%%%

insig = x(Ni:Nf);
L = 200;
p = 2*L;

[outsig Fi Ai sigma_i Phi_i L Ai_dB] = Get_ESPRIT_piano(insig,fs,-50,N,L);

stem(Fi,Ai_dB,'BaseValue',-100); grid on
xlabel('Frequency [Hz]')

ha(end+1) = gca;

linkaxes(ha,'xy');
xlim(minmax(F))
ylim([-48 3])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
