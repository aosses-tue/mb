function outs = demo_dau1996b_with_int_noise(options)
% function outs = demo_dau1996b_with_int_noise(options)
%
% 1. Description:
%       Recreates simulations as presented in Dau1996b. If stimuli are not
%       found, then they are along this script generated.
% 
% 2. Stand-alone example:
%       options.bSave = 1;
%       options.stim_durations = [10 20 40]; % ms
%       demo_dau1996b(options);
% 
%       options.bSave = 1;
%       options.stim_durations = [10 20 40 80 160 320]; % ms
%       demo_dau1996b(options);
% 
% 3. Additional info:
%       Tested cross-platform: No
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 07/10/2014
% Last update on: 22/10/2014 % Update this date manually
% Last use on   : 13/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    close all
    options = [];
end

pathaudio = 'D:\Output\r20141031_update_dau_et_al20150318\';
output_dir = ['D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2015-03-20-update\Figures' delim];

% options = Ensure_field(options, 'nExperiment',2); % Exp. 3 - signal integration
% options = Ensure_field(options, 'bSave', 0);
% options = Ensure_field(options, 'bPlot', 1); % just main plot
options = Ensure_field(options, 'method','dau1996'); % dau1996  uses dau1996preproc
                                                     % dau1996a uses dau1996apreproc
% nExperiment = options.nExperiment;
% bPlot = 0; % This is for getting a lot of plots
% 
% options = Ensure_field(options, 'dB_SPL'      , 85);
% options = Ensure_field(options, 'dB_SPL_noise', 77);
% options = Ensure_field(options, 'output_dir', Get_TUe_paths('outputs'));
% paths.outputs   = options.output_dir;
% 
% h = []; % we initialise handle for Figures

%% II.A Deterministic maskers: simultaneous masking

% % Common parameters:
% dB_SPL_noise    = options.dB_SPL_noise; % reference: Left audio file
% dB_SPL          = options.dB_SPL;

opts.method = options.method;

test_onset_ref  = [-100:40:-20,-15:5:15]*1e-3;

dir = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Masking\dau1996b\';
dirResults = [dir 'Results' delim];
dirStimuli = [dir 'Stimuli' delim];

stim1   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-30-ms.wav']; % -20
stim2   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-35-ms.wav']; % -15
stim3   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-40-ms.wav']; % -10
stim4   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-45-ms.wav']; % -5
stim5   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-50-ms.wav']; % 0
stim6   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-55-ms.wav']; % 5
stim7   = [dirStimuli 'dau1996b_expIB0_stim-10ms-76-onset-60-ms.wav']; % 10

stimnoise = [dir 'Stimuli' delim 'dau1996b_expI_noisemasker.wav'];

[xnoise fs] = Wavread(stimnoise);
Sil         = Gen_silence(50e-3,fs);
xnoise      = [Sil; xnoise(1:end-length(Sil))];
L           = length(xnoise);
t           = ( 0:L-1 )/fs;

x1          = Wavread(stim1); x1 = x1(1:L); % it has one more sample
x2          = Wavread(stim2); x2 = x2(1:L); 
x3          = Wavread(stim3); x3 = x3(1:L); 
x4          = Wavread(stim4); x4 = x4(1:L); 
x5          = Wavread(stim5); x5 = x5(1:L); 
x6          = Wavread(stim6); x6 = x6(1:L); 
x7          = Wavread(stim7); x7 = x7(1:L); 

xref1       = From_dB(9)*x1; % 76 + 9 = 85 dB
xref2       = From_dB(9)*x2; % 76 + 9 = 85 dB
xref3       = From_dB(9)*x3; % 76 + 9 = 85 dB
xref4       = From_dB(9)*x4; % 76 + 9 = 85 dB
xref5       = From_dB(9)*x5; % 76 + 9 = 85 dB
xref6       = From_dB(9)*x6; % 76 + 9 = 85 dB
xref7       = From_dB(9)*x7; % 76 + 9 = 85 dB

lvl = 0;
var = From_dB(lvl);

[s1 xx xx outs]  = Get_internal_representations(x1       ,fs,opts.method,var);
n1  = Get_internal_representations(xnoise   ,fs,opts.method,var);
s1n = Get_internal_representations(x1+xnoise,fs,opts.method,var);

idx = outs.idx(1);
[mean( s1(:,idx)  ) std( s1(:,idx) )]
[mean( n1(:,idx)  ) std( n1(:,idx) )]
[mean( s1n(:,idx) ) std( s1n(:,idx) )]

ytemplate = s1n(:,idx)-n1(:,idx);

figure;
subplot(3,1,1)
plot(t,n1(:,idx)); grid on

subplot(3,1,2)
plot(t,s1n(:,idx)); grid on

subplot(3,1,3)
plot(t,ytemplate); grid on
disp('')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
