function outs = demo_dau1997b(options)
% function outs = demo_dau1997b(options)
%
% 1. Description:
%       Recreates simulations as presented in Dau1996b. If stimuli are not
%       found, then they are along this script generated.
% 
% 2. Stand-alone example:
%       options = [];
%       demo_dau1997b(options);
% 
% 3. Additional info:
%       Tested cross-platform: No
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 01/04/2015
% Last update on: 01/04/2015 % Update this date manually
% Last use on   : 01/04/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    close all
    options = [];
end

options = Ensure_field(options, 'nExperiment',2); % Exp. 3 - signal integration
options = Ensure_field(options, 'bSave', 0);
options = Ensure_field(options, 'bPlot', 1); % just main plot
options = Ensure_field(options, 'method','dau1996'); % dau1996  uses dau1996preproc
                                                     % dau1996a uses dau1996apreproc
nExperiment = options.nExperiment;
bPlot = 0; % This is for getting a lot of plots

options = Ensure_field(options, 'dB_SPL'      , 85);
options = Ensure_field(options, 'dB_SPL_noise', 77);
options = Ensure_field(options, 'output_dir', Get_TUe_paths('outputs'));
paths.outputs   = options.output_dir;

h = []; % we initialise handle for Figures

bListen = 1;
fs = 44100;

%% Dau1997b, Fig. 7:
% Suprathreshold signal...

% Common parameters:
BW = 3; %[3 31 314];
% d = [-40];
mdept = 1;%d2m(d,'dau');

fc = 5000;
finf = fc-BW/2;
fsup = fc+BW/2;
SPL = 65;
dur = 4;
fm  = 20;
insig_NBN   = AM_random_noise_BW(fc,BW,SPL,dur,fs,fm,mdept);
insig_test  = AM_sine(fc,dur,fs,fm,mdept,SPL);
 
t = ( 0:length(insig_NBN)-1 )/fs;

if bListen == 1
    sound(insig_NBN,fs);
end

%%

[outsig1 , fc ,fcm, opts] = dau1997preproc_1Ch(insig_NBN           ,fs,5000);
[outsig2 , fc ,fcm, opts] = dau1997preproc_1Ch(insig_test+insig_NBN,fs,5000);
outsig1 = outsig1{1};
outsig2 = outsig2{1};

opts.step1 = 1;
opts.step2 = 1;
figure;
subplot(3,1,1)
Mesh([1:12],t,outsig1,opts);
ylabel('Time [s]')
xlabel('Modulation filter');
zlabel('[MU]')

subplot(3,1,2)
Mesh([1:12],t,outsig2,opts);
ylabel('Time [s]')
xlabel('Modulation filter');
zlabel('[MU]')

disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [y,t] = Create_sin4this_exp(f,start_phase,dur,fs,win,onset,SPL,total_duration)

    [y, t]= Create_sin_phase(f,start_phase,dur,fs,win);
    
    % (f,dur,fs,win,onset,dB_SPL_above_thr);
    y  = setdbspl(y,SPL);
    y  = [Gen_silence(onset,fs); y]; 
    
    try % Append silence only if total_duration has been specified
    	y = [y; Gen_silence(total_duration-max(t)-onset-1/fs,fs)];
    end

    t = (1:length(y))/fs; % redefine t
    
end
