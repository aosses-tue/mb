function outs = demo_dau1996b(options)
% function outs = demo_dau1996b(options)
%
% 1. Description:
%       Recreates simulations as presented in Dau1996b
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       options.bSave = 1;
%       demo_dau1996b(options);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 07/10/2014
% Last update on: 15/10/2014 % Update this date manually
% Last use on   : 15/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    close all
end

options = Ensure_field(options, 'bSave', 0);
options = Ensure_field(options, 'bPlot', 1);
options = Ensure_field(options, 'dB_SPL', 85);

infilename         = 'dau1996b_expII3_noisemasker';
infilename1        = ['dau1996a_expII3_stim01-' num2str(options.dB_SPL)];
infilename2        = ['dau1996a_expII3_stim02-' num2str(options.dB_SPL)];
infilename3        = ['dau1996a_expII3_stim03-' num2str(options.dB_SPL)];

h = []; % we initialise handle for Figures
paths.outputs   = Get_TUe_paths('outputs');

%% II.A Deterministic maskers: simultaneous masking

% TO DO: 1. Temporal position of short signals in frozen-noise maskers
% TO DO: 2. Relative phase at a fixed temporal position of the signal in the masker
% TO DO: 4. Signal frequency

%% 3. Signal integration
% Stimuli:
%   Noise, white noise, 600 ms, BPF between 20-5000 Hz, 77 dB SPL
%   Tone 1, 10-ms, 3 kHz, hanning windowed
%   Tone 2, 20-ms, 3 kHz, hanning windowed
%   Tone 3, 40-ms, 3 kHz, hanning windowed
%   The onset of the tones was always 100 ms after masker onset

try
    [innoise fs] = Wavread([paths.outputs infilename  '.wav']);
    [instim1 fs] = Wavread([paths.outputs infilename1 '.wav']);
    [instim2 fs] = Wavread([paths.outputs infilename2 '.wav']);
    [instim3 fs] = Wavread([paths.outputs infilename3 '.wav']);
    options.bGenerate   = 0;
catch
    options.bGenerate   = 1;
    fs          = 48000;
end

dB_SPL_noise     = 77; % reference: Left audio file
dB_SPL = options.dB_SPL;

N   = 4096*2; % N-FFT points
K   = N/2;

options.fs = fs;
options.typeplot = 2; % Linear scaled

filename    = [paths.outputs infilename];
filename1   = [paths.outputs infilename1];
filename2   = [paths.outputs infilename2];
filename3   = [paths.outputs infilename3];

if options.bGenerate
    t_silence_bef   = 0e-3;
    t_duration      = 600e-3;
    t_silence_aft   = 0e-3;
    t_total_duration = t_silence_bef + t_duration + t_silence_aft;

    Nsil_bef    = round(options.fs*t_silence_bef);
    Nnoise      = round(options.fs*t_duration);
    Nsil_aft    = round(options.fs*t_silence_aft);

    %% Generating the noise

    title1 = 'White noise';
    ymin = -0.15;
    ymax =  0.15;
    yminMU = -100;
    ymaxMU = 1500;
    
    % Gen1: white noise, band-pass filtered
    y = wgn(Nnoise,1,1);
    
    y   =  y(:); % ensures it is a column vector
    
    Wn = [20 5000]/(options.fs/2); % Normalised cutoff frequency        
    [b,a] = butter(4,Wn); % 8th-order
    y = filtfilt(b,a,y); % Linear-phase implementation
    
    y   = setdbspl(y,dB_SPL_noise);
    innoise = [Gen_silence(t_silence_bef,fs); y; Gen_silence(t_silence_aft,fs)]; % silence at the beginning and at the end
    
    Wavwrite(innoise,fs,filename);
    
    %% Generating the test tones
    
    % Common stim params
    
    onset   = 100e-3;
    f       = 3000;
    win     = 1; % 1 = Hanning window
    
    % Stim 1
    dur     = 10e-3;
    [instim1, t1] = Create_sin4this_exp(f,dur,fs,win,onset,dB_SPL,t_total_duration);
    Wavwrite(instim1,fs,filename1);
    
    % Stim 2
    dur     = 20e-3;
    [instim2, t2] = Create_sin4this_exp(f,dur,fs,win,onset,dB_SPL,t_total_duration);
    Wavwrite(instim2,fs,filename2);
    
    % Stim 3
    dur     = 40e-3;
    [instim3, t3] = Create_sin4this_exp(f,dur,fs,win,onset,dB_SPL,t_total_duration);
    Wavwrite(instim3,fs,filename3);
    
end

out_stim1 = Dau1996compare(innoise,instim1,fs,options.bPlot);

bPlot = 0; % To get values but not to plot
out_stim2 = Dau1996compare(innoise,instim2,fs,bPlot);
out_stim3 = Dau1996compare(innoise,instim3,fs,bPlot);

%%

ha = [];
titlefigure = [];

idx = out_stim1.idx;

if options.bPlot
    figure;
    subplot(3,1,1)
    plot(out_stim1.t,out_stim1.template(:,idx)), grid on
    ha(end+1) = gca;
    title('(a) signal duration: 10 ms')

    subplot(3,1,2)
    plot(out_stim2.t,out_stim2.template(:,idx)), grid on
    ha(end+1) = gca;
    ylabel('normalised amplitude')
    title('(b) signal duration: 20 ms')

    subplot(3,1,3)
    plot(out_stim3.t,out_stim3.template(:,idx)), grid on
    ha(end+1) = gca;
    xlabel('Time relative to masker onset [s]')
    title('(c) signal duration: 40 ms')

    h(end+1) = gcf;
    titlefigure{end+1} = 'dau1996b-fig4-simulated';

    AmpMax = 0.03;
    linkaxes(ha,'xy');
    axis([0 max(out_stim1.t) AmpMax*[-0.25 1]])

    h(end) = Figure2paperfigure(h(end),3,1); % replaces handle
end

outs.out_stim1 = out_stim1;
outs.out_stim2 = out_stim2;
outs.out_stim3 = out_stim3;

if options.bSave
    for i=1:length(h)
        try
            Saveas(h(i), [paths.outputs titlefigure{i}]);
        catch
            Saveas(h(i), [filename '-handle-' num2str(i)]);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [y,t] = Create_sin4this_exp(f,dur,fs,win,onset,SPL,total_duration)

    [y, t]= Create_sin(f,dur,fs,win);
    
    % (f,dur,fs,win,onset,dB_SPL_above_thr);
    y  = setdbspl(y,SPL);
    y  = [Gen_silence(onset,fs); y]; 
    
    try % Append silence only if total_duration has been specified
    	y = [y; Gen_silence(total_duration-max(t)-onset-1/fs,fs)];
    end

end