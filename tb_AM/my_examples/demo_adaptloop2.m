function h = demo_adaptloop2(opts)
% function h = demo_adaptloop2(opts)
%
% 1. Description:
%       Show the effect of adaptation: This script demonstrates the effect 
%       of adaptation applied to a test signal with and without noise.
%
%   The test signal is made of a sinosoidal ramp up and down between 0 and 1.
%
%   Figure 1: Clean test signal
%
%      This figure shows the effect of adaptation on the clean test signal 
%      with and without overshoot limiting.
%
%   Figure 2: Noisy test signal
%
%      This figure shows the effect of adaptation on the noisy test signal
%      with and without overshoot limiting. Notice that in the second plot,
%      the initial spike at the beginning of the signal caused from the sharp
%      transition from complete silence to noise is magnitudes larger than
%      the values in the rest of the output.
%
%   See also: adaptloop
%
% Author        : Alejandro Osses V., HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 24/03/2015
% Adapted from  : dau1996preproc.m (code by Torsten Dau, Morten L. Jepsen, Peter L. Sondergaard)
% Last update on: 24/03/2015 % Update this date manually
% Last use on   : 30/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    close all
    opts = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot save/options:
opts        = Ensure_field(opts, 'bSave', 0);
opts        = Ensure_field(opts, 'bPlot', 1);
h           = []; % we initialise handle for Figures
bPlot       = opts.bPlot;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal options:
opts        = Ensure_field(opts, 'dB_SPL'     , 77);
opts        = Ensure_field(opts, 'dB_SPL_test', 85);
opts.f      = 1500;
fc          = opts.f;

paths.outputs   = Get_TUe_paths('outputs');
% infilename      = sprintf('%.0f-sine-65-dBSPL',opts.f); 
infilename      = ['kohlrausch1992_noisemasker-' num2str(opts.dB_SPL)     ];
infilename2     = ['kohlrausch1992_testsignal-'  num2str(opts.dB_SPL_test)];

fname   = [paths.outputs infilename  '.wav'];
% fname2  = [paths.outputs infilename2 '.wav'];

try
    [insig   fs] = Wavread(fname);
    %[insig2  fs] = Wavread(fname2);
catch
    opts.outputs = paths.outputs;
    tmp = demo_adaptloop2_gen_stim(opts);
    fname = tmp.fname;
    [insig fs] = Wavread(fname);
    % AM_random_noise(Finf,Fsp,SPL,dur,fs,fmod,Mdept,dBFS);
end

title1 = 'White noise';
ymin = -0.15;
ymax =  0.15;
yminMU = -100;
ymaxMU = 1500;
    
t = ( 1:length(insig) )/fs;
t = t(:);

% Dau's options:
opts = Ensure_field(opts,'fc_idx',fc);

[xx, fc , outsig1] = kohlrausch1992preproc(insig,fs);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This next lines were implemented in kohlrausch1992preproc
% tau_LPF = 200e-3;
% mlp_a = exp(-(1/tau_LPF)/fs);
% mlp_b = 1 - mlp_a;
% mlp_a = [1, -mlp_a];
% 
% % Apply the low-pass modulation filter.
% out1 = filter(mlp_b,mlp_a,outsig1.out03_adaptloop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx     = max(find(fc<opts.fc_idx));
fcentre = fc(idx);

idx1    = max(find(fc<1000));
idx2    = max(find(fc< 500));

idx_test = [idx idx1 idx2];

% Dau1996compare(insig,insig2,fs);
if bPlot
    
    for i = idx_test 
        K = 4096*2;
        options.fs = fs;

        %% Similar to Dau1996a, Figure 6
        % Apparently this filter bank is the replaced one (see Dau1997b)
        figure;
        subplot(2,1,1)
        options.bNewFigure = 0;
        freqfft(outsig1.out01_filterbank(:,idx_test),K,options);
        hold on

        % legend(title1)
        xlim([0 8000])
        ylim([-30 20])
        xlabel('Frequency [Hz]')
        ylabel('Amplitude'), hold on

        subplot(2,1,2)
        plot(t,outsig1.out01_filterbank(:,idx_test)); grid on
        xlabel('Time [s]')
    
        legend('1500','1000','500')
        h(end+1)=gcf;

        %% Similar to Dau1996a, Figure 7

        figure
        subplot(2,1,1)
        plot(   t, outsig1.out03_adaptloop(:,idx_test) ); hold on; grid on

        subplot(2,1,2)
        plot(   t, outsig1.out04_LPF(:,idx_test) ); hold on; grid on
        
        legend('1500','1000','500')
        h(end+1)=gcf;
        haxis = gca;
    end
    linkaxes(haxis,'x');
end

% linkaxes(ha,'x');

% h(end+1) = gcf;
% 
% h(end+1) = Figure2paperfigure( [h(end-1) h(end)],3,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end