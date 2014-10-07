function r20141007_FIA_KUL_additional_analysis
% function r20141007_FIA_KUL_additional_analysis
%
% 1. Description:
%       Additional analysis as sent to TF on 07/10/2014
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 07/10/2014
% Last update on: 07/10/2014 % Update this date manually
% Last use on   : 07/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filenames = {'wdz2','wdz5'};

info.bDoNMT = 1;

p.F0mod.F0max = 300;
p.F0mod.F0min = 75;
p.CFG.Fs = 15659;

try
h = []; 
for i = 1:length(filenames)

    [tref  F0ref ] = Get_F0_praat_from_txt(['D:\Databases\dir03-Speech\dutch\list\fda_eval\praat\' filenames{i} '.txt']);
    [tref2 F0ref2] = Get_F0_praat_from_txt(['D:\Databases\dir03-Speech\dutch\list\fda_eval\praat-CP810\CP810' filenames{i} '.txt']);
    tref2 = tref2 - 0.15; % alignment

    idx = find(F0ref2 > p.F0mod.F0max);
    F0ref2(idx) = NaN;

    idx = find(F0ref2 < p.F0mod.F0min);
    F0ref2(idx) = NaN;

    F0ref2  = interp1(tref2, F0ref2, tref);

    %%
    filename = ['D:\Databases\dir03-Speech\dutch\list\alle-zinnen\' filenames{i}];

    [x fs]    = Wavread([filename '.wav']);

    if info.bDoNMT == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NMT Autoc
        %
        %   0 dBFS = 94 dB SPL
        % -25 dBFS = 69 dB
        % -69 dBFS = 25 dB SPL = 3.5481e-4

        % silence_thr = From_dB(-69);
        silence_thr = From_dB(-50);

        pF0m = [];
        pF0m = Ensure_field(pF0m, 'silence_thr', silence_thr); % Should be 0.01
        pF0m.max_f0 = p.F0mod.F0max;
        pF0m.min_f0 = p.F0mod.F0min;
        pF0m = F0mod_map_201401(pF0m);
        pF0m.processes = {'Autoc_opt2_Alt_proc'}; % Just F0 extraction
        pF0m.do_lpf = 1;

        x16k = resample(x,p.CFG.Fs,fs);

        pF0m.audio_sample_rate = p.CFG.Fs;
        out = Autoc_opt2_Alt_proc(pF0m, x16k);
        t_test_NMT = out{4};
        F0test_NMT = out{2};

        idx = find(F0test_NMT > p.F0mod.F0max);
        F0test_NMT(idx) = NaN;

        idx = find(F0test_NMT < p.F0mod.F0min);
        F0test_NMT(idx) = NaN;

        F0test_NMT  = interp1(t_test_NMT, F0test_NMT, tref);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    F0test = F0test_NMT;

    opts        = Get_F0_measures_eTone(tref,F0test,F0ref );
    opts_CP810  = Get_F0_measures_eTone(tref,F0test,F0ref2);

    figure;
    plot(tref,F0ref); hold on
    plot(tref,F0ref2,'b--','LineWidth',4);
    plot(tref,F0test,'r');
    legend(['ref1: ' filenames{i}],['ref2: ' filenames{i} ', CP810'],['test: ' filenames{i} ', AUTOC'])
    grid on
    xlabel('time [s]')
    ylabel('F0 [Hz]')
    h(end+1) = gcf;
    
end

for i = 1:length(h)
    Saveas(h(i),[Get_TUe_paths('outputs') 'F0-test-' filenames{i}]); 
end

catch
    dirs = dirwalk([Get_TUe_paths('MATLAB') 'tb_NMT_4.31']);
    for i = 1:length(dirs)
        addpath(dirs{i});
    end

    dirs = dirwalk([Get_TUe_paths('MATLAB') 'tb_NMT_AddOns']);
    for i = 1:length(dirs)
        addpath(dirs{i});
    end
    error('');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
