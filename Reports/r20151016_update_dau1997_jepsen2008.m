function r20151016_update_dau1997_jepsen2008(bParts)
% function r20151016_update_dau1997_jepsen2008(bParts)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 13/10/2015
% Last update on: 13/10/2015 
% Last use on   : 13/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;

dir_where  = [Get_TUe_paths('outputs') 'audio-20151006' delim];
dir_where7 = [Get_TUe_paths('outputs') 'audio-20151006-forward' delim];
dir_out    = [Get_TUe_paths('lx_Text') 'lx2015-10-16-decision-CASP' delim 'Figures-new' delim];

%%% Creating stimuli for Fig 2
fname1 = [dir_where 'NBN-fc-5000-Hz-BW-3-Hz-65-dB'];
fname2 = [dir_where 'NBN-fc-5000-Hz-BW-31-Hz-65-dB'];
bCreate = 0;
if bCreate
    masker = exp_dau1997a(fs,3);
    Mkdir(dir_where);
    Wavwrite(masker,fs,fname1);
    
    masker = exp_dau1997a(fs,4);
    Wavwrite(masker,fs,fname2);
end

%%% Creating stimuli for Fig 3
if bCreate
    [masker BWs] = exp_dau1997b(fs,2);
    for i = 1:size(masker,2)
        fname = sprintf('NBN-fc-1000-Hz-BW-%.0f-Hz-65-dB',BWs(i));
        Wavwrite(masker(:,i),fs,fname);
    end
    % Wavwrite.m: file NBN-fc-5000-Hz-BW-10-Hz-65-dB created
    % Wavwrite.m: file NBN-fc-5000-Hz-BW-100-Hz-65-dB created
    % Wavwrite.m: file NBN-fc-5000-Hz-BW-800-Hz-65-dB created
    % Wavwrite.m: file NBN-fc-5000-Hz-BW-1000-Hz-65-dB created
    % Wavwrite.m: file NBN-fc-5000-Hz-BW-5000-Hz-65-dB created
end

if nargin == 0
    %         1 2 3 4 5 6 7 8
    bParts = [0 0 0 0 0 0 1 0];
    % CC:     0 1 0 - - 1 0 -
    %                   0       % if 2000-Hz  
end

fs = 44100;

bPart1 = bParts(1); % Calibration of the model: sigma = 1.85 for bDecisionMethod = 2
bPart2 = bParts(2); % Simultaneous masking, deterministic
bPart3 = bParts(3); % Forward masking: on and off-frequency
bPart4 = bParts(4); % plotting results of bPart3
bPart5 = bParts(5); % calibration using 800-Hz tone
bPart6 = bParts(6); % signal integration, 2-kHz tone
bPart7 = bParts(7); % Transition simultaneous-forward masking. Similar to bPart3
bPart8 = bParts(8); % plotting results of bPart7 and all the other new results

opts.nAnalyser      = 103; % 101 = modfilterbank, 103 - jepsen2008

opts.bDecisionMethod = 2; % 2 - cc; 4 - dprime

switch opts.nAnalyser
    case 100
        opts.sigma   = 0.4; % 0.4 for multi-channel, 1.5 for single channel
                             
    case 101
        opts.sigma   = 0.68;    
                                
        opts.modfiltertype = 'dau1997wLP';

    case {103, 104}
        opts.sigma   = 0.38; % 0.59; % 0.615 = 17.86 dB; % 0.45 if gain_after_drnl = 13 dB
    
end
    
opts.var            = opts.sigma.*opts.sigma;
opts.audio.fs       = fs;
if mod(opts.nAnalyser,1)==0
    opts.MethodIntRep   = 1; % 1 - my method; 2 - using casptemplate.m
else
    opts.MethodIntRep   = 2;
end

opts.Reversals4avg  =  6; % 10

opts.do_template    =  1;
opts.do_simulation  =  1;
opts.Nreversals     =  8;

bDebug = 0;
opts.bDebug = bDebug;

close all
count_saved_figures = 1;
if bPart1
    
    tic
    for j = 1 % 0:1
        % To do: save template
        bTones = j;
        bBBN = ~bTones;
        if bTones
            opts.DurRamps   = 125; % additional cosine ramps
            opts.bUseRamp   = 1; % additional cosine ramps
            opts.bUseRampS  = 1; % additional cosine ramps
            erbc2analyse    = freqtoaud([500 2000],'erb'); % 14 for 1000 Hz (approx.)
            opts.filename1 = [Get_TUe_paths('outputs') 'AMTControl-examples' delim 'tone-f-1000-Hz-at-60-dB-dur-800-ms.wav'];
            opts.filename2 = [Get_TUe_paths('outputs') 'AMTControl-examples' delim 'tone-f-1000-Hz-at-42-dB-dur-800-ms.wav'];
            opts.Ntimes     = 1;
            opts.Nsim       = 1;
        end
        if bBBN
            opts.DurRamps   = 0; % additional cosine ramps
            opts.bUseRamp   = 0; % additional cosine ramps
            opts.bUseRampS  = 0;
            erbc2analyse    = freqtoaud([100 8000],'erb'); % 14 for 1000 Hz (approx.)
            % opts.filename1 = [Get_TUe_paths('outputs') 'audio-20151006' delim 'jepsen2008-BW-at-60-dB-dur-500-ms.wav'];
            % opts.filename2 = [Get_TUe_paths('outputs') 'audio-20151006' delim 'jepsen2008-BW-at-42-dB-dur-500-ms.wav']; 
            opts.filename1 = [Get_TUe_paths('outputs') 'audio-20151006' delim 'jepsen2008-BW-at-60-dB-dur-500-ms-2.wav'];
            opts.filename2 = [Get_TUe_paths('outputs') 'audio-20151006' delim 'jepsen2008-BW-at-42-dB-dur-500-ms-2.wav']; 
            opts.Ntimes     = 1;
            opts.Nsim       = 1;
        end

        opts = il_get_freqs(erbc2analyse,opts);

        opts.Gain4supra = 16; % previous results with Gain4supra = 5 dB
        opts.audio.fs   = fs;

        opts.StepdB     = 8; 
        opts.StepdBmin  = 0.5;

        refSPL  = 60; % [20 30 40 50 60 70];

        for i=1:length(refSPL)
            testSPL = refSPL(i)-18; % 42 dB for 60 dB    
            opts.Gain2file1 = refSPL(i) - 60;
            opts.Gain2file2 = refSPL(i) - 60;
            tr_tmp = AMTControl_cl(opts);

            th_JND(j+1,i) = sum_dB_arit([refSPL(i) testSPL+tr_tmp.Threshold]) - refSPL(i);

            disp('')
        end
    end
    toc
    
    % Jepsen2008:
    % 0.6639    0.5314    0.5941    0.6639    0.8279    1.0299 % BBN
    % 0.8279    0.6639    0.5314    0.5314    0.5314    0.5314 % sine
    %
    % Dau1997 - elapsed time = 616 s
    % 0.5941    0.5314    0.6639    0.7623    0.8279    0.8279 % BBN
    % 1.2778    1.0299    0.8279    0.7416    0.6639    0.5941 % sine

    % 0.5:  0.5941    0.4248    0.5314    0.5314    0.6639    0.8279
    %       0.6639    0.5314    0.4752    0.4248    0.4752    0.4752
    
    % 0.55: 0.5941    0.4752    0.5314    0.5941    0.7416    0.9237
    %       0.7416    0.5941    0.5314    0.4752    0.5314    0.4752
    
    % 0.59: 0.6639    0.5314    0.5941    0.6639    0.8279    1.0299
    %       0.8279    0.6639    0.5314    0.5314    0.5314    0.5314
    
    % 0.65: 0.7416    0.5941    0.6639    0.7416    0.9237    1.1476
    %       0.9237    0.6639    0.5941    0.5941    0.5941    0.5941
    
    % 0.7:  0.8279    0.5941    0.7416    0.7416    0.9237    1.1476
    %       0.9237    0.7416    0.6639    0.5941    0.6639    0.6639
    
    % 0.8:   0.9237    0.7416    0.8279    0.8279    1.0299    1.4216
    %        1.1476    0.8279    0.7416    0.7416    0.7416    0.7416
    
    if nargout > 0
        th = th_JND; 
        return;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modulation detection:
close all
opts.StepdB    = 4; 
opts.StepdBmin = 1;
if bPart2
    
    erbc2analyse  = freqtoaud([5000],'erb'); 
    opts = il_get_freqs(erbc2analyse,opts);
    fmodtest = [3 5 10 20 50 100 120];
    
    opts.DurRamps   = 200; % additional cosine ramps
    opts.Gain4supra =  -1; % when creating: lvl = 75 dB; lvl supra = 75 dB + Gain4supra
    opts.audio.fs   =  fs;
    
    opts.bUseRamp  = 1; % additional cosine ramps
    opts.bUseRampS = 1; % additional cosine ramps
    opts.increment_method = 'modulation-depth';
    
    for i = 1:length(fmodtest)
        if strcmp(opts.increment_method,'modulation-depth')
            opts.fmod = fmodtest(i); % Hz
            opts.dur_test = 1; % s
        end

        opts.Ntimes = 25; % 5
        opts.Nsim = 4;
        opts.filename1 = [dir_where 'NBN-fc-5000-Hz-BW-31-Hz-65-dB.wav']; % 'NBN-fc-5000-Hz-BW-31-Hz-65-dB.wav'
        tmp = AMTControl_cl(opts);
        
        th_TMTF(i) = median(tmp.Threshold);
        th_TMTF_Prct75(i) = prctile(tmp.Threshold,75);
        th_TMTF_Prct25(i) = prctile(tmp.Threshold,25);
        
        disp('')
    end
    
    if nargout > 0
        th = th_det;
        return;
    end
    
%%% 3-Hz, Dau1997
% th_TMTF        =  -8.7500  -18.7500  -29.5000  -25.2500  -26.7500  -29.7500  -26.0000
% th_TMTF_Prct75 =  -7.2500  -18.2500  -28.2500  -22.5000  -25.2500  -29.0000  -24.0000
% th_TMTF_Prct25 = -10.7500  -19.2500  -30.2500  -26.7500  -29.5000  -30.2500  -27.7500
    
%%% 3-Hz, Jepsen2008, run 1 (old)
% th_TMTF        =  -0.7500  -13.7500  -25.0000  -16.5000  -19.7500  -21.2500  -23.5000
% th_TMTF_Prct75 =   0.7500  -12.7500  -23.5000  -14.2500  -16.0000  -19.0000  -23.2500
% th_TMTF_Prct25 =  -2.2500  -14.2500  -26.7500  -17.2500  -24.5000  -22.5000  -23.7500

%%% 3-Hz, Jepsen2008, run 2 (old)
% th_TMTF        = -4.0000  -13.2500  -17.2500  -15.2500  -17.0000  -26.2500  -21.7500
% th_TMTF_Prct75 = -3.7500  -11.5000  -17.0000  -14.2500  -17.0000  -24.5000  -20.0000
% th_TMTF_Prct25 = -4.7500  -14.7500  -18.0000  -16.2500  -18.0000  -28.5000  -23.2500

%%% 31-Hz, Dau1997
% th_TMTF        = -3.2500  -14.0000  -11.2500  -10.5000  -15.7500  -16.2500  -15.0000
% th_TMTF_Prct75 = -2.7500  -12.5000  -10.2500   -9.5000  -15.0000  -15.5000  -14.5000
% th_TMTF_Prct25 = -4.2500  -19.5000  -14.2500  -12.5000  -18.0000  -17.5000  -15.0000



%%% 3-Hz, Dau1997. supraT level = -6 dB (run on 01/11/2015)
% th_TMTF        = -11.5000  -20.2500  -30.0000  -28.5000  -27.5000  -27.7500  -26.2500
% th_TMTF_Prct75 = -10.7500  -19.0000  -28.0000  -27.0000  -26.0000  -26.5000  -25.0000
% th_TMTF_Prct25 = -13.2500  -21.7500  -31.7500  -29.2500  -30.0000  -28.7500  -29.0000

%%% 31-Hz, Dau1997. supraT level = -1 dB (run on 01/11/2015)
% th_TMTF        = -5.5000   -9.5000  -13.7500  -12.7500  -19.2500  -18.0000  -19.2500
% th_TMTF_Prct75 = -4.5000   -7.5000  -12.2500  -12.0000  -18.0000  -18.0000  -16.7500
% th_TMTF_Prct25 = -6.0000  -10.0000  -14.0000  -14.0000  -20.5000  -19.0000  -22.7500

%%% 3-Hz, Jepsen2008. supraT level = -6 dB (run on 01/11/2015)
% th_TMTF        = -11.0000  -13.0000  -22.7500  -22.7500  -18.0000  -22.0000  -16.7500
% th_TMTF_Prct75 = -8.2500  -11.0000  -21.2500  -22.0000  -17.5000  -20.0000  -14.5000
% th_TMTF_Prct25 = -13.0000  -15.2500  -24.2500  -26.5000  -19.5000  -24.7500  -18.2500

%%% 31-Hz, Jepsen2008. supraT level = -1 dB, first fmod might be biased by 0 response (run on 01/11/2015)
% th_TMTF        = -1.7500   -6.5000  -11.2500  -11.2500  -17.5000  -17.2500  -17.0000
% th_TMTF_Prct75 = -1.5000   -5.7500  -10.2500   -9.2500  -16.2500  -16.2500  -16.0000
% th_TMTF_Prct25 = -2.5000   -7.5000  -11.7500  -14.2500  -18.2500  -18.5000  -17.0000

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
if bPart3
    
    erbc2analyse    = freqtoaud([3500 5000],'erb'); % 14 for 1000 Hz (approx.)  % 
    opts            = il_get_freqs(erbc2analyse,opts);
        
    fnamesM = {[dir_where7 'jepsen2008-fig7-4000-Hz-tone-masker-60-dB-dur-10-s.wav'], ... % on-frequency masker
               [dir_where7 'jepsen2008-fig7-2400-Hz-tone-masker-60-dB-dur-10-s.wav'] };   % off-frequency masker
    fnames  = {[dir_where7 'jepsen2008-fig7-4000-Hz-tone-60-dB-dur-250-ms-onset-200-ms.wav'], ...
               [dir_where7 'jepsen2008-fig7-4000-Hz-tone-60-dB-dur-250-ms-onset-230-ms.wav']}; 
   
	testlevels      = [40 60 70 80; ... % levels for on-freq
                       60 70 80 85];    % levels for off-freq
                   
	if opts.nAnalyser == 101 | opts.nAnalyser == 103
        opts.resample_intrep = 'resample_intrep';
    end
    opts.Gain4supra =  10; % 10 dB above the masker level 
    opts.audio.fs   =  fs;
    
    opts.StepdB = 8; 
    opts.StepdBmin = 1;
    
    opts.DurRamps  = 2; % [ms] additional cosine ramps
    opts.bUseRamp  = 1; % additional cosine ramps
    opts.bUseRampS = 0; % additional cosine ramps
    opts.bAddSilence2noise = 1;
    opts.Silence2noise = 500e-3; % 200 + silence = 700. If sil = 400e-3, then signal dur = 300e-3 -> simultaneous masking
    opts.increment_method = 'level';
    
    opts.Ntimes = 25;
    opts.Nsim = 4;
    opts.bDebug = 0;
    %%% k = 1: offset-onset of  0 ms
    %%% k = 2: offset-onset of 30 ms
    for k = 1:2
        opts.filename2 = fnames{k};

        for i = 1:length(testlevels) 
            for j = 1:2
                
                opts.filename1 = fnamesM{j}; % 1 is on-freq, 2 is off-freq
                opts.Gain2file1 = testlevels(j,i)-60;
                opts.Gain2file2 = opts.Gain2file1;
                tmp = AMTControl_cl(opts);

                if k == 1
                    TTh_00ms(j,i) = median(tmp.Threshold)+testlevels(j,i);
                    Prct75_00ms(j,i) = prctile(tmp.Threshold,75)+testlevels(j,i);
                    Prct25_00ms(j,i) = prctile(tmp.Threshold,25)+testlevels(j,i);
                elseif k == 2
                    TTh_30ms(j,i) = median(tmp.Threshold)+testlevels(j,i);
                    Prct75_30ms(j,i) = prctile(tmp.Threshold,75)+testlevels(j,i);
                    Prct25_30ms(j,i) = prctile(tmp.Threshold,25)+testlevels(j,i);
                end

                disp('')
                
            end
        end
    end
    
    disp('')
    
%%%
% Dau1997 (first row=on-freq; second row=off-freq):
% TTh_00ms =
%    26.2500   39.0000   50.5000   63.5000
%    14.0000   26.0000   37.0000   42.0000
% 
% Prct75_00ms =
%    26.5000   39.0000   50.5000   64.0000
%    14.0000   26.0000   37.0000   42.0000
% 
% Prct25_00ms =
%    26.0000   39.0000   50.2500   63.0000
%    14.0000   25.5000   37.0000   42.0000
% 
% 
% TTh_30ms =
% 
%     20    34    43    52
%      8    12    17    21
% 
% 
% Prct75_30ms =
% 
%     20    34    43    52
%      8    12    17    21
% 
% 
% Prct25_30ms =
% 
%     20    34    43    52
%      8    12    17    21  

%%% Jepsen2008:
% TTh_00ms =
%    27.0000   30.0000   33.0000   38.0000
%    21.0000   27.0000   34.2500   43.0000
% 
% Prct75_00ms =
%    27.0000   30.0000   33.0000   38.0000
%    21.0000   27.0000   34.7500   43.0000
% 
% Prct25_00ms =
%     27    30    33    38
%     21    27    34    43
% 
% TTh_30ms =
%     25    30    32    35
%     18    20    25    29
% 
% Prct75_30ms =
%     25    30    32    35
%     18    20    25    29
% 
% Prct25_30ms =
%     25    30    32    35
%     18    20    25    29
end


% Only 4000-Hz channel:
% TTh_00ms =    21.0000   25.0000   28.0000   33.2500
%               20.0000   25.0000   33.0000   38.0000
% Prct75_00ms = 21.0000   25.0000   28.2500   33.7500
%               20.0000   25.0000   33.0000   38.0000
% Prct25_00ms = 21    25    28    33
%               20    25    33    38

% TTh_30ms    = 23    28    30    33
%               16    19    23    27
% Prct75_30ms = 23    28    30    33
%               16    19    23    27
% Prct25_30ms = 23    28    30    33
%               16    19    23    27

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart4
    
    hFig4 = []; % handles for the figures generated in this part
    
    %% 4.a Results of intensity discrimination task
    % raw1 - intensity discrimination with 1-kHz tones
    % raw2 - intensity discrimination with BBN (100 Hz - 8 kHz)
    testlevels = [20 30 40 50 60 70];
    minL = 20; maxL = 70;
    minY = 0.3; maxY = 1.5;
    
    raw1_d1997 = [1.3480    1.0299    0.8509    0.7214    0.6281    0.5465]; % sigma = 0.68
    raw2_d1997 = [0.5941    0.5778    0.7214    0.8233    0.8279    0.8279]; % sigma = 0.68
    raw1_j2008 = [0.6109    0.5168    0.4752    0.4752    0.5168    0.5314]; % sigma = 0.615
    raw2_j2008 = [0.6109    0.6639    0.6281    0.7017    0.8279    1.0873]; % sigma = 0.615
    % raw1_j2008 = [0.7623 0.5204]; % sigma = 0.59 on 30/10/2015, SPL = 20, 60
    % raw2_j2008 = [0.7623 1.0023]; % sigma = 0.59; SPL = 60, 70
    
    figure;
    subplot(1,2,1)
    plot(testlevels, raw1_d1997,'rs-','LineWidth',2); hold on
	plot(testlevels, raw1_j2008,'ko-','LineWidth',2); grid on
	xlabel('Standard level [dB]'); xlim([minL maxL])
    ylabel('\Delta at threshold [dB]'); ylim([minY maxY])
    title('A. test 1-kHz tones');
    set(gca,'XTick',testlevels);
    set(gca,'YTick',[minY:0.1:maxY])
    
    subplot(1,2,2)
    plot( testlevels, raw2_d1997,'rs-','LineWidth',2); hold on
    plot( testlevels, raw2_j2008,'ko-','LineWidth',2); grid on
    xlabel('Standard level [dB]'); xlim([minL maxL])
    ylabel('\Delta at threshold [dB]'); ylim([minY maxY])
    title('B. test BBNs');
    set(gca,'XTick',testlevels);
    set(gca,'YTick',[minY:0.1:maxY])
    
    legend('PEMO', 'CASP','Location','NorthWest');
    hFig4(end+1) = gcf;
        
    %%%
    sigmas      = [0.5 0.6 0.615 0.63 0.7 0.85];
    LineColor  = {'bs-';'m>-.';'ko-';'r--';'g>-';'r<-'};
    LineWidth  = [1.5 1 2 1 1 2];
    
    raw3_j2008 = [0.4752    0.4015    0.3796    0.3796    0.4015    0.4248; ... % sigma = 0.5
                  0.5778    0.4887    0.4369    0.4369    0.4621    0.4752; ... % sigma = 0.6
                  0.6109    0.5168    0.4752    0.4752    0.5168    0.5314; ... % sigma = 0.615
                  0.6109    0.5314    0.4752    0.4887    0.5285    0.5465; ... % sigma = 0.63
                  0.6639    0.5619    0.5314    0.5314    0.5619    0.5941; ... % sigma = 0.7
                  0.8279    0.7017    0.6281    0.6639    0.7017    0.7416];    % sigma = 0.85
 
    raw4_j2008 = [0.5026    0.5314    0.5026    0.5619    0.6639    0.8745; ...
                  0.5941    0.6458    0.6281    0.6826    0.8279    1.0583; ...
                  0.6109    0.6639    0.6281    0.7017    0.8279    1.0873; ...
                  0.6281    0.6658    0.6458    0.7214    0.8509    1.1171; ...
                  0.7017    0.7416    0.7214    0.7836    0.9492    1.2441; ...
                  0.8509    0.8988    0.8745    0.9492    1.1476    1.4990];
    
	figure;
    for i = 1:length(LineColor)
        subplot(1,2,1)
        plot(testlevels, raw3_j2008(i,:), LineColor{i},'LineWidth',LineWidth(i)), hold on, grid on
        if i==1;
            title('C. test 1-kHz tones');
            xlabel('Standard level [dB]')
            ylabel('\Delta at threshold [dB]')
            xlim([minL maxL])
            ylim([minY maxY])
            set(gca,'XTick',testlevels);
            set(gca,'YTick',[minY:0.1:maxY]);
        end
        subplot(1,2,2)
        plot(testlevels, raw4_j2008(i,:), LineColor{i},'LineWidth',LineWidth(i)), hold on, grid on
        if i==1;
            title('D. test BBNs');
            xlabel('Standard level [dB]')
            ylabel('\Delta at threshold [dB]')
            xlim([minL maxL])
            ylim([minY maxY])
            set(gca,'XTick',testlevels);
            set(gca,'YTick',[minY:0.1:maxY]);
        end
        Leg{i} = ['\sigma = ' num2str(sigmas(i))];
    end
    legend(Leg,'Location','NorthWest')
    hFig4(end+1) = gcf;
    
    %% 4.b Results of forward-masking experiment with tones
    testlevels = [  40 60 70 80; ...
                    60 70 80 85];
	minL = min(min(testlevels));
    maxL = max(max(testlevels));
    
	% Obtained on 13-oct-2015:
    raw1_d1997 = [17.5000   32.5000   45.0000   59.0000; ...	% on-freq,  0-ms, longer interval
                  17.5000   32.5000   41.2500   51.2500];       % on-freq, 30-ms, longer interval

	raw2_d1997 = [11.2500   22.5000   33.7500   38.7500; ...	% off-freq,   0-ms, longer interval 
                   6.2500    8.7500   15.0000   17.5000];       % off-freq,  30-ms, longer interval     

    raw1_j2008 = [18.75 24.25 28.125 33.75; ... % 	on-freq,  0-ms, longer interval
                  20.00 26.25 28.750 31.25];    %	on-freq, 30-ms, longer interval];
    % 20.0000   26.2500   28.7500   31.2500	on-freq, 30-ms, longer interval
    % 20.0000   26.2500   28.7500   31.3750	on-freq, 30-ms, longer interval (run 2)

    raw2_j2008 = [21.25 27.5 35.0 45.875; ... %	off-freq,  0-ms, longer interval
                  10.00 16.25 25  28.75];     % off-freq, 30-ms, longer interval
    % 60	70	80	85
    % 21.2500   27.5000   35.0000   46.2500	off-freq,  0-ms, longer interval
    % 21.2500   27.5000   35.0000   45.5000	off-freq,  0-ms, longer interval (run 2)
    
    % ref1_d1997;
    % ref2_d1997;
    % ref1_j2008;
    % ref2_j2008;
    
    figure;
    subplot(1,2,1)
    plot(testlevels(1,:), raw1_d1997(1,:),'r>-','LineWidth',2); hold on
	plot(testlevels(1,:), raw1_d1997(2,:),'ro--');
    plot(testlevels(1,:), raw1_j2008(1,:),'k>-');
	plot(testlevels(1,:), raw1_j2008(2,:),'ko--'); grid on
    xlabel('Masker level [dB SPL]'); xlim([minL maxL])
    ylabel('Threshold [dB SPL]')
    title('A. On-freq.');
    
    subplot(1,2,2)
    plot( testlevels(2,:), raw2_d1997(1,:),'r>-','LineWidth',2); hold on
    plot( testlevels(2,:), raw2_d1997(2,:),'ro--')
    plot( testlevels(2,:), raw2_j2008(1,:),'k>-')
	plot( testlevels(2,:), raw2_j2008(2,:),'ko--'); grid on
    xlabel('Masker level [dB SPL]'); xlim([minL maxL])
    ylabel('Threshold [dB SPL]')
    title('B. Off-freq.');
    
    legend('PEMO,  0-ms', ...
           'PEMO, 30-ms', ...
           'CASP,  0-ms', ...
           'CASP, 30-ms','Location','NorthWest');
    hFig4(end+1) = gcf;
       
	Save_all_figures(hFig4,dir_out,count_saved_figures);
    count_saved_figures = count_saved_figures + length(hFig4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart5
    
    erbc2analyse    = freqtoaud([800],'erb'); % 14 for 1000 Hz (approx.) 
    opts            = il_get_freqs(erbc2analyse,opts);
        
    fnamesM = {[dir_where7 'jepsen2008-fig7-800-Hz-tone-masker-60-dB-dur-10-s.wav']}; % on-frequency masker
    fnames  = {[dir_where7 'jepsen2008-fig7-800-Hz-tone-60-dB-dur-250-ms-onset-200-ms.wav']}; 
   
	testlevels      = [40 60 70 80]; % levels for on-freq
                   
	if opts.nAnalyser == 101 | opts.nAnalyser == 103
        opts.resample_intrep = 'resample_intrep';
    end
    opts.Gain4supra =  10; % 10 dB above the masker level 
    opts.audio.fs   =  fs;
    
    opts.StepdB = 10; 
    opts.StepdBmin = 1;
    
    opts.DurRamps  = 2; % [ms] additional cosine ramps
    opts.bUseRamp  = 1; % additional cosine ramps
    opts.bUseRampS = 0; % additional cosine ramps
    opts.bAddSilence2noise = 1;
    opts.Silence2noise = 500e-3; % 200 + silence = 700
    opts.increment_method = 'level';
    
    opts.Ntimes = 8;
    opts.Nsim = 1;
    opts.bDebug = 1;
    %%% k = 1: offset-onset of  0 ms
    k = 1;
    opts.filename2 = fnames{k};
    
    for i = 1:length(testlevels) 
        for j = 1
            opts.filename1 = fnamesM{j}; % 1 is on-freq
            opts.Gain2file1 = testlevels(j,i)-60;
            opts.Gain2file2 = opts.Gain2file1;
            tmp = AMTControl_cl(opts);
     
            TTh(j,i) = median(tmp.Threshold)+testlevels(j,i);
            
            disp('')
        end
    end
    
    disp('')
end

if bPart6
    
    testdurs = [10 20 40 70 150]*1e-3;
    fmaskerbuf  = [dir_where 'masker-buf-j2008.wav'];
    fsignals{1} = [dir_where 'sig1-2kHz.wav'];
    fsignals{2} = [dir_where 'sig2-2kHz.wav'];
    fsignals{3} = [dir_where 'sig3-2kHz.wav'];
    fsignals{4} = [dir_where 'sig4-2kHz.wav'];
    fsignals{5} = [dir_where 'sig5-2kHz.wav'];
    
    try
        Wavread(fmaskerbuf{1});
        bCreate = 0;
    catch
        bCreate = 1;
    end
    
    if bCreate
        
        nFig = 4;
        [masker,insig] = exp_jepsen2008(fs,nFig); 
        Wavwrite(masker,fs,fmaskerbuf);
        
        for i = 1:length(testdurs)
            Wavwrite(insig(:,i),fs,fsignals{i});
        end
        
    end
    
    if opts.nAnalyser == 101 | opts.nAnalyser == 103
        opts.resample_intrep = 'resample_intrep';
    end
    opts.Gain4supra =  10; % 65 dB + 10 = 75 dB
    opts.audio.fs   =  fs;
    
    opts.StepdB    = 8; 
    opts.StepdBmin = 1;
    
    opts.DurRamps  = 10; % [ms] additional cosine ramps
    opts.bUseRamp  = 1; % additional cosine ramps
    opts.bUseRampS = 0; % additional cosine ramps
    opts.bAddSilence2noise = 0;
    opts.increment_method = 'level';
    
    erbc2analyse    = freqtoaud([1000 4000],'erb');
    opts            = il_get_freqs(erbc2analyse,opts);
     
    opts.Ntimes = 25;
    opts.Nsim   = 4;
    opts.bDebug = 0;
    
    opts.filename1 = fmaskerbuf;
        
    for i = 1:length(testdurs) 
        
        opts.filename2 = fsignals{i};
        % opts.Gain2file1 = testlevels(j,i)-60;
        % opts.Gain2file2 = opts.Gain2file1;
        tmp = AMTControl_cl(opts);

        TTh(i) = median(tmp.Threshold)+65;
        Std75(i) = prctile(tmp.Threshold,75)+65;
        Std25(i) = prctile(tmp.Threshold,25)+65;
        
        disp('')
        
    end
    
    disp('')
    % Only 2000-Hz channel:
    % Dau1997   : TTh = 58.0000   56.0000   54.5000   51.5000   51.5000
    % Jepsen2008: TTh = 61.0000   55.5000   51.5000   53.0000   52.5000
    
%%%       
% Jepsen2008 (8 avgs):    
% TTh   = 79.2500   70.7500   65.7500   65.2500   64.2500
% Std75 = 79.7500   71.0000   66.0000   65.5000   64.7500
% Std25 = 78.7500   70.2500   65.2500   64.5000   64.0000

% (25 avgs., 10-ms ramps)
% TTh   = 66.2500   60.5000   58.5000   58.5000   55.7500
% Std75 = 67.5000   60.7500   58.7500   59.5000   56.7500
% Std25 = 65.2500   60.2500   58.2500   57.2500   55.5000     


% Dau1997 (8 avgs): % this is with a 5-ms ramp
% TTh   = 69.2500   67.5000   66.0000   67.0000   64.7500
% Std75 = 70.0000   68.2500   66.5000   67.0000   65.2500
% Std25 = 69.0000   67.0000   65.5000   66.7500   64.5000

% (25 avgs., 10-ms ramps)
% TTh   = 61.2500   58.7500   59.0000   57.2500   58.5000
% Std75 = 61.5000   60.2500   59.0000   57.7500   59.0000
% Std25 = 61.0000   58.0000   58.0000   56.5000   57.7500
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart7
    
    erbc2analyse    = freqtoaud([3500 5000],'erb'); % 14 for 1000 Hz (approx.) 
    opts            = il_get_freqs(erbc2analyse,opts);
        
    fnamesM = {[dir_where7 'jepsen2008-fig7-4000-Hz-tone-masker-60-dB-dur-10-s.wav'], ... % on-frequency masker
               [dir_where7 'jepsen2008-fig7-2400-Hz-tone-masker-60-dB-dur-10-s.wav'] };   % off-frequency masker (NOT USED HERE YET)
	% fnamesM = {[dir_where7 'masker-j2008-60-dB.wav']};   % off-frequency masker
	            
    fnames  = {[dir_where7 'jepsen2008-fig7-4000-Hz-tone-60-dB-dur-250-ms-onset-200-ms.wav'], ...
               [dir_where7 'jepsen2008-fig7-4000-Hz-tone-60-dB-dur-250-ms-onset-230-ms.wav']}; 
   
	testlevels = 80; % [40 60 80]; 
                   
	if opts.nAnalyser == 101 | opts.nAnalyser == 103
        opts.resample_intrep = 'resample_intrep';
    end
    opts.Gain4supra =  10; % 10 dB above the masker level 
    opts.audio.fs   =  fs;
    
    opts.StepdB    = 10; 
    opts.StepdBmin = 1;
    
    opts.DurRamps  = 2; % [ms] additional cosine ramps
    opts.bUseRamp  = 1; % additional cosine ramps
    opts.bUseRampS = 0; % additional cosine ramps
    opts.bAddSilence2noise = 1;
    silences = ([500])*1e-3; %([480 490 495 500:10:530])*1e-3;
    opts.increment_method = 'level';
    
    opts.Ntimes = 25; %8;
    opts.Nsim   = 2;
    opts.bDebug = 0;

    k = 1;
    opts.filename2 = fnames{k};
    
    tic
    for i = 1:length(silences) 
        for j = 1:length(testlevels)
            
            if i == 4
                opts.Nsim = 1;
            end
            opts.Silence2noise = silences(i);
            idxM = 1;
            if idxM == 2
                warning('off-frequency masker being used')
            end
            opts.filename1 = fnamesM{idxM}; % always on-freq masker (choose 2 for off-freq masker)
            opts.Gain2file1 = testlevels(j)-60;
            opts.Gain2file2 = opts.Gain2file1;
            tmp = AMTControl_cl(opts);
     
            TTh(j,i) = median(tmp.Threshold)+testlevels(j);
            Std75(j,i) = prctile(tmp.Threshold,75)+testlevels(j);
            Std25(j,i) = prctile(tmp.Threshold,25)+testlevels(j);
            disp('')
        end
    end
    toc
    
    disp('')
    % Variance: 0.59
    % Script in template/simulations: jepsen2008preproc_multi
    % Elapsed time is 638.396415 seconds. (only 1 simulation per condition)
    % silences = ([480 490 495 500:10:530])*1e-3;
    %    TTh = [  23.2500   31.1250   29.8750   27.5000   27.5000   26.2500   25.0000; ...
    %             55.7500   56.5000   35.7500   30.0000   30.0000   31.2500   30.0000; ...
    %             70.0000   78.0000   54.3750   37.5000   36.2500   36.2500   35.0000];

    % This new simulation was run to obtain Thres. at simultaneous masking conditions:
    % Elapsed time is 1242 seconds
    % silences = ([480 490 495 500])
    % TTh = [ 25.9375   29.5625   29.6250   27.5000
    %         52.8125   56.1875   36.0000   30.0000
    %         73.3750   75.9375   54.6250   37.5000]
    %   
    % Std75 = 
    %    27.9375   29.8125   30.2500    27.5000
    %    55.6250   58.4375   36.2500    30.0000
    %    75.0000   77.8125   55.0000    37.5000
    % 
    % Std25 =
    %    24.6250   29.0625   28.7500    27.5000
    %    49.9375   53.8750   35.6875    30.0000
    %    72.8750   72.8125   54.2500    37.5000
   
% Dau1997, 
% Variance: 0.68
% Script in template   : dau1997preproc
% Elapsed time is 1352 seconds.
% TTh = 28.6875   30.5000   31.6875   26.2500   25.0000   22.5000   20.0000
%       53.1250   54.0625   50.8125   39.4375   40.0000   37.5000   33.7500
%       74.0000   75.6250   72.6250   63.7500   58.7500   56.2500   51.2500
% 
% Std75 = 29.4375   32.7500   32.4375   26.2500   25.0000   22.5000   20.0000
%    54.0000   55.9375   51.2500   39.9375   40.0000   37.5000   33.7500
%    75.3125   78.3125   72.9375   63.7500   58.7500   56.2500   51.2500
% 
% Std25 = 26.5625   28.0625   30.8125   26.2500   25.0000   22.5000   20.0000
%    52.1250   53.1250   50.1250   38.9375   40.0000   37.5000   33.7500
%    73.3750   72.8125   71.7500   63.7500   58.7500   56.2500   51.2500


% Dau1997 (25 Avgs):
% TTh =   25.5000   32.0000   31.2500   26.2500   25.0000   22.5000   20.0000
%         49.3750   53.9375   50.7500   39.7500   40.0000   37.5000   33.7500
%         76.3125   71.5625   71.8750   63.7500   58.7500   56.2500   51.2500
% 
% Std75 = 29.5000   34.0625   31.5625   26.2500   25.0000   22.5000   20.0000
%         53.1250   55.1875   51.3750   40.0000   40.0000   37.5000   33.7500
%         77.2500   74.6875   72.5000   63.7500   58.7500   56.2500   51.2500
% 
% Std25 = 24.0000   30.4375   31.1250   26.2500   25.0000   22.5000   20.0000
%         48.7500   53.2500   49.2500   39.3750   40.0000   37.5000   33.7500
%         73.6875   70.6250   70.8750   63.7500   58.7500   56.2500   51.2500
%
% Jepsen2008 (25 Avgs, 10-ms cos ramps)
% TTh   = 24.0000   29.0000   28.4375   27.5000   27.5000   26.2500   25.0000
%         55.4375   57.0000   36.0000   30.0000   30.0000   31.2500   30.0000
%         76.1875   78.1250   54.9375   37.5000   36.2500   36.2500   35.0000
% 
% Std75 = 26.2500   30.1250   30.0000   27.5000   27.5000   26.2500   25.0000
%         56.9375   57.6875   36.2500   30.0000   30.0000   31.2500   30.0000
%         78.4375   79.3750   55.0000   37.5000   36.2500   36.2500   35.0000
% 
% Std25 = 23.3125   25.8750   27.8125   27.5000   27.5000   26.2500   25.0000
%         54.4375   56.2500   35.6250   30.0000   30.0000   31.2500   30.0000
%         74.9375   74.5625   54.3125   37.5000   36.2500   36.2500   35.0000




% Dau1997, off frequency
% TTh = 5.0000    5.0000    5.0000    5.0000    5.0000    5.0000    5.0000
%       11.2500   11.2500   17.1250   13.7500   10.0000    8.7500    7.5000
%       31.2500   31.2500   37.5625   36.2500   25.0000   21.2500   17.5000
%
% Std75 =  5.0000    5.0000    5.0000    5.0000    5.0000    5.0000    5.0000
%         11.2500   11.2500   17.3750   13.7500   10.0000    8.7500    7.5000
%         31.2500   31.2500   37.6250   36.2500   25.0000   21.2500   17.5000
% 
% Std25 =  5.0000    5.0000    5.0000    5.0000    5.0000    5.0000    5.0000
%         11.2500   11.2500   16.6875   13.7500   10.0000    8.7500    7.5000
%         31.2500   31.2500   37.5000   36.2500   25.0000   21.2500   17.5000
%
% Jepsen2008, off frequency
% 
% TTh = 13.7500   13.7500   15.0000   15.0000   15.0000   15.0000   15.0000
%       24.7500   24.8750   22.5000   21.2500   18.7500   17.5000   17.5000
%       53.7500   52.6250   42.5000   34.2500   28.7500   27.5000   25.0000
% 
% 
% Std75 = 13.7500   13.7500   15.0000   15.0000   15.0000   15.0000   15.0000
%         25.0000   25.0000   22.5000   21.2500   18.7500   17.5000   17.5000
%         54.0000   53.1875   42.5000   34.2500   28.7500   27.5000   25.0000
% 
% Std25 = 13.7500   13.7500   15.0000   15.0000   15.0000   15.0000   15.0000
%         24.3750   24.6250   22.5000   21.2500   18.7500   17.5000   17.5000
%         53.7500   52.5625   42.5000   34.2500   28.7500   27.5000   25.0000

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jepsen 2008, on-frequency but 13 dB: (on 02/11/2015)
% ([490 495 500 510])*1e-3;
% TTh = [73.0000   51.7500   38.1875   40.0000]; 
% Std75 = [74.2500   52.2500   38.7500   40.0000];
% Std25 = [71.7500   51.2500   37.6250   40.0000];
end

% ([490 495 500 510 520 530])*1e-3;
% Jepsen2008 (25 Avgs, 10-ms cos ramps)
% TTh   = 29.0000   28.4375   27.5000   27.5000
%         57.0000   36.0000   30.0000   30.0000

%         78.1250   54.9375   37.5000   36.2500 % Jepsen2008

%         71.5625   71.8750   63.7500   58.7500 % Dau1997

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart8
    
    % jepsen2008:
    %   -var = 0.59, gain_after_drnl = 17.86 dB
    
    toffsetonset = [-20 -10 -5 0 10 20 30]; % corresponding to ([480 490 495 500:10:530])*1e-3;
    testlevels = [40; 60; 80];
    
    TTh_j2008 = [25.9375   29.5625   29.6250   27.5000   27.5000   26.2500   25.0000; ...
                 52.8125   56.1875   36.0000   30.0000   30.0000   31.2500   30.0000; ...
                 73.3750   75.9375   54.6250   37.5000   36.2500   36.2500   35.0000];
    Prct75_j2008 = [27.9375   29.8125   30.2500    27.5000  27.5000   26.2500   25.0000; ...
                    55.6250   58.4375   36.2500    30.0000  30.0000   31.2500   30.0000; ...
                    75.0000   77.8125   55.0000    37.5000  36.2500   36.2500   35.0000];
    
    Prct25_j2008 = [24.6250   29.0625   28.7500    27.5000  27.5000   26.2500   25.0000; ...
                    49.9375   53.8750   35.6875    30.0000  30.0000   31.2500   30.0000; ...
                    72.8750   72.8125   54.2500    37.5000  36.2500   36.2500   35.0000];
            
    errorL_j2008 = TTh_j2008-Prct25_j2008;
    errorU_j2008 = Prct75_j2008-TTh_j2008;
    
     TTh_d1997 =   [28.6875   30.5000   31.6875   26.2500   25.0000   22.5000   20.0000; ...
                   53.1250   54.0625   50.8125   39.4375   40.0000   37.5000   33.7500; ...
                   74.0000   75.6250   72.6250   63.7500   58.7500   56.2500   51.2500];
    Prct75_d1997 = [29.4375   32.7500   32.4375   26.2500   25.0000   22.5000   20.0000; ...
                    54.0000   55.9375   51.2500   39.9375   40.0000   37.5000   33.7500; ...
                    75.3125   78.3125   72.9375   63.7500   58.7500   56.2500   51.2500];
    Prct25_d1997 = [26.5625   28.0625   30.8125   26.2500   25.0000   22.5000   20.0000; ...
                    52.1250   53.1250   50.1250   38.9375   40.0000   37.5000   33.7500; ...
                    73.3750   72.8125   71.7500   63.7500   58.7500   56.2500   51.2500];
    errorL_d1997 = TTh_d1997-Prct25_d1997;
    errorU_d1997 = Prct75_d1997-TTh_d1997;
    
	figure;
    subplot(1,2,1)
    errorbar(toffsetonset,TTh_j2008(3,:),errorL_j2008(3,:),errorU_j2008(3,:),'ko-'); grid on; hold on
    errorbar(toffsetonset,TTh_j2008(2,:),errorL_j2008(2,:),errorU_j2008(2,:),'k>-');
    errorbar(toffsetonset,TTh_j2008(1,:),errorL_j2008(1,:),errorU_j2008(1,:),'ks-'); 
    xlabel('Offset-onset interval [ms]');
    ylabel('Masked threshold [dB SPL]')
    legend(sprintf('M at %.0f dB',testlevels(1)), ...
           sprintf('M at %.0f dB',testlevels(2)), ...
           sprintf('M at %.0f dB',testlevels(3)) );
    ha = gca;
    title('Forward masking (jepsen2008)');
    
    subplot(1,2,2)
    errorbar(toffsetonset,TTh_d1997(3,:),errorL_d1997(3,:),errorU_d1997(3,:),'ro-'); grid on; hold on
    errorbar(toffsetonset,TTh_d1997(2,:),errorL_d1997(2,:),errorU_d1997(2,:),'r>-');
    errorbar(toffsetonset,TTh_d1997(1,:),errorL_d1997(1,:),errorU_d1997(1,:),'rs-'); 
    xlabel('Offset-onset interval [ms]');
    ylabel('Masked threshold [dB SPL]')
    legend(sprintf('M at %.0f dB',testlevels(1)), ...
           sprintf('M at %.0f dB',testlevels(2)), ...
           sprintf('M at %.0f dB',testlevels(3)) );
    ha(end+1) = gca;
    title('Forward masking (dau1997)');
    
    linkaxes(ha,'xy');
    
    disp('')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF

function opts = il_get_freqs(erbc2analyse, opts)

if length(erbc2analyse) == 1
    opts.fc2plot_idx    = round(erbc2analyse(1))-2;
    opts.fc2plot_idx2 = opts.fc2plot_idx;
else
    opts.fc2plot_idx    = ceil(erbc2analyse(1))-2;
    opts.fc2plot_idx2   = floor(erbc2analyse(end))-2;
end
