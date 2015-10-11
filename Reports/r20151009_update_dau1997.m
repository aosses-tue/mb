function r20151009_update_dau1997(bParts)
% function r20151009_update_dau1997(bParts)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 06/10/2015
% Last update on: 06/10/2015 
% Last use on   : 06/10/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;

dir_where = [Get_TUe_paths('outputs') 'audio-20151006' delim];

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
    bParts = [1 0 1 0 0 0 1];
end

fs = 44100;

bPart1 = bParts(1); % Calibration of the model: sigma = 1.85 for bDecisionMethod = 2
bPart2 = bParts(2); % Simultaneous masking, deterministic
bPart3 = bParts(3); % Simultaneous masking, stochastic

opts.nAnalyser      = 101; % 99 - dau1996a, 99.1, 100 - dau1996, my template estimation

opts.bDecisionMethod = 2; % 2 - cc; 4 - dprime
switch opts.bDecisionMethod
    case 2
        switch opts.nAnalyser
            case 99
                opts.sigma   = 3.25; % tone: 3.25 = band 13-15; 2.7 = band 13-14; 1.9 = band 14; 
            case 100
                opts.sigma   = 1.35; % tone: 3 = band 13-15; % 2.45 = band 13-14; % 1.72 = band 14;
                                     %   BW:                                        1.35 = band 14; (target = 0.83 = -2 dB);
            case 101
                opts.sigma   = 1.08;    % 11/10/2015, BBN: bands 1-30 = 1.08; band 14 = 0.95 (target = 0.83 dB)
                                        % old -- tone: 1.23 = band 14 (all); 1.30 = band 14 (1,2); 1.35 = band 14 (1)
                opts.modfiltertype = 'dau1997wLP';
            case 103
                opts.sigma   = 21;  %   BW:  21 = band 2-33 (all); 1.95 = band 14; (target = 0.83 = -2 dB);
        end
    case 3
        opts.sigma   = 0.85; % dprime NOT GIVING RELIABLE RESULTS
    case 4
        opts.sigma   = 1.12; % dprime
end
    
opts.var            = opts.sigma.*opts.sigma;
opts.audio.fs       = fs;
opts.MethodIntRep   = 1; % 1 - my method; 2 - using casptemplate.m

opts.Reversals4avg  =  6; % 10

opts.do_template    =  1;
opts.do_simulation  =  1;
opts.Nreversals     =  8;

erbc2analyse        = freqtoaud([500 2000],'erb'); % 14 for 1000 Hz (approx.)  
opts.fc2plot_idx    = ceil(erbc2analyse(1))-2;
if length(erbc2analyse) > 1
    opts.fc2plot_idx2 = floor(erbc2analyse(end))-2;
else
    opts.fc2plot_idx2 = opts.fc2plot_idx;
end

bDebug = 0;
opts.bDebug = bDebug;

count_saved_figures = 1;
if bPart1
    % To do: save template
    bTones = 1;
    bBBN = ~bTones;
    if bTones
        opts.DurRamps   = 125; % additional cosine ramps
        opts.bUseRamp   = 1; % additional cosine ramps
        opts.bUseRampS  = 1; % additional cosine ramps
        opts.filename1 = [Get_TUe_paths('outputs') 'AMTControl-examples' delim 'tone-f-1000-Hz-at-60-dB-dur-800-ms.wav'];
        opts.filename2 = [Get_TUe_paths('outputs') 'AMTControl-examples' delim 'tone-f-1000-Hz-at-42-dB-dur-800-ms.wav'];
    end
    if bBBN
        opts.DurRamps   = 0; % additional cosine ramps
        opts.bUseRamp   = 0; % additional cosine ramps
        opts.bUseRampS  = 0;
        opts.filename1 = [Get_TUe_paths('outputs') 'audio-20151006' delim 'jepsen2008-BW-at-60-dB-dur-500-ms.wav'];
        opts.filename2 = [Get_TUe_paths('outputs') 'audio-20151006' delim 'jepsen2008-BW-at-42-dB-dur-500-ms.wav'];    
    end
    
    opts.Gain4supra = 5; % dB
    opts.audio.fs   = fs;
    
    opts.StepdB = 2; 
    opts.StepdBmin = 0.2;
    
    refSPL  = 60; % [20 30 40 50 60 70 80];
    
    for i=1:length(refSPL)
        testSPL = refSPL(i)-18; % 42 dB for 60 dB    
        opts.Gain2file1 = refSPL(i) - 60;
        opts.Gain2file2 = refSPL(i) - 60;
        tr_tmp = AMTControl_cl(opts);

        th_JND(i) = sum_dB_arit([refSPL(i) testSPL+tr_tmp.Threshold]) - refSPL(i);
    end
    
    if nargout > 0
        th = th_JND; 
        return;
    end
end

close all
opts.StepdB    = 4; 
opts.StepdBmin = 1;
if bPart2
    
    opts.DurRamps   = 200; % additional cosine ramps
    opts.Gain4supra =  -3; % when creating: lvl = 75 dB; lvl supra = 75 dB + Gain4supra
    opts.audio.fs   =  fs;
    
    opts.bUseRamp  = 1; % additional cosine ramps
    opts.bUseRampS = 1; % additional cosine ramps
    opts.increment_method = 'modulation-depth';
    
    if strcmp(opts.increment_method,'modulation-depth')
        opts.fmod = 120; % Hz
        opts.dur_test = 1; % s
    end

    opts.Ntimes = 5;
    opts.Nsim = 5;
    opts.filename1 = [dir_where 'NBN-fc-5000-Hz-BW-31-Hz-65-dB.wav'];
    tmp = AMTControl_cl(opts);
    
    if nargout > 0
        th = th_det;
        return;
    end
end

%%%
if bPart3
    
    fnames = {[dir_where 'NBN-fc-5000-Hz-BW-10-Hz-65-dB.wav'], ...
              [dir_where 'NBN-fc-5000-Hz-BW-100-Hz-65-dB.wav'], ...
              [dir_where 'NBN-fc-5000-Hz-BW-250-Hz-65-dB.wav'], ...
              [dir_where 'NBN-fc-5000-Hz-BW-500-Hz-65-dB.wav'], ...
              [dir_where 'NBN-fc-5000-Hz-BW-800-Hz-65-dB.wav'], ...
              [dir_where 'NBN-fc-5000-Hz-BW-1000-Hz-65-dB.wav'], ...
              [dir_where 'NBN-fc-5000-Hz-BW-2500-Hz-65-dB.wav'], ...
              [dir_where 'NBN-fc-5000-Hz-BW-5000-Hz-65-dB.wav'], ...
              [dir_where 'NBN-fc-5000-Hz-BW-10000-Hz-65-dB.wav']}
   
    opts.DurRamps   =  50; % additional cosine ramps
    opts.Gain4supra =  0; 
    opts.audio.fs   =  fs;
    
    opts.bUseRamp  = 1; % additional cosine ramps
    opts.bUseRampS = 1; % additional cosine ramps
    opts.increment_method = 'modulation-depth';
    
    if strcmp(opts.increment_method,'modulation-depth')
        opts.fmod = 5; % Hz
        opts.dur_test = 0.5; % s
    end

    opts.Ntimes = 25;
    opts.Nsim = 5;
    
    for i = 1:length(fnames)
        opts.filename1 = fnames{i};
        tmp = AMTControl_cl(opts);
        
        TTh(i) = median(tmp.Threshold);
        disp('')
    end
    figure; plot(TTh);
    
    if nargout > 0
        th = th_det;
        return;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
