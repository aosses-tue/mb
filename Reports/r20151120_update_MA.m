function r20151120_update_MA
% function r20151120_update_MA
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%       See also r20151106_update_MA
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 20/11/2015
% Last update on: 20/11/2015 
% Last use on   : 20/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

% f1 = 'meas-ac-4-dist-ane-HP.wav';
f1 = 'meas-ac-2-dist-ane.wav';

dirout = [Get_TUe_paths('outputs') 'lx20151111' delim];
Mkdir(dirout);

bSave = 0;
bPart1 = 0; % Creates HWN Hummer weighted noise
bPart2 = 1;
bPart3 = 0;
bPart4 = 0; % comparing DRNL with Gammatone output...
bPart5 = 1;

if bPart1
    
    dir_files = [Get_TUe_paths('lx_Text') 'lx2015-10-30-update-to-music-instruments' delim 'Audio' delim];
    
    files_in  = {   'meas-ac-2-dist-ane-HP+LP.wav', ... % Measured Hummer, ac mode 2
                    'meas-ac-4-dist-ane-HP+LP.wav'};    % Measured Hummer, ac mode 4
                
    files_HWN = {   'hummer_noise_ac2_20151112.wav', ... % Hummer-weighted noise, ac mode 2
                    'hummer_noise_20151104.wav'};        % Hummer-weighted noise, ac mode 4
    for i = 1:length(files_HWN)
        [x fs] = Wavread([dir_files files_in{i}]);
        [Pxx(:,i)  F] = Get_PSD_analysis_arrays(x,fs,20,0); % it pools buffered signals and gets averaged power
        
        [y fs] = Wavread([dir_files files_HWN{i}]);
        [Pxxn(:,i) F] = Get_PSD_analysis_arrays(y,fs,20,0); % it pools buffered signals and gets averaged power
    end
    
    figure;
    plot(F,10*log10(Pxx)); grid on; hold on
    plot(F,10*log10(Pxxn),'LineWidth',2); grid on
    xlim([350 1800])
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');

    % title('Spectra of the Hummer sounds and the ''new'' noise');
    legend({'ac-mode-2';'ac-mode-4'});
    Save_figure_as(gcf,[dirout 'HWN-fig'],'epsc');
    
end

if bPart2
    
    close all;
    dir2check = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Hummer\demo-HWN+meas-sim-hummer\Results-test-1\';
    
    %%%
    opts = [];
    opts.N4mean = 4;
    opts.mode = 'median';
    opts.filter = 'hummer-ac-2-*T0*';
    opts.bPlot = 1;
    [SRT2 Stdev2 fns] = quick_staircases(dir2check,opts);
    
    opts.filter = 'hummer-ac-2-*T2*';
    [SRT2(3,1) Stdev2(3,1) fns] = quick_staircases(dir2check,opts);
    
    opts.filter = 'hummer-ac-2-*T3*';
    [SRT2(4:5,1) Stdev2(4:5,1) fns] = quick_staircases(dir2check,opts);
       
    %%%
    opts = [];
    opts.N4mean = 4;
    opts.mode = 'median';
    opts.filter = 'hummer-ac-4-*T0*';
    opts.bPlot = 1;
    [SRT4 Stdev4 fns] = quick_staircases(dir2check,opts);
    
    opts.filter = 'hummer-ac-4-*T2*';
    [SRT4(3,1) Stdev4(3,1) fns] = quick_staircases(dir2check,opts);
    
    opts.filter = 'hummer-ac-4-*T3*';
    [SRT4(4:5,1) Stdev4(4:5,1) fns] = quick_staircases(dir2check,opts);
        
    SRT2avg = [mean(SRT2(1:2)) SRT2(3) mean(SRT2(4:5))];
    SRT4avg = [mean(SRT4(1:2)) SRT4(3) mean(SRT4(4:5))];
    
    GreenNice = [0 0.5 0];
    offsetx = 0.1;
    figure;
    stem([1:3]-offsetx,SRT2avg,'bo','LineWidth',2); hold on
    stem([1:3]+offsetx,SRT4avg,'s','Color',GreenNice,'LineWidth',2)
    grid on
    
    subjects = {'AO','AS','GS'};
    
    ymin = min(min([SRT2avg SRT4avg]));
    ymax = max(max([SRT2avg SRT4avg]));
    ystep = 2;
    
    ylim([ymin-2 ymax+2])
    set(gca,'YTick',[ymin:ystep:ymax]);
    set(gca,'XTick',1:3);
    set(gca,'XTickLabel',subjects);
    xlabel('Subject')
    ylabel('SNR (70.7%-correct point) [dB]')
    legend('ac-mode-2','ac-mode-4')
    
    Save_figure_as(gcf,[dirout 'exp-res'],'epsc');
end

if bPart3
    mu = 0;
    sigma = [0.6 1 2 5 10 20 50 100];
    
    load([Get_TUe_paths('outputs') 'template-20151104.mat']);
    load([Get_TUe_paths('outputs') 'calc-20151104.mat']);
    
    idx1 = 8500;
    idx2 = 28500;
    sigint1   = sigint1(idx1:idx2,:);
    sigint1s2 = sigint1s2(idx1:idx2,:);
    sigint2   = sigint2(idx1:idx2,:);
    template = template(idx1:idx2,:);
    intM = intM(idx1:idx2,:);
    xx = xx(idx1:idx2,:);
    yy = yy(idx1:idx2,:);
    zz = zz(idx1:idx2,:);
    
    N = size(template,1);
    Ntimes = 100;
    
    % intM
    intN1 = sigint1 - xx;
    intN2 = sigint1s2 - yy;
    intSN = sigint2 - zz;
    
    for k = 1:length(sigma)
        
        for i = 1:Ntimes

            noise = normrnd(mu,sigma(k),N,6);
            intN1t = intN1 + noise(:,1); 
            intN2t = intN2 + noise(:,2);
            intSNt = intSN + noise(:,3);
            intMt1 = intM + noise(:,4);
            intMt2 = intM + noise(:,5);
            intMt3 = intM + noise(:,6);

            decision(i,1) = optimaldetector( intN1t-intMt1 ,template,fs_intrep);
            decision(i,2) = optimaldetector( intN2t-intMt2,template,fs_intrep);
            % decision(i,3) = optimaldetector( intSNt-intMt3,template,fs_intrep); 
        end
        sig(k,:) = prctile(decision(:),[75 95]); %[prctile(decision(:),25) prctile(decision(:),50) prctile(decision(:),75)];
        
    end
        
    disp('')
end

if bPart4
    
    dirs1 = 'D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2015-10-30-update-to-music-instruments\Audio-drnl\';
    dirs2 = 'D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2015-10-30-update-to-music-instruments\Audio-gamma\';
    
    files = {   'sum-audfilter-fc-3983-Hz.wav'; ...
                'file1-audfilter-fc-3983-Hz.wav'; ...
                'file2-audfilter-fc-3983-Hz.wav'}
    for i = 1:3
        [y(:,i) fs] = Wavread([dirs1 files{i}]); 
    end
    y = gaindb(y,0); % approx. 2*17 dB
    
    files = {   'sum-audfilter-fc-3990-Hz.wav'; ...
                'file1-audfilter-fc-3990-Hz.wav'; ...
                'file2-audfilter-fc-3990-Hz.wav'}
    for i = 1:3
        x(:,i) = Wavread([dirs2 files{i}]);
    end
    
    t = (1:size(x,1))/fs;
    
    close all
    figure;
    plot(t,x(:,1),t,y(:,1)+0.08), grid on
    legend('gamma','drnl');
    
    ysign = sign(y(:,1));
    yexpansion = gaindb( ysign.*y(:,1).*y(:,1),25 );
    
    figure;
    plot(t,x(:,1),t,y(:,1)+0.08), grid on
    legend('gamma','drnl');
    
    figure;
    plot(t,yexpansion,t+0.05,x(:,1)), grid on
    % legend('gamma','drnl');
    legend('drnl','gamma');
end

if bPart5
    
    
    
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Taken from: r20151016_update_dau1997_jepsen2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_where  = [Get_TUe_paths('outputs') 'audio-20151006' delim];
dir_where7 = [Get_TUe_paths('outputs') 'audio-20151006-forward' delim];
dir_out    = [Get_TUe_paths('lx_Text') 'lx2015-10-16-decision-CASP' delim 'Figures-new' delim];

if nargin == 0
    %         1 2 3 4 5 6 7 8
    bParts = [0 0 1 0 0 0 1 0];
    % CC:     0 1 0 - - 1 0 -
    %                   0       % if 2000-Hz  
end

fs = 44100;


opts.nAnalyser      = 103; % 101 = modfilterbank, 103 - jepsen2008
opts.bDecisionMethod = 5; % 2 - cc; 4 - dprime; 5 - cc updated

switch opts.nAnalyser
    
    case 101
        opts.sigma   = 0; % 0.8; % 0.68;    
                                
        opts.modfiltertype = 'dau1997wLP';

    case {103, 104}
        opts.sigma   = 0; %0.88; 
    
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
opts.Nreversals     =  10;
opts.nDown          =  2;
% opts.experiment_type = 'AFC'; % opts.experiment_type = 'constant';
opts.experiment_type = 'constant';
opts.Ntimes     = 1;
opts.Nsim       = 1;

bDebug = 0;
opts.bDebug = bDebug;

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
    
    erbc2analyse    = freqtoaud([4000],'erb'); % 14 for 1000 Hz (approx.)  % 
    opts            = il_get_freqs(erbc2analyse,opts);
        
    fnamesM = {[dir_where7 'jepsen2008-fig7-4000-Hz-tone-masker-60-dB-dur-10-s.wav'], ... % on-frequency masker
               [dir_where7 'jepsen2008-fig7-2400-Hz-tone-masker-60-dB-dur-10-s.wav'] };   % off-frequency masker
    fnames  = {[dir_where7 'jepsen2008-fig7-4000-Hz-tone-60-dB-dur-250-ms-onset-200-ms.wav'], ...
               [dir_where7 'jepsen2008-fig7-4000-Hz-tone-60-dB-dur-250-ms-onset-230-ms.wav']}; 
   
	testlevels      = 60;
    
	if opts.nAnalyser == 101 | opts.nAnalyser == 103 | opts.nAnalyser == 104
        opts.resample_intrep = 'resample_intrep';
    end
    
    opts.Gain4supra=  10; % 10 dB above the masker level 
    
    opts.audio.fs  =  fs;
    
    opts.DurRamps  = 2; % [ms] additional cosine ramps
    opts.bUseRamp  = 1; % additional cosine ramps
    opts.bUseRampS = 0; % additional cosine ramps
    opts.bAddSilence2noise = 1;
    opts.Silence2noise = 500e-3; % 200 + silence = 700. If sil = 400e-3, then signal dur = 300e-3 -> simultaneous masking
    opts.increment_method = 'level';
    
    opts.Ntimes = 1;
    opts.Nsim = 1;
    opts.bDebug = 0;
    %%% k = 1: offset-onset of  0 ms
    %%% k = 2: offset-onset of 30 ms
    for k = 1 % :2
        opts.filename2 = fnames{k};

        for j = 1% :2

            opts.filename1 = fnamesM{j}; % 1 is on-freq, 2 is off-freq
            opts.Gain2file1 = testlevels-60;
            opts.Gain2file2 = opts.Gain2file1;
            [tmp outsIR] = AMTControl_cl(opts);

           disp('')

        end
    end
    
    switch opts.nAnalyser
        case 101
            L = length(outsIR.out_SN_interval);
            error('')
            Lreal = round(L/3/12);
            outsIR.out_SN_interval = reshape(outsIR.out_SN_interval,round(L/3),3);
            
            outsIR.out_SN_interval = outsIR.out_SN_interval(1:Lreal,:);
            t = ( 1:length(outsIR.out_SN_interval(:,3)) )/(fs/4);
            
        case 103
            NCh = 1;
            L = length(outsIR.out_SN_interval);
            error('')
            Lreal = round(L/NCh/12);
            outsIR.out_SN_interval = reshape(outsIR.out_SN_interval,round(L/NCh),NCh);
            
            outsIR.out_SN_interval = outsIR.out_SN_interval(1:Lreal,:);
            t = ( 1:length(outsIR.out_SN_interval(:,NCh)) )/(fs/4);
    end
    
    figure;
    plot(t,outsIR.out_SN_interval(:,NCh)); grid on
    title(sprintf('Analyser: %.0f',opts.nAnalyser))
    xlabel('Time [s]')
    ylabel('Output [MU]')
    
    filename = sprintf('%sS+N-IR-%.0f-no-limit',dirout,opts.nAnalyser);
    Saveas(gcf,filename,'epsc');
    
    disp('')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nieuwe resultaten op 6/11
    %% Dau1997 (first row=on-freq; second row=off-freq):
    % TTh_00ms =    27.5000   35.0000   48.0000   62.5000
    %               14.2500   25.5000   36.2500   41.7500
    % Prct75_00ms = 28.0000   38.0000   50.7500   63.0000
    %               14.7500   26.0000   36.5000   42.0000
    % Prct25_00ms = 25.7500   33.7500   44.7500   61.0000
    %               13.5000   25.0000   35.7500   41.2500
    % 
    % TTh_30ms =    20.0000   34.0000   42.2500   51.5000
    %                7.5000   12.0000   17.2500   20.5000
    % Prct75_30ms = 20.2500   34.2500   42.5000   51.7500
    %                8.0000   12.0000   18.2500   21.5000
    % Prct25_30ms = 19.5000   34.0000   41.7500   51.2500
    %                7.0000   11.5000   17.0000   19.7500
    %----------------------------------------------------------------------
    %% Jepsen2008
    % TTh_00ms =    29.0000   31.0000   32.0000   37.7500
    %               21.0000   26.7500   34.2500   43.0000
    % Prct75_00ms = 29.2500   31.5000   33.0000   38.5000
    %               21.5000   27.2500   35.2500   43.2500
    % Prct25_00ms = 28.7500   29.7500   31.7500   37.5000
    %               21.0000   26.5000   33.2500   42.0000
    % 
    % TTh_30ms    = 26.0000   29.7500   32.2500   35.5000
    %               17.2500   20.0000   25.0000   29.0000
    % Prct75_30ms = 26.0000   30.0000   32.5000   36.5000
    %               17.7500   20.5000   25.7500   29.2500 
    % Prct25_30ms = 26.0000   29.0000   32.0000   34.7500
    %               16.5000   19.7500   24.7500   29.0000
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end


if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF

function il_get_avg_power(file_mat)

try
    load(file_mat); 
catch
    error('Go to r20151106_update_MA')
end

B = fir2(2^11,F/22050,sqrt(Pin));       %% sqrt because Pxx is energy
w = randn(11*44100,1);

y = filter(B,1,w);

Pxxn = Get_PSD_analysis_arrays(y,fs,20,0); % it pools buffered signals and gets averaged power


function opts = il_get_freqs(erbc2analyse, opts)

if length(erbc2analyse) == 1
    opts.fc2plot_idx    = round(erbc2analyse(1))-2;
    opts.fc2plot_idx2 = opts.fc2plot_idx;
else
    opts.fc2plot_idx    = ceil(erbc2analyse(1))-2;
    opts.fc2plot_idx2   = floor(erbc2analyse(end))-2;
end