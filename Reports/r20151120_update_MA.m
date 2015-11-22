function r20151120_update_MA
% function r20151120_update_MA
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
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
bPart2 = 1; % Plot results of hummer pilot

bPart3 = 0;
bPart4 = 0; % comparing DRNL with Gammatone output...
bPart5 = 1;

if bPart1
    
    dir_files = [Get_TUe_paths('lx_Text') 'lx2015-11-20-update-3AFC' delim 'Audio' delim];
    
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
    dir2check = [Get_TUe_paths('ex_Experiments') '2015-APEX-my-experiments' delim ...
                 'Hummer' delim 'demo-HWN+meas-sim-hummer' delim 'Results-test-1' delim];
    
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