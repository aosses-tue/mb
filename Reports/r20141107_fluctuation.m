function [FSAM, FSFM, h, afiles] = r20141107_fluctuation(N_blocks)
% function [FSAM, FSFM, h, afiles] = r20141107_fluctuation(N_blocks)
%
% 1. Description:
%       Implements the Fluctuation Strength model as reported in the update 
%       on 14/11/2014.
% 
% 2. Stand-alone example:
%       r20141107_fluctuation;
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 05/11/2014
% Last update on: 05/11/2014 % Update this date manually
% Last use on   : 18/11/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO:
%
% Correct displayed values at the end (see 'output').

if nargin == 0
    N_blocks = 1;
end

pathaudio   = Get_TUe_paths('outputs');
pathfigures = [Get_TUe_paths('outputs') 'Figures-20141114' delim];
Mkdir(pathfigures);

filename    = [pathaudio 'ref_fluct'];

N           = 8192*6;
bDebug      = 1;
count_afiles = 1;

try
    
    Wavread(filename);
    
catch % if reference does not exist, then it is created...
    
    opts.bDoFluct = 1;
    opts.bDoRough = 0;
    opts.dur = (20*N/44100); %  approx. 20 times 200e-3 ms;
    opts.bGen_test_tones = 1;
    
    opts.bPsySound = 1;
    opts.bDoRamp = 0; 
    opts.bDoZeroPadding = 0;
    
    outs = Generate_reference_sounds_Zwicker2007(opts);
end

filenamesAM = { 'test_fluct_fc_1000_AM_m_100_fmod_001Hz_70_dBSPL', ...
                'test_fluct_fc_1000_AM_m_100_fmod_002Hz_70_dBSPL', ...
                'test_fluct_fc_1000_AM_m_100_fmod_004Hz_70_dBSPL', ...
                'test_fluct_fc_1000_AM_m_100_fmod_008Hz_70_dBSPL', ...
                'test_fluct_fc_1000_AM_m_100_fmod_016Hz_70_dBSPL', ...
                'test_fluct_fc_1000_AM_m_100_fmod_032Hz_70_dBSPL'};

% filenamesFM = { 'test_fluct_fc_1000_FM_dev_700_fmod_001Hz_60_dBSPL', ...
%                 'test_fluct_fc_1000_FM_dev_700_fmod_002Hz_60_dBSPL', ...
%                 'test_fluct_fc_1000_FM_dev_700_fmod_004Hz_60_dBSPL', ...
%                 'test_fluct_fc_1000_FM_dev_700_fmod_008Hz_60_dBSPL', ...
%                 'test_fluct_fc_1000_FM_dev_700_fmod_016Hz_60_dBSPL', ...
%                 'test_fluct_fc_1000_FM_dev_700_fmod_032Hz_60_dBSPL'};

filenamesFM = { 'test_fluct_fc_1500_FM_dev_700_fmod_001Hz_70_dBSPL', ...
                'test_fluct_fc_1500_FM_dev_700_fmod_002Hz_70_dBSPL', ...
                'test_fluct_fc_1500_FM_dev_700_fmod_004Hz_70_dBSPL', ...
                'test_fluct_fc_1500_FM_dev_700_fmod_008Hz_70_dBSPL', ...
                'test_fluct_fc_1500_FM_dev_700_fmod_016Hz_70_dBSPL', ...
                'test_fluct_fc_1500_FM_dev_700_fmod_032Hz_70_dBSPL'};
            
for k = 1:length(filenamesAM)
    
    filename = [pathaudio filenamesAM{k} '.wav'];
    
    [x Fs] = Wavread(filename);
    
    afiles{count_afiles} = filenamesAM{k};
    count_afiles = count_afiles + 1;
    
    %%%%
    Fsnew = 44100;
    x = resample(x,Fsnew,Fs);
    Fs = Fsnew;
    %%%%
    
    starti  = 1;
    endi    = starti + N-1 + 2048*(N_blocks-1); % 5 additional analysis blocks
    
    insig = x(starti:endi);
    t = ( 1:length(insig) )/Fs;     

    if bDebug % Similar to Fig. 7.2
        h(1) = figure(1);
        subplot(6,1,k)
        plot(t,insig);
        xlabel('Time [s]')
        ylabel('Amplitude')
        title(sprintf('Sig: %s, %.1f [dB]',name2figname(filenamesAM{k}(12:end)),rmsdb(insig)+90))
        grid on
        
        optsfig.format = 'epsc';
        Saveas(gcf,[pathfigures 'AM-test-' num2str(k)],optsfig);
    end

    out = FluctuationStrength_offline_debug(insig,Fs,N, bDebug); %No padding needed for off-line version
    FSAM(:,k) = out{1}; 
    % disp(sprintf('FS=%.3f [vacils]\t test signal: %s\n',out{1},name2figname(filenamesAM{k})))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:length(filenamesFM)
    
    filename = [pathaudio filenamesFM{k} '.wav'];
    
    [x Fs] = Wavread(filename);
    
    afiles{count_afiles} = filenamesAM{k};
    count_afiles = count_afiles + 1;
    
    %%%%
    Fsnew = 44100;
    x = resample(x,Fsnew,Fs);
    Fs = Fsnew;
    %%%%
    
    starti = 1;
    endi    = starti + N-1 + 2048*(N_blocks-1); % 6 additional analysis blocks
    
    insig = x(starti:endi);
    t = ( 1:length(insig) )/Fs;     

    if bDebug % Similar to Fig. 7.2
        h(2) = figure(101);
        subplot(6,1,k)
        plot(t,insig);
        xlabel('Time [s]')
        ylabel('Amplitude')
        title(sprintf('Sig: %s, %.1f [dB]',name2figname( filenamesFM{k}(12:end) ),rmsdb(insig)+90))
        grid on
        ylim(1.5*[minmax(insig')])
        % optsfig.format = 'epsc';
        % Saveas(gcf,[pathfigures 'FM-test-' num2str(k)],optsfig);
    end

    out = FluctuationStrength_offline_debug(insig,Fs,N);%, bDebug); %No padding needed for off-line version
    FSFM(:,k) = out{1}; 
    % disp(sprintf('FS=%.3f [vacils]\t test signal: %s\n',out{1},filenamesFM{k}))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
