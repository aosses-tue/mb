function outs = r20141107_roughness
% function outs = r20141107_roughness
%
% 1. Description:
%       Implement the model of Roughness (off-line)
% 
% 2. Stand-alone example:
%       r20141107_roughness;
% 
%       outs = r20141107_roughness;
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%       See also r20150706_comparing_roughness_FS.m, r20141126_roughness_validation.m
% 
% Programmed by Alejandro Osses V., HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 05/11/2014
% Last update on: 14/11/2014 % Update this date manually
% Last use on   : 05/07/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathaudio   = Get_TUe_paths('outputs');
pathfigures = [Get_TUe_paths('outputs') 'Figures-20141113' delim];
Mkdir(pathfigures);

filename    = [pathaudio 'ref_rough'];

N           = 8192;
bDebug      = 1;

filenamesAM = { ... %'Audacity-1kHz-amp-08',...
                'test_rough_fc_1000_AM_m_100_fmod_030Hz_60_dBSPL', ...
                'test_rough_fc_1000_AM_m_100_fmod_050Hz_60_dBSPL', ...
                'test_rough_fc_1000_AM_m_100_fmod_070Hz_60_dBSPL', ...
                'test_rough_fc_1000_AM_m_100_fmod_090Hz_60_dBSPL', ...
                'test_rough_fc_1000_AM_m_100_fmod_110Hz_60_dBSPL', ...
                'test_rough_fc_1000_AM_m_100_fmod_130Hz_60_dBSPL', ...
                'test_rough_fc_1000_AM_m_100_fmod_150Hz_60_dBSPL'};

filenamesFM = { 'test_rough_fc_1000_FM_dev_800_fmod_020Hz_60_dBSPL',...
                'test_rough_fc_1000_FM_dev_800_fmod_040Hz_60_dBSPL', ...
                'test_rough_fc_1000_FM_dev_800_fmod_060Hz_60_dBSPL', ...
                'test_rough_fc_1000_FM_dev_800_fmod_070Hz_60_dBSPL', ...
                'test_rough_fc_1000_FM_dev_800_fmod_090Hz_60_dBSPL', ...
                'test_rough_fc_1000_FM_dev_800_fmod_110Hz_60_dBSPL', ...
                'test_rough_fc_1000_FM_dev_800_fmod_150Hz_60_dBSPL', ...
                'test_rough_fc_1000_FM_dev_800_fmod_200Hz_60_dBSPL'};


if nargout == 1
    outs.filenames = [filenamesAM'; filenamesFM'];
    outs.dir_where = pathaudio;
    return;
end

try
    
    Wavread(filename);
    
catch % if reference does not exist, then it is created...
    
    opts.bDoRough = 1;
    opts.dur = (N/44100); %  approx. 200e-3;
    opts.bGen_test_tones = 1;
    
    opts.bPsySound = 1;
    opts.bDoRamp = 0; 
    % opts.dur_ramp_ms = 25;
    opts.bDoZeroPadding = 0;
    % opts.dur_zero_samples = N/2; % 4096 + 4096
    
    outs = Generate_reference_sounds(opts);
end

for k = 1:length(filenamesAM)
    
    close all
    
    filename = [pathaudio filenamesAM{k} '.wav'];
    
    [x Fs] = Wavread(filename);
    
    %%%%
    Fsnew = 44100;
    x = resample(x,Fsnew,Fs);
    Fs = Fsnew;
    %%%%
    
    starti = 1;
    insig = x(starti:starti + N-1);
    t = ( 1:length(insig) )/Fs;     

    if bDebug % Similar to Fig. 7.2
        figure(1)
        plot(t,insig);
        xlabel('Time [s]')
        ylabel('Amplitude')
        title(sprintf('Test signal %s.wav, level = %.2f [dB]',name2figname(filenamesAM{k}),rmsdb(insig)+90))
        grid on
        
        optsfig.format = 'epsc';
        Saveas(gcf,[pathfigures 'AM-test-' num2str(k) '-1-time-series'],optsfig);
    end

    out = Roughness_offline_debug(insig,Fs,N, bDebug); %No padding needed for off-line version
    % out = FluctuationStrength_offline_debug(insig,Fs,N, bDebug); %No padding needed for off-line version
    RAM(k)  = out{1};
    h       = out{4};
    
    for i = 1:length(h)
        try 
            Maximise_figure(h(i));
        catch
            warning('No maximisation of the figure');
        end
        
        optsfig.format = 'epsc';
        optsfig.bScale = 0;
        Saveas(h(i),[pathfigures 'AM-test-' num2str(k) '-' num2str(i+1)],optsfig);
        
    end
    disp(sprintf('R=%.3f [aspers]\t test signal: %s\n',out{1},name2figname(filenamesAM{k})))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:length(filenamesFM)
    
    close all
    filename = [pathaudio filenamesFM{k} '.wav'];
    
    [x Fs] = Wavread(filename);
    
    %%%%
    Fsnew = 44100;
    x = resample(x,Fsnew,Fs);
    Fs = Fsnew;
    %%%%
    
    starti = 1;
    insig = x(starti:starti + N-1);
    t = ( 1:length(insig) )/Fs;     

    if bDebug % Similar to Fig. 7.2
        figure(1)
        plot(t,insig);
        xlabel('Time [s]')
        ylabel('Amplitude')
        title(sprintf('Test signal %s.wav, level = %.2f [dB]',name2figname(filenamesFM{k}),rmsdb(insig)+90))
        grid on
        
        optsfig.format = 'epsc';
        optsfig.bScale = 0;
        Saveas(gcf,[pathfigures 'FM-test-' num2str(k)],optsfig);
    end

    out = Roughness_offline_debug(insig,Fs,N, bDebug); %No padding needed for off-line version
    % out = FluctuationStrength_offline_debug(insig,Fs,N, bDebug); %No padding needed for off-line version
    RFM(k) = out{1}; 
    h       = out{4};
    
    for i = 1:length(h)
        try 
            Maximise_figure(h(i));
        catch
            warning('No maximisation of the figure');
        end
        
        optsfig.format = 'epsc';
        Saveas(h(i),[pathfigures 'FM-test-' num2str(k) '-' num2str(i+1)],optsfig);
        
    end
    disp(sprintf('R=%.3f [aspers]\t test signal: %s\n',out{1},filenamesFM{k}))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
