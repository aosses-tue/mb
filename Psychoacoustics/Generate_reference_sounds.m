function outs = Generate_reference_sounds(options)
% function outs = Generate_reference_sounds(options)
%
% 1. Description:
%       Generates reference and/or test tones for psychoacoustic metrics
% 
% 2. Additional info:
%       Tested cross-platform: Yes
%
% 3. Stand-alone example:
%       % To generate roughness files:
%       options.bDoFluct = 0;
%       options.bDoRough = 1;
%       Generate_reference_sounds(options);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 14/08/2014
% Last update on: 26/11/2014 % Update this date manually
% Last use on   : 26/11/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    options = [];
end

bSave   = 1;

options = Ensure_field(options,'bDoLoud',0);
options = Ensure_field(options,'bDoSharp',0);
options = Ensure_field(options,'bDoFluct',0);
options = Ensure_field(options,'bDoRough',1);
options = Ensure_field(options,'bPsySound',0); % important if fluctuation strength is to be determined
options = Ensure_field(options,'dur',4); % Test tone duration, 4 seconds by default
options = Ensure_field(options,'bDoRamp',0);

if options.bDoRamp
    options = Ensure_field(options,'dur_ramp_ms',75);
end

options = Ensure_field(options,'bDoZeroPadding',0); % previous default value: 200e-3 s
if options.bDoZeroPadding
    options = Ensure_field(options,'dur_zero_samples',0);
end

p = Get_date;
pathoutput = [Get_TUe_paths('outputs') 'Daniel1997_test_' p.date4files delim];
Mkdir(pathoutput);

options = Ensure_field(options,'bGen_test_tones',1);

bDoLoud     = options.bDoLoud;
bDoSharp    = options.bDoSharp;
bDoFluct    = options.bDoFluct;
bDoRough    = options.bDoRough;

bPsySound       = options.bPsySound;
bGen_test_tones = options.bGen_test_tones;
bDoRamp         = options.bDoRamp;
bDoZeroPadding  = options.bDoZeroPadding;

if bDoRamp;         dur_ramp_ms = options.dur_ramp_ms;              end;
if bDoZeroPadding;  dur_zero_samples = options.dur_zero_samples;    end;

outs.filename = {};
% Common params
fc  = 1000;
T   = 1/fc;
fs  = 44100;
dur = options.dur;
sig = Create_sin(fc,dur,fs,0);

fc125  = 125;
T125   = 1/fc125;
sig125 = Create_sin(fc125,dur,fs,0);

fc500  = 500;
T500   = 1/fc500;
sig500 = Create_sin(fc500,dur,fs,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Loudness

if bDoLoud
    lvl = 40;
    lvlAMT  = lvl + 10; % Zwicker's correction
    y = setdbspl(sig,lvlAMT);

    if bDoZeroPadding
        yzero = Zero_padding(y,dur_zero_samples/fs,fs);
    end

    if bSave
        filename = [Get_TUe_paths('outputs') 'ref_loud'];
        Wavwrite(yzero,fs,filename);
        outs.filename{end+1} = filename;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Sharpness % 
if bDoSharp
    fmod    = 70;
    m       = 100;
    option  = 'm';
    lvl     = 70;

    if isunix
        sigbp = Wavread('~/Documenten/MATLAB/outputs/track_04_t09_bp.wav'); % temporal location on 17/08/2014
    end

    start_phase = pi/2; % begin in maximum. Use -pi/2 to begin in minimum
    % ch_am(sig,fmod,m,option,fs,start_phase); % uncomment this to plot
    y = ch_am(sigbp,fmod,m,option,fs,start_phase);

    lvlAMT  = lvl + 10; % Zwicker's correction
    y = setdbspl(y,lvlAMT);

    if bDoZeroPadding
        y = Zero_padding(y,dur_zero_samples/fs,fs);
    end

    if bSave
        filename = [Get_TUe_paths('outputs') 'ref_sharp'];
        Wavwrite(y,fs,filename);
        outs.filename{end+1} = filename;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Fluctuation strength
% 3.1 AM tones
if bDoFluct
    
    fmod    = 4;
    m       = 100;
    option  = 'm';
    lvlref  = 60;
    lvl     = 60;
    
    start_phase = -pi/2; % pi/2;
    y = ch_am(sig,fmod,m,option,fs,start_phase);
    y2 = ch_am(sig,fmod,0,option,fs,start_phase);
    
    lvlAMT  = lvlref + 10; % Zwicker's correction
    
    if bPsySound
        y = setdbspl(y,lvlAMT);
    else
        y = setdbspl(y,lvlref);
    end

    if bDoRamp
        ramp2apply = cos_ramp(length(y),fs,dur_ramp_ms);
        y       = ramp2apply'.*y;
    end
    
    if bDoZeroPadding
        y   = Zero_padding(y,dur_zero_samples/fs,fs);
    end
    
    if bSave
        filename = [Get_TUe_paths('outputs') 'ref_fluct'];
        Wavwrite(y,fs,filename);
        outs.filename{end+1} = filename;
    end
    
    if bGen_test_tones
        % Modulated tones:
        fi = 0.5; % 0,5 Hz to start
        for k = 0:6
            lvlAMT  = lvl + 10;
            fmod = fi*2^k;
            d       = 40;
            option = 'd';
            y = ch_am(sig,fmod,m,option,fs,start_phase);

            if bPsySound
                y = setdbspl(y,lvlAMT);
            else
                y = setdbspl(y,lvl);
            end
           
            if bDoRamp
                ramp2apply = cos_ramp(length(y),fs,dur_ramp_ms);
                y    = ramp2apply'.*y;
            end
            
            if bDoZeroPadding
                y   = Zero_padding(y,dur_zero_samples/fs,fs);
            end
            
            filename = [Get_TUe_paths('outputs') 'test_fluct_fc_' Num2str(fc   ,4) '_AM_m_' Num2str(m,3) '_fmod_' Num2str(floor(fmod),3) 'Hz_' num2str(lvl) '_dBSPL'];
            Wavwrite(y   ,fs,filename);
            outs.filename{end+1} = filename;

        end
        
        % 3.2 FM tones
        fc      = 1000;
        deltaf  = 700;
        % tzero   = 0; % time to Zero-pad
        lvl     = 60;

        %% Modulated tones:
        fi = 0.5; % 0,5 Hz to start
        for k = 0:6
            lvlAMT  = lvl + 10;
            fmod    = fi*2^k;
            yfm     = fm(fc, dur, fs, fmod, deltaf);
            
            if bDoRamp
                ramp2apply = cos_ramp(length(yfm),fs,dur_ramp_ms);
                yfm     = ramp2apply'.*yfm;
            end

            if bPsySound
                yfm = setdbspl(yfm,lvlAMT);
            else
                yfm = setdbspl(yfm,lvl);
            end
            filename = [Get_TUe_paths('outputs') 'test_fluct_fc_' Num2str(fc   ,4) '_FM_dev_' Num2str(deltaf,3) '_fmod_' Num2str(floor(fmod),3) 'Hz_' num2str(lvl) '_dBSPL'];
            Wavwrite(yfm   ,fs,filename);
            outs.filename{end+1} = filename;
        end
    else
        disp('Test tones for ''Fluctuation Strength'', FM, not generated')
    end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Roughness

if bDoRough
    fmod    = 70;
    m       = 100;
    option = 'm';
    lvl     = 60;
    
    start_phase = pi/2;
    y    = ch_am(sig   ,fmod,m,option,fs,start_phase);
    
    lvlAMT  = lvl + 10; % Zwicker's correction
    if bPsySound
        y    = setdbspl(y   ,lvlAMT);
    else
        y    = setdbspl(y   ,lvl);
    end
    
    if bDoZeroPadding
        y = Zero_padding(y,dur_zero_samples/fs,fs);
    end
    
    if bDoRamp
        ramp2apply = cos_ramp(length(y),fs,dur_ramp_ms);
        y = ramp2apply'.*y;
        disp('Ramp applied')
    end
    
    % if bSave
    %     filename = [Get_TUe_paths('outputs') 'ref_rough'];
    %     Wavwrite(y,fs,filename);
    %     outs.filename{end+1} = filename;
    % end
    
    % fmod_test = [30 50 70 100 150];
    % for k = 1:length(fmod_test)
    %     y    = ch_am(sig   ,fmod_test(k),m,option,fs,start_phase);
    % 
    %     lvlAMT  = lvl + 10; % Zwicker's correction
    % 
    %     if bPsySound
    %         y    = setdbspl(y   ,lvlAMT);
    %     else
    %         y    = setdbspl(y   ,lvl);
    %     end
    % 
    %     if bDoRamp
    %         ramp2apply = cos_ramp(length(y),fs,dur_ramp_ms);
    %         y = ramp2apply'.*y;
    %         disp('Ramp applied')
    %     end
    % 
    %     if bDoZeroPadding
    %         y = Zero_padding(y,dur_zero_samples/fs,fs);
    %     end
    % 
    %     if bSave
    %         filename = [Get_TUe_paths('outputs') 'test_rough_fc_' Num2str(fc   ,4) '_AM_m_' Num2str(m,3) '_fmod_' Num2str(fmod_test(k),3) 'Hz'];
    %         Wavwrite(y   ,fs, filename);
    %         outs.filename{end+1} = filename;
    %     end
    % 
    % end
    % 
    % % Modulated tones:
    % for k = 10:20:170
    %     fmod = k;
    %     y       = ch_am(sig,fmod,m,option,fs,start_phase);
    %     y125    = ch_am(sig125,fmod,m,option,fs,start_phase);
    %     y500    = ch_am(sig500,fmod,m,option,fs,start_phase);
    % 
    %     lvlAMT  = lvl + 10; % Zwicker's correction
    %     if bPsySound
    %         y    = setdbspl(y   ,lvlAMT);
    %         y125 = setdbspl(y125,lvlAMT);
    %         y500 = setdbspl(y500,lvlAMT);
    %     else
    %         y    = setdbspl(y   ,lvl);
    %         y125 = setdbspl(y125,lvl);
    %         y500 = setdbspl(y500,lvl);
    %     end
    % 
    %     if bDoRamp
    %         y    = ramp2apply'.*y;
    %         y125 = ramp2apply'.*y125;
    %         y500 = ramp2apply'.*y500;
    %     end
    % 
    %     filename = [Get_TUe_paths('outputs') 'test_rough_fc_' Num2str(fc   ,4) '_AM_m_' Num2str(m,3) '_fmod_' Num2str(fmod,3) 'Hz_' num2str(lvl) '_dBSPL'];
    %     Wavwrite(y   ,fs,filename);
    %     outs.filename{end+1} = filename;
    % 
    %     filename = [Get_TUe_paths('outputs') 'test_rough_fc_' Num2str(fc125,4) '_AM_m_' Num2str(m,3) '_fmod_' Num2str(fmod,3) 'Hz_' num2str(lvl) '_dBSPL'];
    %     Wavwrite(y125,fs,filename);
    %     outs.filename{end+1} = filename;
    % 
    %     filename = [Get_TUe_paths('outputs') 'test_rough_fc_' Num2str(fc500,4) '_AM_m_' Num2str(m,3) '_fmod_' Num2str(fmod,3) 'Hz_' num2str(lvl) '_dBSPL'];
    %     Wavwrite(y500,fs,filename);
    %     outs.filename{end+1} = filename;
    % 
    % end

    %% 4.3 FM tones
    % Daniel1997, Fig.9
    fc      = 1600;
    fmod    = [0 10 20 40 50 60 70 80 100 120 150 200 500];
    deltaf  = 800;
    lvl     = 60;

    for k = 1:length(fmod)

        lvlAMT  = lvl + 10;
        yfm     = fm(fc, dur, fs, fmod(k), deltaf);

        if bDoRamp
            ramp2apply = cos_ramp(length(yfm),fs,dur_ramp_ms);
            yfm     = ramp2apply'.*yfm;
        end

        if bPsySound
            yfm = setdbspl(yfm,lvlAMT);
        else
            yfm = setdbspl(yfm,lvl);
        end

        if bDoZeroPadding
            y = Zero_padding(y,dur_zero_samples/fs,fs);
        end

        filename = [pathoutput 'rough_test_fc_' Num2str(fc   ,4) '_FM_dev_' Num2str(deltaf,3) '_fmod_' Num2str(floor(fmod(k)),3) 'Hz_' num2str(lvl) '_dBSPL'];
        Wavwrite(yfm   ,fs,filename);
        outs.filename{end+1} = filename;
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
