function outs = Generate_reference_sounds_Zwicker2007(options)
% function outs = Generate_reference_sounds_Zwicker2007(options)
%
% 1. Description:
%       Generates reference and/or test tones for psychoacoustic metrics
% 
%       1. Reference sound
% 
% 2. Additional info:
%       Tested cross-platform: Yes
%
% 3. Stand-alone example:
%       % To generate roughness files:
%       options.bDoFluct = 0;
%       options.bDoRough = 1;
%       Generate_reference_sounds_Zwicker2007(options);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 14/08/2014
% Last update on: 14/11/2014 % Update this date manually
% Last use on   : 25/11/2014 % Update this date manually
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
pathoutput = [Get_TUe_paths('outputs') 'Fastl2007_test_' p.date4files delim];
Mkdir(pathoutput);

path_src = Get_TUe_data_paths('db_Fastl2007');

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
Fs  = 44100;
dur = options.dur;
sig = Create_sin(fc,dur,Fs,0);

fc125  = 125;
T125   = 1/fc125;
sig125 = Create_sin(fc125,dur,Fs,0);

fc500  = 500;
T500   = 1/fc500;
sig500 = Create_sin(fc500,dur,Fs,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Loudness

if bDoLoud
    lvl = 40;
    lvlAMT  = lvl + 10; % Zwicker's correction
    y = setdbspl(sig,lvlAMT);

    if bDoZeroPadding
        yzero = Zero_padding(y,dur_zero_samples/Fs,Fs);
    end

    if bSave
        filename = [pathoutput 'ref_loud'];
        Wavwrite(yzero,Fs,filename);
        outs.filename{end+1} = filename;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Fluctuation strength
if bDoFluct
    
    %% 3.0 Reference
    fmod    = 4;
    m       = 100;
    option  = 'm';
    lvlref  = 60;
    lvl     = 70;
    
    start_phase = -pi/2; % pi/2;
    y = ch_am(sig,fmod,m,option,Fs,start_phase);
    
    lvlAMT  = lvlref + 10; % Zwicker's correction
    
    if bPsySound
        y = setdbspl(y,lvlAMT);
    else
        y = setdbspl(y,lvlref);
    end

    if bDoRamp
        ramp2apply = cos_ramp(length(y),Fs,dur_ramp_ms);
        y       = ramp2apply'.*y;
    end
    
    if bDoZeroPadding
        y   = Zero_padding(y,dur_zero_samples/Fs,Fs);
    end
    
    if bSave
        filename = [pathoutput 'ref_fluct'];
        Wavwrite(y,Fs,filename);
        outs.filename{end+1} = filename;
    end
    
    %% 3.2 Fastl2007, Fig. 10.2 (from track_04)
    % Modulated tones:
    filename = [path_src 'track_04_t08_white.wav'];
    [sig Fs] = Wavread(filename);
    m       = [0 6 10 20 40 50 60 70 80 90 95 98 100];
    option   = 'm';
    start_phase = -pi/2;
    lvl = 60;
    
    fmod = 4; 
    
    fc      = 1000;
    dur     = length(sig)/Fs;
    sigTone = Create_sin(fc,dur,Fs,0);
    
    for k = 1:length(m)
        
        lvlAMT  = lvl + 10;
        y  = ch_am(sig    ,fmod,m(k),option,Fs,start_phase);
        yT = ch_am(sigTone,fmod,m(k),option,Fs,start_phase);
        
        if bPsySound
            y = setdbspl(y,lvlAMT);
            yT = setdbspl(yT,lvlAMT);
        else
            y = setdbspl(y,lvl);
            yT = setdbspl(yT,lvl);
        end

        if bDoRamp
            ramp2apply = cos_ramp(length(y),Fs,dur_ramp_ms);
            y    = ramp2apply'.*y;
            ramp2apply = cos_ramp(length(yT),Fs,dur_ramp_ms);
            yT    = ramp2apply'.*yT;
        end

        if bDoZeroPadding
            y   = Zero_padding(y,dur_zero_samples/Fs,Fs);
            yT   = Zero_padding(yT,dur_zero_samples/Fs,Fs);
        end

        filename = [pathoutput 'fluct_test_bbn_AM_m_' Num2str(m(k),3) '_fmod_' Num2str(floor(fmod),3) 'Hz_' num2str(lvl) '_dBSPL'];
        Wavwrite(y   ,Fs,filename);
        outs.filename{end+1} = filename;
        
        filename = [pathoutput 'fluct_test_fc_' Num2str(fc,3) '_AM_m_' Num2str(m(k),3) '_fmod_' Num2str(floor(fmod),3) 'Hz_' num2str(lvl) '_dBSPL'];
        Wavwrite(yT  ,Fs,filename);
        outs.filename{end+1} = filename;

    end

    %% 3.3 Fastl2007, Fig. 10.3
    % Stimuli can be reconstructed Fig.10.2-stimuli
    
    %% 3.4 Fastl2007, Fig. 10.4
    fc      = [125 250 500  1000 2000 4000 8000];
    d       = 40; % dB
    option  = 'd';
    lvl     = 70;
    fmod    = 4;
    
    %   Fig. 10.4.a: AM tones
    for k = 1:length(fc) 
        
        lvlAMT = lvl + 10;
        sig = Create_sin(fc(k),dur,Fs,0);
        
        y  = ch_am(sig,fmod,d,option,Fs,start_phase);
        
        if bPsySound
            y = setdbspl(y,lvlAMT);
        else
            y = setdbspl(y,lvl);
        end

        if bDoRamp
            ramp2apply = cos_ramp(length(y),Fs,dur_ramp_ms);
            y    = ramp2apply'.*y;
        end

        if bDoZeroPadding
            y   = Zero_padding(y,dur_zero_samples/Fs,Fs);
        end

        filename = [pathoutput 'fluct_test_fc_' Num2str(fc(k),4) '_AM_d_' Num2str(d,2) '_fmod_' Num2str(floor(fmod),3) 'Hz_' num2str(lvl) '_dBSPL'];
        Wavwrite(y  ,Fs,filename);
        outs.filename{end+1} = filename;
        
    end
    
    %   Fig. 10.4.b: FM tones
    deltaf  = 200;
    lvl     = 70;
    fc      = [500  1000 2000 4000 8000];
    
    for k = 1:length(fc)
        
        lvlAMT  = lvl + 10;
        yfm     = fm(fc(k), dur, Fs, fmod, deltaf);

        if bDoRamp
            ramp2apply = cos_ramp(length(yfm),Fs,dur_ramp_ms);
            yfm     = ramp2apply'.*yfm;
        end

        if bPsySound
            yfm = setdbspl(yfm,lvlAMT);
        else
            yfm = setdbspl(yfm,lvl);
        end
        filename = [pathoutput 'fluct_test_fc_' Num2str(fc(k),4) '_FM_dev_' Num2str(deltaf,3) '_fmod_' Num2str(floor(fmod),3) 'Hz_' num2str(lvl) '_dBSPL'];
        Wavwrite(yfm   ,Fs,filename);
        outs.filename{end+1} = filename;
        
    end
    
    %% 3.5 Fastl2007, Fig. 10.5
    fc      = 1500;
    deltaf  = [0 16 30 60 120 150 400 600 800 1000];
    lvl     = 70;
    fmod    = 4;
    
    for k = 1:length(deltaf)
        
        lvlAMT  = lvl + 10;
        yfm     = fm(fc, dur, Fs, fmod, deltaf(k));

        if bDoRamp
            ramp2apply = cos_ramp(length(yfm),Fs,dur_ramp_ms);
            yfm     = ramp2apply'.*yfm;
        end

        if bPsySound
            yfm = setdbspl(yfm,lvlAMT);
        else
            yfm = setdbspl(yfm,lvl);
        end
        
        if bDoZeroPadding
            y = Zero_padding(y,dur_zero_samples/Fs,Fs);
        end
        
        filename = [pathoutput 'fluct_test_fc_' Num2str(fc   ,4) '_FM_dev_' Num2str(deltaf(k),3) '_fmod_' Num2str(floor(fmod),3) 'Hz_' num2str(lvl) '_dBSPL'];
        Wavwrite(yfm   ,Fs,filename);
        outs.filename{end+1} = filename;
    end
    
    %% 3.6 Fastl2007, Fig. 10.6 (from track_36)
    %                           BW       
    %       'track_36_t01.wav'  2      
    %       'track_36_t02.wav'  6
    %       'track_36_t03.wav'  50
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Roughness

if bDoRough
    
    %% 4.0 Reference
    fmod    = 70;
    m       = 100;
    option = 'm';
    lvl     = 60;
    
    start_phase = pi/2;
    y    = ch_am(sig   ,fmod,m,option,Fs,start_phase);
    
    lvlAMT  = lvl + 10; % Zwicker's correction
    if bPsySound
        y    = setdbspl(y   ,lvlAMT);
    else
        y    = setdbspl(y   ,lvl);
    end
    
    if bDoZeroPadding
        y = Zero_padding(y,dur_zero_samples/Fs,Fs);
    end
    
    if bDoRamp
        ramp2apply = cos_ramp(length(y),Fs,dur_ramp_ms);
        y = ramp2apply'.*y;
        disp('Ramp applied')
    end
    
    if bSave
        filename = [pathoutput 'rough_ref'];
        Wavwrite(y,Fs,filename);
        outs.filename{end+1} = filename;
    end
    
    %% 4.1 Fastl2007, Fig. 11.1 (from track_38)
    %                           m       rms [dB]    dB SPL
    %       'track_38_t01.wav'  1       -23.24      
    %       'track_38_t02.wav'  0.7
    %       'track_38_t03.wav'  0.4
    %       'track_38_t04.wav'  0.25
    %       'track_38_t05.wav'  0.125
    %       'track_38_t06.wav'  0.1
    %       'track_38_t07.wav'  0
    
    %% 4.2 Fastl2007, Fig. 11.2
    fc_test = [125 250 500 1000 2000 4000 8000];
    fmod_test = [20 30 50 70 100 150];
    m       = 100;
    option  = 'm';
    start_phase = pi/2;
    lvl = 60;
    
    for k = 1:length(fc_test)
        
        fc      = fc_test(k);
        sigt    = Create_sin(fc,dur,Fs,0);
        
        for j = 1:length(fmod_test)
            fmod    = fmod_test(j);
            
            y    = ch_am(sigt  ,fmod,m,option,Fs,start_phase);

            lvlAMT  = lvl + 10; % Zwicker's correction
            if bPsySound
                y    = setdbspl(y   ,lvlAMT);
            else
                y    = setdbspl(y   ,lvl);
            end

            if bDoZeroPadding
                y = Zero_padding(y,dur_zero_samples/Fs,Fs);
            end

            if bDoRamp
                ramp2apply = cos_ramp(length(y),Fs,dur_ramp_ms);
                y = ramp2apply'.*y;
                disp('Ramp applied')
            end

            if bSave
                filename = [pathoutput 'rough_test_fc_' Num2str(fc   ,4) '_AM_m_' Num2str(m,3) '_fmod_' Num2str(fmod,3) 'Hz'];
                Wavwrite(y,Fs,filename);
                outs.filename{end+1} = filename;
            end
        end
        
    end
    
    %% 4.3 Fastl2007, Fig. 11.3
    %                               fmod
    % track_39_t01.wav  AM BBN      20 Hz, d = 40 dB 
    % track_39_t02.wav              70
    % track_39_t03.wav              200
    % track_39_t04.wav  AM Tone     20
    % track_39_t05.wav              70  
    % track_39_t06.wav              200
    % track_39_t07.wav  FM Tone     20
    % track_39_t08.wav              70
    % track_39_t09.wav              200
    
    %% 4.4 Fastl2007, Fig. 11.4
    % Adjust levels from 40:10:80
    
    %% 4.5 Fastl2007, Fig. 11.5 - FM tones
    fc      = 1500;
    deltaf  = [30 60 120 150 400 600 800 1000];
    lvl     = 60;
    fmod    = 70;
    
    for k = 1:length(deltaf)
        
        lvlAMT  = lvl + 10;
        yfm     = fm(fc, dur, Fs, fmod, deltaf(k));

        if bDoRamp
            ramp2apply = cos_ramp(length(yfm),Fs,dur_ramp_ms);
            yfm     = ramp2apply'.*yfm;
        end

        if bPsySound
            yfm = setdbspl(yfm,lvlAMT);
        else
            yfm = setdbspl(yfm,lvl);
        end
        
        if bDoZeroPadding
            y = Zero_padding(y,dur_zero_samples/Fs,Fs);
        end
        
        filename = [pathoutput 'rough_test_fc_' Num2str(fc   ,4) '_FM_dev_' Num2str(deltaf(k),3) '_fmod_' Num2str(floor(fmod),3) 'Hz_' num2str(lvl) '_dBSPL'];
        Wavwrite(yfm   ,Fs,filename);
        outs.filename{end+1} = filename;
    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
