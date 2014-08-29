function Generate_reference_sounds(options)
% function Generate_reference_sounds(options)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: Yes
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 14/08/2014
% Last update on: 22/08/2014 % Update this date manually
% Last use on   : 22/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    options = [];
end

bZeroPad = 1;
bSave = 1;

bDoLoud = 0;
bDoSharp = 0;
options = Ensure_field(options,'bDoFluct',1);
bDoRough = 0;

% Common params
fc  = 1000;
T   = 1/fc;
fs  = 44100;
dur = 4;
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

    if bZeroPad
        yzero = Zero_padding(y,200e-3,fs);
    end

    if bSave
        Wavwrite(yzero,fs,[Get_TUe_paths('outputs') 'ref_loud']);
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

    if bZeroPad
        yzero = Zero_padding(y,200e-3,fs);
    end

    if bSave
        Wavwrite(yzero,fs,[Get_TUe_paths('outputs') 'ref_sharp']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Fluctuation strength
% 3.1 AM tones
if options.bDoFluct
    tzero   = 0; % time to Zero-pad
    fmod    = 4;
    m       = 100;
    option  = 'm';
    lvlref  = 60;
    lvl     = 70;
    bForPsySound = 0;
    
    start_phase = pi/2;
    y = ch_am(sig,fmod,m,option,fs,start_phase);

    lvlAMT  = lvlref + 10; % Zwicker's correction
    
    if bForPsySound
        y = setdbspl(y,lvlAMT);
    else
        y = setdbspl(y,lvlref);
    end

    ramp2apply = cos_ramp(length(y),fs,75);
    y       = ramp2apply'.*y;
    
    if bZeroPad
        y   = Zero_padding(y,tzero,fs);
    end
    
    if bSave
        Wavwrite(y,fs,[Get_TUe_paths('outputs') 'ref_fluct']);
    end
    
    % Modulated tones:
    fi = 0.5; % 0,5 Hz to start
    for k = 0:6
        lvlAMT  = lvl + 10;
        fmod = fi*2^k;
        d       = 40;
        option = 'd';
        y = ch_am(sig,fmod,m,option,fs,start_phase);
        % y125 = ch_am(sig125,fmod,m,option,fs,start_phase);
        % y500 = ch_am(sig500,fmod,m,option,fs,start_phase);
        
        if bForPsySound
            y = setdbspl(y,lvlAMT);
        else
            y = setdbspl(y,lvl);
        end
        %y125 = setdbspl(y125,lvlAMT);
        %y500 = setdbspl(y500,lvlAMT);
        
        y    = ramp2apply'.*y;
        %y125 = ramp2apply'.*y125;
        %y500 = ramp2apply'.*y500;
        
        Wavwrite(y   ,fs,[Get_TUe_paths('outputs') 'test_fluct_fc_' Num2str(fc   ,4) '_AM_m_' Num2str(m,3) '_fmod_' Num2str(floor(fmod),3) 'Hz']);
        %Wavwrite(y125,fs,[Get_TUe_paths('outputs') 'test_fluct_fc_' Num2str(fc125,4) '_AM_m_' Num2str(m,3) '_fmod_' Num2str(fmod,3) 'Hz']);
        %Wavwrite(y500,fs,[Get_TUe_paths('outputs') 'test_fluct_fc_' Num2str(fc500,4) '_AM_m_' Num2str(m,3) '_fmod_' Num2str(fmod,3) 'Hz']);
        
    end
    
end
% 3.2 FM tones
if options.bDoFluct
    fc      = 1500;
    deltaf  = 700;
    % tzero   = 0; % time to Zero-pad
    lvl     = 70;
    bForPsySound = 0;
    
    % Modulated tones:
    fi = 0.5; % 0,5 Hz to start
    for k = 0:6
        lvlAMT  = lvl + 10;
        fmod    = fi*2^k;
        yfm     = fm(fc, dur, fs, fmod, deltaf);
        
        ramp2apply = cos_ramp(length(yfm),fs,75);
        yfm     = ramp2apply'.*yfm;
        
        if bForPsySound
            yfm = setdbspl(yfm,lvlAMT);
        else
            yfm = setdbspl(yfm,lvl);
        end
        
        Wavwrite(yfm   ,fs,[Get_TUe_paths('outputs') 'test_fluct_fc_' Num2str(fc   ,4) '_FM_dev_' Num2str(deltaf,3) '_fmod_' Num2str(floor(fmod),3) 'Hz']);
        
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
    y    = setdbspl(y   ,lvlAMT);
    
    ramp2apply = cos_ramp(length(y),fs,75);
    yzero    = ramp2apply'.*y;
    
    if bSave
        Wavwrite([yzero] ,fs,[Get_TUe_paths('outputs') 'ref_rough']);
    end
    
    % Modulated tones:
    for k = 10:20:170
        fmod = k;
        y = ch_am(sig,fmod,m,option,fs,start_phase);
        y125 = ch_am(sig125,fmod,m,option,fs,start_phase);
        y500 = ch_am(sig500,fmod,m,option,fs,start_phase);
        
        lvlAMT  = lvl + 10; % Zwicker's correction
        y = setdbspl(y,lvlAMT);
        y125 = setdbspl(y125,lvlAMT);
        y500 = setdbspl(y500,lvlAMT);
        
        y    = ramp2apply'.*y;
        y125 = ramp2apply'.*y125;
        y500 = ramp2apply'.*y500;
        
        Wavwrite(y   ,fs,[Get_TUe_paths('outputs') 'test_rough_fc_' Num2str(fc   ,4) '_AM_m_' Num2str(m,3) '_fmod_' Num2str(fmod,3) 'Hz']);
        Wavwrite(y125,fs,[Get_TUe_paths('outputs') 'test_rough_fc_' Num2str(fc125,4) '_AM_m_' Num2str(m,3) '_fmod_' Num2str(fmod,3) 'Hz']);
        Wavwrite(y500,fs,[Get_TUe_paths('outputs') 'test_rough_fc_' Num2str(fc500,4) '_AM_m_' Num2str(m,3) '_fmod_' Num2str(fmod,3) 'Hz']);
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
