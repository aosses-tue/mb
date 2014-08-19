function Generate_reference_sounds
% function Generate_reference_sounds
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
% Last update on: 18/08/2014 % Update this date manually
% Last use on   : 18/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bZeroPad = 1;
bSave = 1;

bDoLoud = 0;
bDoSharp = 0;
bDoFluct = 0;
bDoRough = 1;

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
if bDoFluct
    fmod    = 4;
    m       = 100;
    option = 'm';
    lvl     = 60;

    start_phase = pi/2;
    y = ch_am(sig,fmod,m,option,fs,start_phase);

    lvlAMT  = lvl + 10; % Zwicker's correction
    y = setdbspl(y,lvlAMT);

    if bZeroPad
        yzero = Zero_padding(y,200e-3,fs);
    end

    if bSave
        Wavwrite(yzero,fs,[Get_TUe_paths('outputs') 'ref_fluct']);
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
