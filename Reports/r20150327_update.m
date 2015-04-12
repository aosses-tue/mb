function r20150327_update
% function r20150327_update
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 25/03/2015
% Last update on: 25/03/2015 % Update this date manually
% Last use on   : 25/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

bDiary = 1;
Diary(mfilename,bDiary);

bPart1 = 0; % Generation of audio files at the proper sampling rate
bPart2 = 0; % Dau1997
bPart3 = 1; % Stimulus generation to determine variance of internal noise

path = Get_TUe_paths('db_instruments');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart1
    
    f1 = [path 'Piano' delim '00-Original-files' delim 'pressionexpeCd5.wav']; % at 48 kHz
    f2 = [path 'Piano' delim '00-Original-files' delim 'pressionsimuCd5.wav']; % at 50 kHz
    
    [x1 fs]     = Wavread(f1); 
    [x2 fs2]    = Wavread(f2);
    
    fsnew = 48000;
    
    fo1 = [path 'Piano' delim 'pressionexpeCd5.wav']; 
    fo2 = [path 'Piano' delim 'pressionsimuCd5.wav']; 
    
    if fs ~= fsnew
        y1 = resample(x1,fsnew,fs);
    else
        y1 = x1;
    end
    
    if fs2 ~= fsnew
        y2 = resample(x2,fsnew,fs2);
    else
        y2 = x2;
    end
    
    Wavwrite(y1,fsnew,fo1);
    Wavwrite(y2,fsnew,fo2);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart2
    
    outputLevel = [];
    pref    = 2e-5;
    dur     = 1000e-3;
    Li      = 0;
    Lf      = 100;
    fs      = 44100;
    t       = 1/fs:1/fs:dur;
    tmin    = 800e-3;
    tmax    = 1000e-3;
    Lp      = Li:10:Lf;
    
    cl = {'b','r'};
    count = 1;
    figure;
    for Lpidx = Lp 
        A = 0.5*pref*10^( Lpidx/20 ); % A signal producing 100 dB
        insig = A*ones(fs*dur,1);
        
        [outsig outputLevel(count)] = Il_get_adaptloop( insig,fs );
        
        if mod(count,2)==0
            plot(t,outsig, cl{1}), hold on
        else
            plot(t,outsig, cl{2}), hold on
        end
        
        count = count + 1;
    end
    
    str1 = [num2str(Li) ' : 20' ];
    str2 = [num2str(Li+10) ' : 20' ];
    
    legend(str1,str2)
    xlim([tmin tmax])
    xlabel('Time [s]')
    ylabel('Output [MU]')
    grid on
    
    % Dau1996a, Fig. 3
    figure;
    plot(Lp,outputLevel,'o-'), grid on, hold on
    xlabel('Input level [dB]')
    ylabel('Output level [MU]')
    
    plot(Lp,Lp,'r--');
    
    jnd61 = interp1(Lp,outputLevel,[60 61]);
    jnd61 = diff(jnd61);
    
    count = 1;
    for Lpidx = Lp +1
        A = 0.5*pref*10^( Lpidx/20 ); % A signal producing 100 dB
        insig = A*ones(fs*dur,1);
        
        [outsig outputLevelp1(count)] = Il_get_adaptloop( insig,fs );
        
        count = count + 1;
    end
    diffOut = outputLevelp1 - outputLevel;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f       = 1000;
    [y, t]  = Create_sin(f,dur,fs);
    A       = setdbspl(60+3);
    ycal    = A*y;
    
    out = [];
    out1 = [];
    for lvl = 0:-20:-40
        
        var = From_dB(lvl);
        [ytmp, outLvl, idxs]    = Il_get_int_repr( ycal           ,fs,var );
        
        [ytmp1, outLvl1, idxs1] = Il_get_int_repr( ycal*From_dB(1),fs,var );
        out     = [out  ytmp(:,idxs(1))];
        out1    = [out1 ytmp1(:,idxs1(1))];
        
        sprintf('JND = %.3f [dB], int.noise var of %.3f [dB]',outLvl1-outLvl,lvl)
    end
    
    figure;
    plot(t,out), grid on
    xlim([tmin tmax])
    
    % legend('-100','-95','-90','-85','-80')
    % legend('-0','-50','-100')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bPart3
    
    bListen = 1;
    bSaveWave = 1;
    
    %% Dau1997b, Fig. 7:
    % Suprathreshold signal...

    % Common parameters:
    fs = 44100;
    BW = 3; %[3 31 314];
    % d = [-40];
    mdept = 0;%d2m(d,'dau');
    
    % For BBN:
    finf = 20;
    fsup = 10000;
    
    fc_tone = 5000;
    
    SPL = 65;
    SPLtest = SPL;
    dur = 3;
    fm  = 20;
    
    insig_BBN   = AM_random_noise(finf,fsup,SPL,dur,fs);
    insig_NBN   = AM_random_noise_BW(fc_tone,BW,SPL,dur,fs,fm,mdept);
    insig_test  = AM_sine(fc_tone,dur,fs,fm,mdept,SPLtest);

    % insig_NBN = [insig_NBN];
    % insig_test = [insig_test];

    t = ( 0:length(insig_NBN)-1 )/fs;

    if bListen == 1

        sound(insig_BBN,fs);
        
        sound(insig_NBN,fs);

        sound(insig_test,fs);
    end
    
    if bSaveWave
        
        fname = sprintf('%sBBN-BW-%.0f-SPL-%.0fdB',Get_TUe_paths('outputs'),fsup-finf,SPL);
        Wavwrite(insig_BBN,fs,fname);
        
        fname = sprintf('%stest-fc-%.0f-SPL-%.0fdB',Get_TUe_paths('outputs'),fc_tone,SPLtest);
        Wavwrite(insig_test,fs,fname);
        
        att = -20;
        fname = sprintf('%sBBN-BW-%.0f-SPL-%.0fdB',Get_TUe_paths('outputs'),fsup-finf,SPL+att);
        Wavwrite(From_dB(att)*insig_BBN,fs,fname)
        
        att = -40;
        fname = sprintf('%sBBN-BW-%.0f-SPL-%.0fdB',Get_TUe_paths('outputs'),fsup-finf,SPL+att);
        Wavwrite(From_dB(att)*insig_BBN,fs,fname)
        
    end

end
    
disp('')

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inline signals:
function [outsig outLevel] = Il_get_adaptloop(insig,fs)

limit   = 10;
minlvl  = 1e-5;
outsig  = adaptloop(insig,fs,limit,minlvl);

N       = size(insig,1);
len     = round(N/10);
tmp     = outsig([end-len+1:end],:);
outLevel = mean(tmp); % last element, should be most stable value

function [ytmp,outLevel,idxs] = Il_get_int_repr( ycal,fs,var )
  
model   = 'dau1996';
[ytmp xx xx opts] = Get_internal_representations(ycal,fs,model,var);
idxs = opts.idx;

N       = size(ycal,1);
len     = round(N/10);
tmp     = ytmp([end-len+1:end],idxs(1));
outLevel = mean(tmp);

