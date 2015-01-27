function y = r20150123_update_VoD
% function y = r20150123_update_VoD
%
% 1. Description:
%       Data presented in report generated in GUI PsySoundControl.fig. Use
%       the screenshots saved inside the folder of each condition to see
%       the details of each configuration.
% 
%       In parallel, I started to generate the test battery for Fluctuation
%       Strength (nothing written about this by my side)
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 23/01/2015
% Last update on: 23/01/2015 % Update this date manually
% Last use on   : 23/01/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);
close all

bPart1 = 0;
bPartFluct = 0;
bPartFluct2= 1;
bPartRough = 0;

if bPart1
% directory = '/home/alejandro/Documenten/Databases/dir01-Instruments/Mech-Hummer/';
directory = '~/Documenten/Databases/dir01-Instruments/new-ac-mode2-close+far/';

f1 = [directory 'model-ac-mode-2.wav'];
f2 = [directory 'meas-ac-mode-2.wav']; % file already aligned
f3 = [directory 'model-ac-mode-2-far.wav'];
f4 = [directory 'meas-ac-mode-2-far.wav']; 

acmode  = 5;

% model
[x fs]  = Wavread(f1);
x2      = Wavread(f2);
x3      = Wavread(f3);
x4      = Wavread(f4);

t       = ( 1:length(x) )/fs;
misc    = Get_VoD_params(0);
T       = misc.Tmodel(acmode-1); % misc.ti_model(acmode-1);

% 8 ac mode 2
% 3 ac mode 5
ti = 3*T; % starting from 8th period
tf = ti + 4*T; % ending 4 periods later

tis = ti*fs
tfs = tf*fs

sprintf('ti = %.0f; tf = %.0f',tis,tfs) % Put this in PsySoundControl.fig

% measured
offset = 0.04; % just to plot

figure;
plot(t,x +3*offset), hold on, grid on
plot(t,x2+1*offset,'r')

plot(t,x3-1*offset)
plot(t,x4-3*offset,'r')

xlim([ti tf])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cal file: white noise:

dir_where_ref = Get_TUe_paths('db_fastl2007');
dir_output  = [Get_TUe_paths('outputs') 'tones' delim];

file    = [dir_where_ref 'track_03.wav'];

if bPartFluct
    
    % For Figure 10.a (Fastl 2007)
    %   - AM-white noise
    
    t       = 2;
    
    [insig fs] = Wavread(file);
    options.callevel = 60;
    
    fmod    = [1 2 4 8 16 32];
    m       = 100;
    option  = 'm';
    % lvlref  = 60;
    % lvl     = 60;
    
    start_phase = -pi/2; % pi/2;
    
    for i = 1:length(fmod)
        [y env] = ch_am(insig,fmod(i),m,option,fs,start_phase);
        % Wavwrite( y(1:t*fs),fs, sprintf('%sbbn_AM_m_%s_fmod_%sHz',dir_output,Num2str(m,3),Num2str(fmod(i),3)) );
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   - AM-sine-tone
    fmod    = 4;
    m       = [0 10 50 70 100];
    option  = 'm';
    % lvlref  = 60;
    % lvl     = 60;
    
    start_phase = -pi/2; % pi/2;
    
    fc  = 1000;
    T   = 1/fc;
    fs  = 44100;
    t   = 2;
    sig = Create_sin(fc,t,fs,0);
    
    for i = 1:length(m)
        
        [y env] = ch_am(insig,fmod,m(i),option,fs,start_phase);
        y = setdbspl(y,70);
        Wavwrite( y(1:t*fs),fs, sprintf('%sbbn_AM_m_%s_fmod_%sHz_%.0fdB',dir_output,Num2str(m(i),3),Num2str(fmod,3),rmsdb(y)+90) );
        
        [yt envt] = ch_am(sig,fmod,m(i),option,fs,start_phase);
        yt = setdbspl(yt,70);
        Wavwrite( yt(1:t*fs),fs, sprintf('%stone_%s_AM_m_%s_fmod_%sHz_%.0fdB',dir_output,Num2str(fc,4), Num2str(m(i),3),Num2str(fmod,3),rmsdb(yt)+90) );
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   - FM-sine-tone
    fs = 44100;
    deltaf = [100 400 700]; % Hz
    
    for i = 1:length(deltaf)
        yf = fm(fc,t,fs,fmod,deltaf(i));
        yf = setdbspl(yf,70);
        Wavwrite( yf(1:t*fs),fs, sprintf('%stone_%s_FM_dev_%s_fmod_%sHz_%.0fdB',dir_output,Num2str(fc,4), Num2str(deltaf(i),3),Num2str(fmod,3),rmsdb(yf)+90) );
    end
end

close all
if bPartFluct2
   
    fl = {  'tone_1000_AM_m_100_fmod_004Hz_60dB.wav', ...
            'tone_1000_AM_m_050_fmod_004Hz_60dB.wav', ...
            'tone_1000_AM_m_010_fmod_004Hz_60dB.wav'};
    
    res1 = [];
    res2 = [];
    k = 17;
    for i = 1:length(fl)
        [x fs] = Wavread([dir_output fl{i}]);
        t = ( 1:length(x) )/fs;
        figure;
        plot(t,x);
        
        out = fs_offline(x,fs,0);
        res1 = [res1 transpose( out{1}(k,:) )];
        res2 = [res2 transpose( out{2}(k,:) )];
    end
    
    fl = {  'tone_1000_FM_dev_700_fmod_004Hz_60dB.wav', ...
            'tone_1000_FM_dev_400_fmod_004Hz_60dB.wav', ...
            'tone_1000_FM_dev_100_fmod_004Hz_60dB.wav'};
    
    res1f = [];
    res2f = [];
    k = 17;
    for i = 1:length(fl)
        [x fs] = Wavread([dir_output fl{i}]);
        t = ( 1:length(x) )/fs;
                
        out = fs_offline(x,fs,0);
        res1f = [res1f transpose( out{1}(k,:) )];
        res2f = [res2f transpose( out{2}(k,:) )];
    end

    fl = {  'bbn_AM_m_100_fmod_004Hz_60dB.wav', ...
            'bbn_AM_m_050_fmod_004Hz_60dB.wav', ...
            'bbn_AM_m_010_fmod_004Hz_60dB.wav'};
    
    res1b = [];
    res2b = [];
    for i = 1:length(fl)
        [x fs] = Wavread([dir_output fl{i}]);
        t = ( 1:length(x) )/fs;
        % figure;
        % plot(t,x);
        % title('BBN')
        
        out = fs_offline(x,fs,0);
        res1b = [res1b transpose(out{1}(k,:))];
        res2b = [res2b transpose(out{2}(k,:))];
    end
    
end

%--
figure; % AM
subplot(1,2,1)
plot(res1) % eim1
title('AM')
hold on
plot(res1) % eim1

subplot(1,2,2)
plot(res2) % eim2
legend('1','2','3')

%--
figure; 
subplot(1,2,1)
plot(res1b)
title('BBN')

subplot(1,2,2)
plot(res2b)
legend('1','2','3')

%--
figure; 
subplot(1,2,1)
plot(res1f)
title('FM')

subplot(1,2,2)
plot(res2f)
legend('1','2','3')

%--

if bPartRough
    
    % For Figure 11.4.a (Fastl2007)
    t       = 1;
    
    [insig fs] = Wavread(file);
    % options.callevel = 60;
    
    fmod    = [30 50 70 90 110];
    m       = 100;
    option  = 'm';
    lvlref  = 60;
    lvl     = 60;
    
    start_phase = -pi/2; % pi/2;
    
    for i = 1:length(fmod)
        [y env] = ch_am(insig,fmod(i),m,option,fs,start_phase);
        Wavwrite( y(1:t*fs),fs, sprintf('%sbbn_AM_m_%s_fmod_%sHz',dir_output,Num2str(m,3),Num2str(fmod(i),3)) );
    end
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
