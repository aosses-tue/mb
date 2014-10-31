function r20141031_noise
% function r20141031_noise
%
% 1. Description:
%
% 2. Stand-alone example:
%       r20141031_noise;
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 28/10/2014
% Last update on: 28/10/2014 % Update this date manually
% Last use on   : 28/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

%%
dira = Get_TUe_data_paths('db_audacity'); % dir audacity files
diro = Get_TUe_paths('outputs');

[xA fs ] = Wavread([dira 'white-5-sec-onset-1-sec-20-5000-Hz.wav']);
t = ( 1:length(xA) )/fs;

%%

nTag = 99; % 5-s white noise, 1-s onset; 1-s offset
filename = 'test';
options.dB_SPL_noise = 77;
options.fs = fs;
options.bSave_noise = 0;
innoise = Create_noise_dau1996(nTag,filename,options);

%%
ti = 1;
tf = 6;

rmsa = rmsdb(xA,fs,ti,tf)+100;
diff2apply = options.dB_SPL_noise - rmsa;

xA = From_dB(diff2apply)*xA;
rmsa = rmsdb(xA,fs,ti,tf)+100;

disp( num2str(rmsa) )

rmsn = rmsdb(innoise,fs,ti,tf)+100;
disp( num2str(rmsn) )

figure;
plot( t,20*log10(abs(xA)),t,20*log10(abs(innoise)) ), legend('audacity','MATLAB')
ylim([-20 0])
xlim([0 tf+1])

hn = freqz(innoise,1,4096);

ha = freqz(xA,1,4096);

K = 4096;
figure;
freqfft([hn ha],K,options);

%%
[xA1 fsA1] = Wavread([diro 'dau1996b_expI1_noisemasker.wav']); % nTag = 1
[xA3 fsA3] = Wavread([diro 'dau1996b_expI3_noisemasker.wav']); % nTag = 3
[xB0 fsB0] = Wavread([diro 'dau1996b_expIB0_noisemasker.wav']); % nTag = 20
tB= ( 1:length(xB0) )/fsB0;

nTag = 1;
options.fs = fsB0;
innoiseA1 = Create_noise_dau1996(nTag,filename,options);

nTag = 3;
innoiseA3 = Create_noise_dau1996(nTag,filename,options);

nTag = 20;
innoiseB0 = Create_noise_dau1996(nTag,filename,options);


figure;
subplot(3,1,1)
plot(   tB, xA1, ...
        tB, innoiseA1)
% legend('experiment','just generated')
ha = gca;

subplot(3,1,2)
plot(   tB, xA3, ...
        tB, innoiseA3)
% legend('experiment','just generated')
ha(end+1) = gca;

subplot(3,1,3)
plot(   tB, xB0, ...
        tB, innoiseB0)
% legend('experiment','just generated')
ha(end+1) = gca;

linkaxes(ha,'xy')

xlim([-0.2 0.7])

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
