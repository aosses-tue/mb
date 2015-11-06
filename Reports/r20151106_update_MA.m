function r20151106_update_MA
% function r20151106_update_MA
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 04/11/2015
% Last update on: 04/11/2015 
% Last use on   : 04/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

dir_files = ['D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2015-10-30-update-to-music-instruments\Audio' delim];
f1 = 'meas-ac-4-dist-ane-HP.wav';

bSave = 0;
bPart1 = 0;
bPart2 = 1;
bPart3 = 0;

if bPart1
try
    load([dir_files 'Spectra-hummer.mat']);
catch
    [insig fs] = Wavread([dir_files f1]);
    N = round(500e-3 * fs);
    K = N/2;
    insigbuf = buffer(insig,N,K,'nodelay');

    [Pin F] = Get_PSD_analysis_arrays(insigbuf,fs);

    save([dir_files 'Spectra-hummer'],'F','Pin','fs');
end

Pxx = Pin;
% 2^12 = 4096
% 2^11 = 2048
% 2^10 = 1024
%                 (order,norm-f ,amplitudes, grid to interpolate)
[H1,Fn]=freqz( fir2(2^12,F/22050,sqrt(Pxx)), 1);
[H2,Fn]=freqz( fir2(2^11,F/22050,sqrt(Pxx)), 1);
[H3,Fn]=freqz( fir2(2^10,F/22050,sqrt(Pxx)), 1);

figure;
semilogx(Fn, 10*log10(abs([H1 H2 H3])));

legend({'4096 taps' '2048 taps' '1024 taps'});

B = fir2(2^11,F/22050,sqrt(Pxx));       %% sqrt because Pxx is energy
w = randn(11*44100,1);

y = filter(B,1,w);
filename = [dir_files 'hummer_noise_20151104.wav'];
if bSave == 1
    Wavwrite(y,44100,16,filename);
end
Pxxn = Get_PSD_analysis_arrays(y,fs,20,0);

lvl = rmsdb(y);
w2save = setdbspl(w,lvl+100);

filename = [dir_files 'whitenoise_20151104.wav'];
if bSave == 1
    Wavwrite(w2save,44100,16,filename);
end
Pxxnt = Get_PSD_analysis_arrays(w2save,fs,20,0);

figure;
semilogx(F,10*log10(Pxxnt),F,10*log10(Pxxn),F,10*log10(Pxx));
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');

title('Spectra of the Hummer sounds and the ''new'' noise');
legend({'White noise';'Hummer-shaped noise (hummer_noise.wav)';'Average hummer spectrum (all segments)'});
end

if bPart2
    
    dir_files = ['D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Hummer\demo-HWN+meas-sim-hummer\Stimuli\' delim];
    files = {   'hummer_noise_20151104.wav', ...
                'meas-ac-4-dist-ane-HP+LP-1-seg.wav', ...
                'model-ac-4-dist-ane-1-seg.wav'};
    for i = 1:length(files)
        [x fs] = Wavread([dir_files files{i}]);
        meanvals(i) = rmsdb(x)+100;
    end
    
end

if bPart3
    mu = 0;
    sigma = [0.6 1 2 5 10 20 50 100];
    
    load([Get_TUe_paths('outputs') 'template-20151104.mat']);
    load([Get_TUe_paths('outputs') 'calc-20151104.mat']);
    
    idx1 = 8500;
    idx2 = 28500;
    sigint1   = sigint1(idx1:idx2,:);
    sigint1s2 = sigint1s2(idx1:idx2,:);
    sigint2   = sigint2(idx1:idx2,:);
    template = template(idx1:idx2,:);
    intM = intM(idx1:idx2,:);
    xx = xx(idx1:idx2,:);
    yy = yy(idx1:idx2,:);
    zz = zz(idx1:idx2,:);
    
    N = size(template,1);
    Ntimes = 100;
    
    % intM
    intN1 = sigint1 - xx;
    intN2 = sigint1s2 - yy;
    intSN = sigint2 - zz;
    
    for k = 1:length(sigma)
        
        for i = 1:Ntimes

            noise = normrnd(mu,sigma(k),N,6);
            intN1t = intN1 + noise(:,1); 
            intN2t = intN2 + noise(:,2);
            intSNt = intSN + noise(:,3);
            intMt1 = intM + noise(:,4);
            intMt2 = intM + noise(:,5);
            intMt3 = intM + noise(:,6);

            decision(i,1) = optimaldetector( intN1t-intMt1 ,template,fs_intrep);
            decision(i,2) = optimaldetector( intN2t-intMt2,template,fs_intrep);
            % decision(i,3) = optimaldetector( intSNt-intMt3,template,fs_intrep); 
        end
        sig(k,:) = prctile(decision(:),[75 95]); %[prctile(decision(:),25) prctile(decision(:),50) prctile(decision(:),75)];
        
    end
        
    disp('')
end
if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
