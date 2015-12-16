function CreateEsMatrixManNoise()
% function CreateEsMatrixManNoise()
%
% Created by TF
% Adapted/Edited by Alejandro Osses V.
% Last used on 18/06/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_matrix = 'D:\Databases\dir03-Speech\spanish\Matrix\';
dir_analysis = [dir_matrix '02-Analysis' delim];
try
    load([dir_analysis 'EsMatrixSpectra.mat']);
catch
    [filename fileplusdir] = Get_filenames(dir_matrix,['*.wav']);
    filename = fileplusdir(1:260);
    [Pxx_average, F ,Pxx] = AnalysePSDofWAV(filename);
    save([dir_analysis 'EsMatrixSpectra'],'F','Pxx','Pxx_average');
end

Pxx = Pxx_average;
% 2^12 = 4096
% 2^11 = 2048
% 2^10 = 1024
%                 (order,norm-f ,amplitudes, grid to interpolate)
[H1,Fn]=freqz( fir2(2^12,F/22050,sqrt(Pxx)), 1);
[H2,Fn]=freqz( fir2(2^11,F/22050,sqrt(Pxx)), 1);
[H3,Fn]=freqz( fir2(2^10,F/22050,sqrt(Pxx)), 1);

semilogx(Fn, 10*log10(abs([H1 H2 H3])));

legend({'4096 taps' '2048 taps' '1024 taps'});

% return;

%B = fir2(2^11,F/22050,sqrt(mean([Pxxnogaps2'; Pxxnogaps3'])'));
B = fir2(2^11,F/22050,sqrt(Pxx));       %% sqrt want Pxx is vermogen
w = randn(11*44100,1);

y = filter(B,1,w);
filename = [dir_analysis 'aonoise_AO20150618.wav'];
Wavwrite(y,44100,16,filename);
Pxxn = AnalysePSDofWAV(filename,20,0);

lvl = rmsdb(y);
w2save = setdbspl(w,lvl+100);
filename = [dir_analysis 'genwhitenoise_AO20150618.wav'];
Wavwrite(w2save,44100,16,filename);
Pxxnt = AnalysePSDofWAV(filename,20,0);

figure;semilogx(F,10*log10(Pxxnt),F,10*log10(Pxxn),F,10*log10(Pxx));
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');

title('Spectra of the EsMatrix Lists and of the ''new'' noise');
legend({'White noise';'Speech shaped noise (aonoise.wav)';'Speech without gaps (all lists)'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end