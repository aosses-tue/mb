function r20150617_SpanishMATRIX
% function r20150617_SpanishMATRIX
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 17/06/2015
% Last update on: 17/06/2015 % Update this date manually
% Last use on   : 17/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
bDiary = 0;
Diary(mfilename,bDiary);

bCalculate = 0;

if bCalculate
    freq_one_third = nan(260,28);

    dir_speech = 'D:\Databases\dir03-Speech\spanish\Matrix\';
    type = '';
    % f1 = [dir_speech '00823.wav'];
    fnames = Get_filenames(dir_speech,[type '*.wav']);

    N = 260;
    fn1 = [dir_speech 'VlMatrixnoise_ltass.wav'];
    fn2 = [dir_speech 'EsMatrixnoise_ltass.wav'];

    options.nAnalyser = 10;
    options.CalMethod = 5;
    options.bPlot = 0;

    for i = 1:N
        f1 = [dir_speech fnames{i}];
        out_file1 = PsySoundCL(f1,options);

        freq_one_third(i,:) = out_file1.DataSpecOneThirdAvg;
        disp(['counter: ' num2str(i)]);
    end
    % Saved MANUALLY TILL HERE
else
    load('D:\Output\matlab-SPANISH-MATRIX.mat')
end

out_noise1 = PsySoundCL(fn1,options);
out_noise2 = PsySoundCL(fn2,options);
f = out_file1.f;

SNR2plot = 10;

freq_mean   = mean(freq_one_third);
freq_std    = std(freq_one_third);

figure; % Figure 2.1
semilogx(   f, freq_mean, 'bo--', ...
            f, out_noise1.DataSpecOneThirdAvg-SNR2plot, 'kx--', ...
            f, out_noise2.DataSpecOneThirdAvg-SNR2plot,'r>-');

grid on, hold on
% legend('wdz2, 65 dB SPL','white noise, 55 dB SPL','SSN, 55 dB SPL')
legend('ES-MATRIX','noise1','noise2')
xlabel('Frequency [Hz]')
ylabel('Sound Pressure Level [dB]')
% xlim([50 22050])

%%%

diff_noise1 = freq_mean - (out_noise1.DataSpecOneThirdAvg-SNR2plot);
diff_noise2 = freq_mean - (out_noise2.DataSpecOneThirdAvg-SNR2plot);

figure; % Figure 2.2
semilogx(   f, diff_noise1, 'kx--', ...
            f, diff_noise2, 'r>-');
legend('SNR using white noise','SNR using SSN')
xlabel('Frequency [Hz]')
ylabel('Estimated SNR [dB]')
% ylim([-20 30])
% xlim([50 22050])
grid on, hold on
semilogx([10 20000],[SNR2plot SNR2plot],'k--','LineWidth',2)

[xx fs] = Wavread(fn1);

options.nAnalyser = 1;
out_file1 = PsySoundCL(f1,options);
freq_mean1 = out_file1.Data2;
f1 = out_file1.f;

% CreateEsMatrixManNoise.m

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
