function my_answers
% function my_answers

close all

tmax = 5; % just going to process 5 seconds
nSents = 5; % use a number between 1 and 10, % number of concatenated sentences 

bPart1 = 0;
bPart2 = 1;

fnoise1 = 'SSN_HINT_22kHz.wav';
fnoise2 = 'SSN_HINT_modReduced_-20dB.wav';
testSNRs = -4:4:4; %-6:2:6;

if bPart1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

[x1 fs1] = Wavread(fnoise1);
[x2 fs2] = Wavread(fnoise2);

x1 = x1(1:fs1*tmax);
x2 = x2(1:fs2*tmax);

%% Answer 2: 'Plot their long-term spectra'
AudioBands  = [ 125  250  500  1000 2000 4000 8000 ];
% AudioBands  = [25 32 40 50 63 80 100 125 160 200 250 315 400 500 630 800 1000  2000 4000 8000 1250 1600 2000 2500 3150 4000 5000 6300 8000];

ymean1 = Get_values_OB(x1,fs1,AudioBands);
ymean2 = Get_values_OB(x2,fs2,AudioBands);

figure;
plot( AudioBands,[ymean1; ymean2])
legend(name2figname(fnoise1),name2figname(fnoise2))
grid on

%% Answer 3: 'Would you expect the AI/SII to predict a difference between these two noise tyoes?'
% No
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if bPart2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Answer 4: How does the STI change as a function of the input SNR?
load Danish_HINT_10sentence_samples_22kHz % sentences stored as sentenceArray{i} with i from 1 to 10
fs = 22050;

insig = [];
for q = 1:nSents
    insig = [insig sentenceArray{q}'];
end

[STI1 MTF_mean1 MTF_std1 out1] = Get_STI(insig,fnoise1,fs,testSNRs);
[STI2 MTF_mean2 MTF_std2]      = Get_STI(insig,fnoise2,fs,testSNRs);
modBands = out1.modBands;

figure
plot(testSNRs,[STI1 STI2],'s', 'markersize',10)
xlabel('Input SNR (dB)')
ylabel('STI')
ylim([0 1])
xlim([-7 7])
legend('SSN', 'Modulation filtered SSN');
grid on
% STI increases as the noise decreases (i.e, when SNR increase)

%% Answer 5: Plot MTF for the different input SNRs

figure
for i = 1:7
    subplot(3,3,i)
    errorbar(testSNRs-0.1,[MTF_mean1(:,i)],[MTF_std1(:,i)]), hold on
    errorbar(testSNRs+0.1,[MTF_mean2(:,i)],[MTF_std2(:,i)],'r')
    title(sprintf('Mod-band %.2f [Hz]',modBands(i)))
    grid on
end
legend('SSN', 'Mod-filtered SSN','Location','EastOutside');

figure;
for i = 8:14
    subplot(3,3,i-7)
    errorbar(testSNRs-0.1,[MTF_mean1(:,i)],[MTF_std1(:,i)]), hold on
    errorbar(testSNRs+0.1,[MTF_mean2(:,i)],[MTF_std2(:,i)],'r')
    title(sprintf('Mod-band %.2f [Hz]',modBands(i)))
    grid on
end
legend('SSN', 'Mod-filtered SSN','Location','EastOutside');

% xlabel('Input SNR (dB)')
% ylabel('STI')
% ylim([0 1])
% xlim([-7 7])

% grid on

disp('')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function ymean = Get_values_OB(x,fs,AudioBands)

[out_time out_spec] = OctFilterBank(x,fs,AudioBands);

% sound(out_time(:,1),fs) % to listen to  band 1
% sound(out_time(:,2),fs) % to listen to  band 2

N = size(out_spec,1);
sp = 20*log10( abs(out_spec(1:N/2,:)) +eps);
[xx ymean] = sum_db(sp,20*log10(eps));
