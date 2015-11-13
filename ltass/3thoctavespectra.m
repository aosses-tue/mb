clear;
clc;
load numbers;

F3 = 250*2.^([-1:1/3:6])/2^(1/6);    %Cut-off frequencies for 1/3th octave bands
F3c = 250*2.^([-1:1/3:6]);              %Center frequencies of 1/3th octave bands
F3c = F3c(1:end-1);

for i=1:length(F3c)
    ind = find(F >= F3(i) & F < F3(i+1));
    P3(i) = 10*log10(sum(Pxx(ind)));
    Pn3(i) = 10*log10(sum(Pxxn(ind)));
    PnCD3(i) = 10*log10(sum(PxxnCD(ind)));
    N(i) = length(ind);
end


figure;
semilogx(F3c,P3,F3c,Pn3,F3c,PnCD3);
title('1/3th octave band spectra of numbers');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend({'The 400 numbers';'The calibration noise (numbers\_calibration.wav)';'The calibration noise (CD track)'});

figure;
semilogx(F3c,P3-Pn3);
title('Difference of 1/3th octave band spectra between numbers and numbers-noise');
xlabel('Frequency (Hz)');
ylabel('Power difference (dB)');

clear;
clc;
load wivine;

F3 = 250*2.^([-1:1/3:6])/2^(1/6);    %Cut-off frequencies for 1/3th octave bands
F3c = 250*2.^([-1:1/3:6]);              %Center frequencies of 1/3th octave bands
F3c = F3c(1:end-1);

for i=1:length(F3c)
    ind = find(F >= F3(i) & F < F3(i+1));
    P3(i) = 10*log10(sum(Pxx350(ind)));
    Pn3(i) = 10*log10(sum(Pxxn(ind)));
    PnCD3(i) = 10*log10(sum(PxxnCD(ind)));
end


figure;
semilogx(F3c,P3,F3c,Pn3,F3c,PnCD3);
title('1/3th octave band spectra of sentences');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend({'The 350 sentences';'The calibration noise (astrid\_wom.wav)';'The calibration noise (CD-track)'});

figure;
semilogx(F3c,P3-Pn3);
title('Difference of 1/3th octave band spectra between sentences and sentences-noise');
xlabel('Frequency (Hz)');
ylabel('Power difference (dB)');

AIgain = [0 0.0083 0.0095 0.0150 0.0289 0.0440 0.0578 0.0653 0.0711 0.0818 0.0844 0.0882 0.0898 0.0868 0.0844 0.0771 0.0527 0.0364 0.0185 0 0];