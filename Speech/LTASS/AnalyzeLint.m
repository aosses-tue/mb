function AnalyzeLint

listpath='/mnt/l/speechmaterials/dutch/list/';
lintpath='/mnt/l/speechmaterials/dutch/lint/';
for i=1:10
    lintfiles1_10{i}=[lintpath 'wd' num2str(i) '.wav'];
end

lintfiles_1234568={};
for i=[1:6 8]
    lintfiles_1234568{end+1}=[lintpath 'wd' num2str(i) '.wav'];
end

[Pxxn,F] = AnalyzePSDofWAV([lintpath 'noise.wav'],20);
[Pxxlist,F] = AnalyzePSDofWAV([listpath 'noise/wivineruis.wav'],20);
[Pxx1_10,F] = AnalyzePSDofWAV(lintfiles1_10,20);
[Pxx1234568,F] = AnalyzePSDofWAV(lintfiles_1234568,20);


% figure;
% semilogx(F,20*log10(Pxx),F,20*log10(Pxxn));
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% legend({'LINT wd 1-10', 'LINT full'});

figure;
semilogx(F,smooth(20*log10(Pxx1_10)),F,smooth(20*log10(Pxxn)),F,smooth(20*log10(Pxxlist)), F,smooth(20*log10(Pxx1234568)));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend({'LINT wd 1-10', 'LINT full', 'LIST', 'LINT 1234568'});

save LintSpectra;