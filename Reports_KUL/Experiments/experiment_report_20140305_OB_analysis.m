function experiment_report_20140305_OB_analysis
% function experiment_report_20140305_OB_analysis
%
% Programmed by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc

filename1 = '/home/alejandro/Documenten/fda_eval_LIST/wav/wivineruis.wav';
filename2 = '/home/alejandro/Documenten/fda_eval/wav/whitenoise-f.wav';
outputname = '/home/alejandro/Documenten/LaTeX_Docs/proefschrift/Figures/OB-analysis-noise';

[f P1 extra1] = One_third_OB_analysis(filename1);
[f P2 extra2] = One_third_OB_analysis(filename2);

semilogx(f,20*log10(P1),'bo-', ... 
         f,20*log10(P2),'ro--');

xlim([100 8000])

grid on, hold on

legend(['SWN, RMS = ' num2str(extra1.rms) ' dB'], ...
       ['White noise, RMS = ' num2str(extra2.rms) ' dB'])
ylabel('Relative amplitude [dB]')
xlabel('Frequency [Hz]')

h = gcf;

Saveas(h, outputname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end