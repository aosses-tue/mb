function h = quick_Sentence_adaptive(List, Strat)

if nargin == 0
    List = [-4 -6 -4 -2 -4 -2 0 -2 -4 -6 -8];
    Strat = 0;
end

t = 1:length(List);

nReversals = 6;
SNR4SRT = List(end-nReversals+1:end);
t4SRT = t(end-nReversals+1:end);

figure
plot(t(1:end-nReversals)   , List(1:end-nReversals),'bo'); hold on
plot(t4SRT, SNR4SRT,'ro', 'LineWidth',3)

ylabel('SNR (dB)')
xlabel('Trial Number')
switch Strat
    case 0
        title(['xPC ACE: SRT =' Num2Str(mean(SNR4SRT)) ' dB; SD =' Num2Str(std(SNR4SRT)) ' dB'])
    case 1
        title(['xPC F0mod: SRT =' Num2Str(mean(SNR4SRT)) ' dB; SD =' Num2Str(std(SNR4SRT)) ' dB'])
    case 10
        title(['ACE Own: SRT =' Num2Str(mean(SNR4SRT)) ' dB; SD =' Num2Str(std(SNR4SRT)) ' dB'])
end
ylim(minmax(List)+[-2 2]); grid on

h = gcf;

end