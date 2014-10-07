function stLB = getLBValues(test_regs)
% function stLB = getLBValues(test_regs)
%
% Script for getting manual Loudness Balance values, this is done by means
% of a linear interpolation
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bPlot = 1;
if nargin == 0
    test_regs       = [104, 147, 208, 294];
end

Freqs = [];

if ~isfield(Freqs, 'Fmin')
    Freqs.Fmin      = min(test_regs);
end

if ~isfield(Freqs, 'Fmax')
    stFmax          = freq2note(max(test_regs)); % 4 semitones above Fmax
    Fmax            = sumSemitones2note( stFmax, 4); % 4 semitones above Fmax
    Freqs.Fmax      = round(note2freq(Fmax));
end

k = 1;
Freqs.F(k) = Freqs.Fmin;
stFmin = freq2note(Freqs.Fmin);
newFreq = Freqs.Fmin; % Initial value

while newFreq < Freqs.Fmax
    [xx xx newFreq] = sumSemitones2note(stFmin,k);
    k = k + 1;
    Freqs.F(k) = newFreq;
end

aceData = input([mfilename '.m - Give ACE   values for frequencies ' num2str(test_regs) ' [Hz]: ']);
f0mData = input([mfilename '.m - Give F0mod values for frequencies ' num2str(test_regs) ' [Hz]: ']);

Freq_ref = 131;
count_less_than_ref = find(test_regs < Freq_ref);
test_regs_ACE   = [test_regs(1:count_less_than_ref)  Freq_ref   test_regs(count_less_than_ref+1:end)];
aceData_ACE     = [aceData(1,1:count_less_than_ref)  0          aceData(1,count_less_than_ref+1:end)];
display(['Inside ' mfilename ' Caution, reference tone of 131 Hz is being used for interpolation'])

LB_ACE(1,:) = interp1(test_regs_ACE, aceData_ACE ,Freqs.F,'linear','extrap');
LB_F0m(1,:) = interp1(test_regs    , f0mData(1,:),Freqs.F,'linear','extrap');

if bPlot == 1
    figure
    plot(   test_regs, aceData(1,:), 'LineWidth', 2, 'Color', 'r', 'Marker', 'o', 'LineStyle','o'), hold on
    plot(   test_regs, f0mData(1,:), 'LineWidth', 2, 'Color', 'b', 'Marker', 'x', 'LineStyle','x')

    plot(   Freqs.F, LB_ACE(1,:), 'r-')
    plot(   Freqs.F, LB_F0m(1,:), 'b-')

    legend('ACE', 'F0m'), grid on
    title(['Loudness Balance for Subject, generated inside: ' mfilename '.m'])
    ylabel('\Delta dB, ref 131 Hz @ 60 dB(A)')
    xlabel('Frequency [Hz]')
    ylim([-25 25])
end

stLB.F = Freqs.F;
stLB.LB_ACE = LB_ACE;
stLB.LB_F0m = LB_F0m;

end