% out = Adaptlevel(in, level, switches);
%
% Scales signal ’in’ to the demanded digital ’level’ in dB SPL
% with respect to the selected module implementations defined
% in array ’switches’.
%
function out = AdaptLevel(in, level, switches)
% Calculate current rms value

[rows cols] = size(in);
rmsSig = sqrt(sum(in.^2, 2) / cols);
% Calculate the factor to use for adapting the signal,
% depending on the adaptation model
if switches(3) == 2
    factor = 10.^(level / 20 - log10(rmsSig));
else
    factor = 10.^((level - 30) / 20 - log10(rmsSig));
end

% Adapt the signal
out = zeros(rows, cols);
for i = 1:rows,
    out(i,:) = in(i,:) * factor(i);
end
