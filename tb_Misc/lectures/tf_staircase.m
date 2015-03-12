function tf_staircase
% function tf_staircase
%
% 1. Description:
%       Tested: yes
% 
% Created by Tom Francart
% Edited by Alejandro Osses
% Last use: 10/03/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=[-4;-2;0;2;4;2;0;2;0;2;0;2;0;2;0];

figure;
hold on
plot(s, 'b', 'LineWidth', 2);
xlabel('Trial');
ylabel('SNR (dB)');
title('Adaptive procedure, LIST in noise');
srt=mean(s(end-5:end));
plot(length(s)-[5 0], [1 1] * srt, '--k');
text(length(s)-2, srt, ['SRT = ' num2str(srt) ' dB']);

end