function Plot_specific_loudness(Ns)
% function Plot_specific_loudness(Ns)
%
% 1. Description:
%   Plots frame by frame the Ns function. Each frame is separated in Dt time
% 
%   Ns - specific loudness
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%       res = Zwicker_dynamic_loudness_model(ymeas, fs, info);
%       Plot_specific_loudness( res.InstantaneousSpecificLoudness )
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 3/7/2014
% Last update on: 3/7/2014 % Update this date manually
% Last used on: 3/7/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a,b] = size(Ns);
max_global = max(max(Ns));

for i = 1:a
    plot(Ns(i,:));
    ylim([0 max_global])
    pause(0.1)
    fprintf('%.0f frames out of %.0f frames completed\n',i,a);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end