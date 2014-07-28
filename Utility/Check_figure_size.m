function Check_figure_size(h)
% function Check_figure_size(h)
%
% 1. Description:
%       Constrains figure width to A4 Paper size
% 
%   h - Figure handle
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 2/6/2014
% Last update: 2/6/2014 % Update this date manually
% Last used: 2/6/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    h = gcf;
end

Position = get(h,'Position');

Figure_size_cm = px2cm([Position(3)-Position(1)  Position(4)-Position(2)]);
Paper_size_cm = get(h,'PaperSize');
Margins = 0.30; % asuming 30% of margins


if Paper_size_cm(1) < Figure_size_cm(1)
    Rate = (1-Margins)*Paper_size_cm(1)/Figure_size_cm(1);
    Position = floor(Position*Rate);
    set(h,'Position',Position);
    disp([mfilename '.m: Figure re-scaled'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])