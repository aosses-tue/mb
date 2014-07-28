function y = Set_XTick(ha,lbl_num)
% function y = Set_XTick(ha,lbl_num)
%
% 1. Description:
%   If label has more than 10 elements, decreases the number of XTicks
% 
%       ha - axis handle
%       lbl_num - Numeric labelµ
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
    ha = gca;
end

N = length(lbl_num);
step = 1;

while N/step > 10
    step = step + 1;
end
idx = 1:step:N;

set(ha,'XTick'     ,idx);
set(ha,'XTickLabel',lbl_num(idx));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
