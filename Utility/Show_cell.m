function Show_cell(in)
% function Show_cell(in)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 12/04/2016
% Last update on: 12/04/2016 
% Last use on   : 12/04/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N, M] = size(in);

for j = 1:M
    for i = 1:N
        if isstr( in{i,j} )
            fprintf('(%.0f,%.0f): %s\n', i,j,in{i,j});
        elseif isnumeric( in{i,j} )
            fprintf('(%.0f,%.0f): %.4f\n', i,j,in{i,j});
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
