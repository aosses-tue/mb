function str = Num2str(num, numDecimals) %, bEqualLenght)
% function str = Num2str(num, numDecimals) %, bEqualLenght)
%
% 1. Description:
%       num           - number to be rounded
%       numDecimals   - number of decimals
%
%       if num is an integer number then numDecimals represents the number 
%       of characters of num when converted to string 
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       Num2str(5.2345634)  % answer = 5.23
%       Num2str(5,3)        % answer = 005
%
% Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium, 2012-2013
% Created in    : 2012-2013
% Last update on: 19/08/2014 % Update this date manually
% Last use on   : 19/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1
    numDecimals = 2;
end

if nargin < 3
    % bEqualLength = 1;
    Length = 6;
end

str = [];


if mod(num,1) ~= 0 % then it is a decimal...
    K = min( size(num,1),size(num,2) );

    for k = 1:K
        for i = 1:length(num)
            str_tmp = num2str( round(num(i)*10^(numDecimals))/10^(numDecimals) );

            if length(str_tmp) < Length
                Remaining = Length - length(str_tmp);

                if length(num) > 1
                    str_tmp = [' ' str_tmp];
                end

                for j = 1:Remaining
                    str_tmp = [' ' str_tmp];
                end
            end
            str = [str str_tmp];
        end
    end
else % then it is an integer
    str = num2str(num);
    
    tmp = [];
    while length(tmp)+length(str) < numDecimals
        tmp = ['0' tmp];
    end
    str = [tmp str];
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
