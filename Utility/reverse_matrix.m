function y = reverse_matrix(x)
% function y = reverse_matrix(x)
%
%                           Tested
% x - row vector (1 x N)    OK
% x - column vector (M x 1) NOT TESTED
% x - Matrix (M x N)        NOT TESTED
% 
% Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium, 2012-2013
% Created on  : 2012-2013
% Last updated: 19/6/2014 % Update this date manually
% Last used   : 19/6/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a, b] = size(x);

Minimum = min(a,b);
Maximum = max(a,b);

if Minimum == a
    
    if Minimum == 1 % row vector
        for i=1:Maximum
            y(i) = x(end-i+1);
        end 
    else % Matrix
        for i=1:a
            y(a+1-i,:) = x(i,:);
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end