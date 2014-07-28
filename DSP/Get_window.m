function [win,wtype] = Get_window(nwtype,N,M)
% function [win,wtype] = Get_window(nwtype,N,M)
%
% 1. Description:
%       Generates a window in a column vector with length N. Default value
%       of M is 1. If you have a signal stored in a column vector use M to 
%       get the window immediately repeated in columns.
%
%       nwtype - window type (number arbitrarily assigned)
%       wtype  - window type (window name)
%       win    - array containing window
% 
%       nwtype      wtype                   Tested  Tested on
%       0           'rectangular'           No
%       1           'hanning'               Yes     24/06/2014
%       2           'triangular'            No
%       3           'bartlett'              No
%       4           'hamming'               No
%       5           'blackman'              No
%       6           'blackman-harris'       No
%       7           'gaussian(alpha=2.5)'   No
%       8           'gaussian(alpha=3.5)'   No
%       9           'gaussian(alpha=4.5)'   No
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 16/6/2014
% Last update: 24/6/2014 % Update this date manually
% Last used: 26/6/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    M = 1;
end

if length(N) ~= 1 % then N is assumed to be a 'y'-vector
    N = length(N); 
end

switch nwtype
    case 0
        win = ones(N,1); % 'rectangular window'
        wtype = 'rectangular';
    case 1
        win = hanning(N,'periodic'); % periodic for analysis purposes, symmetric for filtering purposes
        wtype = 'hanning';
    case 2 
        win = triang(N);
        wtype = 'triangular';
    case 3
        win = bartlett(N);
        wtype = 'bartlett';
    case 4
        win = hamming(N,'periodic');
        wtype = 'hamming';
    case 5
        win = blackman(N);
        wtype = 'blackman';
    case 6
        H     = sigwin.blackmanharris(N);
        win   = generate(H);
        wtype = 'blackman-harris';
    case 7
        win   = gausswin(N,2.5);
        wtype = 'gaussian(alpha=2.5)';
    case 8
        win   = gausswin(N,3.5);
        wtype = 'gaussian(alpha=3.5)';
    case 9
        win   = gausswin(N,4.5);
        wtype = 'gaussian(alpha=4.5)';
    otherwise
        error('Not recognised window')
end

win = repmat(win,1,M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end