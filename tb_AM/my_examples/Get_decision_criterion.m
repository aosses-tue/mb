function out = Get_decision_criterion(ir1,ir2,param,opts)
% function out = Get_decision_criterion(ir1,ir2,param,opts)
%
% 1. Description:
%
%       'cross-correlation', using optimal detector:
%           The cross correlation between the expected signal sm and the 
%           received signal xm is monotonic with likelihood ratio. The optimal
%           detection scheme for the signal specified exactly is therefore
%           the cross correlation between the expected signal and the received
%           waveform. If the receiver computes the cross correlation and
%           accepts the hypothesis that the signal was present only when the
%           cross correlation is above some criterion value, it will be the
%           optimal detector (Green & Swets 1966, pp. 163.)
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 30/04/2015
% Last update on: 03/05/2015 % Update this date manually
% Last use on   : 03/05/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    opts = [];
end

opts = Ensure_field(opts,'idx',1:size(ir1,2));
opts = Ensure_field(opts,'fs',44100);
idx = opts.idx;
fs = opts.fs;
deltat = 1/fs;

switch param
    case 'dprime'
        out = mean(ir2(:,idx))-mean(ir1(:,idx));
    case 'cross-correlation'
        out = deltat*sum( ir2(:,idx).*ir1(:,idx) );
    case 'cross-correlation-non-normalised' % as implemeted in AMT, it is not robust to fs
        out = optimaldetector(ir2(:,idx),ir1(:,idx));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
