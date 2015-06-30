function [val, idx] = subharm_detect(p, acf, val, idx)

% function [val, idx] = subharm_detect(p, acf, val, idx)
%
% Sub_harm_detect: Function that detects dominant subharmonics of the 
% originally detected fundamental frequency F0
%
% Programmed by Matthias, edited by AOV

if nargin == 0
    return;
end

harmonic = 1;
tmpIdx  = idx;
dpos = p.delta_lag_subh;
while tmpIdx/(harmonic + 1) > p.min_lag    
    harmonic = harmonic + 1;
    lag = tmpIdx/harmonic;
    pos  = round(lag);
    % Index is local thus ...
    if((pos - dpos) >= 1)
        [maxVal, maxIdx] = max(acf(pos - dpos:pos + dpos));
    else
        continue;
    end
    % ... (re)locate index in original autocorrelation array
    maxIdx = pos - dpos + maxIdx - 1;
   
    if maxVal > p.subharm_thr
        val = maxVal;
        idx = maxIdx;
    end
end
