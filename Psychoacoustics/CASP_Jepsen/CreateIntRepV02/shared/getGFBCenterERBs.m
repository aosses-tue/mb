function [info, cerbs] = getGFBCenterERBs(lowf, uppf, baseF, density)
% function [info, cerbs] = getGFBCenterERBs(lowf, uppf, baseF, density)
% 
% gammaFB.m - applies gammatone auditory filterbank
%
% Usage: [inf, cerbs] = getGFBCenterERBs(lowf, uppf, density);
%
% in            = input vector
% lowf          = lower frequency boundary for filters
% uppf          = upper frequency boundary for filters
%                 (if lowf == uppf, only one filter is chosen)
% fs            = sampling rate in Hz
% alignfreq     : alignment frequency; corresponds to the center frequency 
%                 of the first filter set
% 
% out           = output array (columns = channels, rows = time)
% cfs           = vector of center frequencies in ERB (use erbtofreq.m 
%                 to convert to Hz)

% modified version
    % 09/01/12, C.T. Iben

% 09/01/12, ci: alignment of a filter around a certain frequency baseF = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp = freqtoerb([lowf baseF uppf]);
ci = freqtoern(baseF); % added by AO

if ( tmp(1) == tmp(3) )

	cerbs = tmp(1);

else 
	
	tmp2    = [0:density:100];
	tmp2    = [-fliplr(tmp2) tmp2(2:end)]+ci;
	
	i_start = min(find(tmp2>=tmp(1)));
	i_end   = max(find(tmp2<=tmp(3)));

	cerbs   = tmp2(i_start:i_end);
end

info = length(cerbs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
