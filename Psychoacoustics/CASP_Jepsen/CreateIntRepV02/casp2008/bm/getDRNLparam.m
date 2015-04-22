function [linDRNLparOut,NlinDRNLparOut] = getDRNLparam(CF);
% script for getting the (DRNL) filter parameters (Lopez-Poveda, meddis 2001)
% Author: Morten Løve Jepsen, 2.nov 2005, rev. 15 feb 2006, 19 feb 2007
%
% usage:  [linDRNLparOut,nlinDRNLparOut] = getDRNLparam(CF);
%
% The returned DRNLparam is a strucures containing the parameter name and
% values
% The parameters are from the article "A human nonlinear cochlear
% filterbank", table II, and table III, and modified in 2007

% revised version(*)
%   31/05/2011, C.T. Iben
% (*) set of parameters differs from the one published in Jepsen et al.
% (2008), but it reproduces the iso-intensity curve (Fig 2, p.425); after
% correspondence with Morten we agreed on this parameter set

%% DRNL for normal hearing, Jepsen 2007
% init structure
linDRNLstruct  = struct('parname',{},'vals',{}); % initialize linear paramater vector
NlinDRNLstruct = struct('parname',{},'vals',{}); % initialize nonlinear paramater vector
linDRNLparOut  = struct('parname',{},'vals',{});
NlinDRNLparOut = struct('parname',{},'vals',{});

linDRNLstruct(1).parname = 'CF_lin';
linDRNLstruct(2).parname = 'nGTfilt_lin';
linDRNLstruct(3).parname = 'BW_lin';
linDRNLstruct(4).parname = 'g';
linDRNLstruct(5).parname = 'LP_lin_cutoff';
linDRNLstruct(6).parname = 'nLPfilt_lin';
linDRNLparOut=linDRNLstruct;

NlinDRNLstruct(1).parname = 'CF_nlin';
NlinDRNLstruct(2).parname = 'nGTfilt_nlin';
NlinDRNLstruct(3).parname = 'BW_nlin';
NlinDRNLstruct(4).parname = 'a';
NlinDRNLstruct(5).parname = 'b';
NlinDRNLstruct(6).parname = 'c';
NlinDRNLstruct(7).parname = 'LP_nlin_cutoff';
NlinDRNLstruct(8).parname = 'nLPfilt_nlin';
NlinDRNLparOut=NlinDRNLstruct;



%% DRNL for normal hearing, Jepsen et al. (2008)
% assign values
linDRNLstruct(1).vals = 10^(-0.06762+1.01679*log10(CF)); % Hz, CF_lin
% 06/04/2011 CI: modified (this is 3 in Lopex-Poveda and Meddis (2001), table II)
linDRNLstruct(2).vals = 2; % number of cascaded gammatone filters
linDRNLstruct(3).vals = 10^(.03728+.75*log10(CF)); % Hz, BW_lin
linDRNLstruct(4).vals = 10^(4.20405 -.47909*log10(CF)); %g
linDRNLstruct(5).vals = 10^(-0.06762+1.01*log10(CF)); % Hz, LP_lin cutoff
linDRNLstruct(6).vals = 4; % no. of cascaded LP filters
NlinDRNLstruct(1).vals = 10^(-0.05252+1.01650*log10(CF)); % Hz, CF_nlin
% 06/04/2011 CI: modified (this is 3 in Lopex-Poveda and Meddis (2001), table II)
NlinDRNLstruct(2).vals = 2; % number of cascaded gammatone filters
% 06/04/2011 CI: modified version from Morten (this is +.77*log10(CF)
% in Lopex-Poveda and Meddis (2001), table II)
NlinDRNLstruct(3).vals = 10^(-0.03193+.70*log10(CF)); % Hz, BW_nlin
if CF<=1000
    % SE 03.02.2011, the 2008 paper states <= 1500 Hz
    % 06/04/2011 CI: answer from Morten regarding the discontinuity:
    % This is imprecisely described in the paper. It was simulated as
    % described with parameter a, having the value for 1500 Hz, for CFs
    % above 1000 Hz. I do recognize the discontinuity in the derived
    % parameter, but I think this is not critical
    NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(CF)); % a, the 1500 assumption is no good for compressionat low freq filters
    NlinDRNLstruct(5).vals = 10^(1.61912-.81867*log10(CF)); % b [(m/s)^(1-c)]
else
    NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(1500)); % a, the 1500 assumption is no good for compressionat low freq filters
    NlinDRNLstruct(5).vals = 10^(1.61912-.81867*log10(1500)); % b [(m/s)^(1-c)]
    % SE 03.02.2011, the current values are
    % NlinDRNLstruct(4).vals equals 4.0047069 (same as in paper)
    % NlinDRNLstruct(5).vals equals -0.98105063 (paper says -0.98015)
    % typo in paper?
end
NlinDRNLstruct(6).vals = 10^(-.60206); % c, compression coeff
% 06/04/2011 CI: modified (this is 1.01650 in Lopex-Poveda and Meddis
% (2001), table II)
NlinDRNLstruct(7).vals = 10^(-0.05252+1.01*log10(CF)); % LP_nlincutoff
% 31/05/2011 CI: modified (this is 3 in Lopex-Poveda and Meddis (2001),
% table II)
NlinDRNLstruct(8).vals = 1; % no. of cascaded LP filters in nlin path


for k=1:6
    linDRNLparOut(k).vals = linDRNLstruct(k).vals;
end
for k=1:8
    NlinDRNLparOut(k).vals = NlinDRNLstruct(k).vals;
end

%eof

