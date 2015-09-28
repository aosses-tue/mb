function outs = Dau1996compare(insig1,insig2,fs,opts)
% function outs = Dau1996compare(insig1,insig2,fs,opts)
%
% 1. Description:
%
% 2. Stand-alone example:
%       options.bSave = 0;
%       demo_dau1996b(options);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 08/10/2014
% Last update on: 24/10/2014 % Update this date manually
% Last use on   : 24/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = [];

t = ( 1:length(insig1) )/fs;
t = t(:);

if nargin < 4
    opts = [];
    if nargout == 0
        opts.bPlot = 1;
    else
        opts.bPlot = 0;
    end
end

opts = Ensure_field(opts,'fc_idx',3000);
opts = Ensure_field(opts,'method','dau1996'); % dau1996  uses dau1996preproc
                                              % dau1996a uses dau1996apreproc
bPlot = opts.bPlot;
 

%% Processing 

insig2 = insig1 + insig2;

exp1 = sprintf('[out1, fc , outsig1] = %spreproc(insig1,fs);',opts.method);
exp2 = sprintf('[out2, fc , outsig2] = %spreproc(insig2,fs);',opts.method);
eval(exp1);
eval(exp2);

idx     = max(find(fc<opts.fc_idx));
fcentre = fc(idx);

haxis = [];

% Normalisation of the template:
% outs.template_no_norm   = outsig2.out04_LPF-outsig1.out04_LPF;
outs.template_no_norm   = out2(:,idx)-out1(:,idx);
outs.template           = Normalise_signal(outs.template_no_norm,fs);

outs.idx                = idx; % plotted/to be plotted idx

outs.t = t;
outs.outsig1 = outsig1;
outs.outsig2 = outsig2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
