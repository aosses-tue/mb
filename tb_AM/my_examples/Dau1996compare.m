function outs = Dau1996compare(insigM,insigS,fs,opts,method)
% function outs = Dau1996compare(insigM,insigS,fs,opts,method)
%
% 1. Description:
%       method can be: 'all','difference','template','internal-representations'
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
% Last update on: 29/09/2015 
% Last use on   : 29/09/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = ( 1:length(insigM) )/fs;
t = t(:);

if nargin < 5
    method = 'all';
end

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
opts = ef(opts,'bAddNoise',0);

bAddNoise = opts.bAddNoise;
if bAddNoise
    opts = ef(opts,'mu',0);
    opts = ef(opts,'sigma',0.8);
end

bPlot = opts.bPlot;
 

%% Processing 

insigMS = insigM + insigS;

exp1 = sprintf('[out1, fc , outsig1] = %spreproc(insigM ,fs);',opts.method);
exp2 = sprintf('[out2, fc , outsig2] = %spreproc(insigMS,fs);',opts.method);
eval(exp1);
eval(exp2);

idx     = max(find(fc<opts.fc_idx));
fcentre = fc(idx);

switch opts.calc_method
    case 1
        tmpmethod = 0;
    case 2
        tmpmethod = 2;
end

if bAddNoise
    out1(:,idx) = Add_gaussian_noise(out1(:,idx),opts.mu,opts.sigma);
    out2(:,idx) = Add_gaussian_noise(out2(:,idx),opts.mu,opts.sigma);
end

curr_diff = out2(:,idx)-out1(:,idx);
switch method
    case {'all','difference'}
        outs.curr_diff = curr_diff;
    case {'all','template'}
        % Normalisation of the template:
        outs.template= Normalise_signal(curr_diff,fs,tmpmethod);
    case {'all','internal-representations'}
        outs.outsig1 = outsig1;
        outs.outsig2 = outsig2;
end

outs.idx     = idx; % plotted/to be plotted idx
outs.t       = t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
