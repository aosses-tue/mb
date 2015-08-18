function [outsigs, fc, t, opts] = Get_internal_representations_deterministic(insigs,fs,model,opts)
% function [outsigs, fc, t, opts] = Get_internal_representations_deterministic(insigs,fs,model,opts)
%
% 1. Description:
%       LvNoise - level of the internal noise
%       Same noise is added to all input signals specified by insigs
% 
% 2. Stand-alone example:
%       [insig t]   = Create_sin; % default sine tone
%       fs          = 44100;
%       var         = 0;
%       model       = 'dau1996';
%       Get_internal_representations_deterministic(insig,fs,model,var);
% 
% 3. Additional info:
%       Tested cross-platform: No
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Original file : Get_internal_representations.m
% Created on    : 30/04/2015
% Last update on: 18/08/2015 
% Last use on   : 18/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4
    opts = [];
end

if nargin < 3
    model = 'dau1996';
end

if nargin < 2
    error('Specify the sample rate of insig');
end

opts    = Ensure_field(opts,'fc',1000);
fc      = opts.fc;

insig   = insigs(:,1);
try
    insig2  = insigs(:,2);
    insig3  = insigs(:,3);
    insig4  = insigs(:,4);
end
    
if strcmp(model,'dau1996a') % No overshoot limit
    [outsig  fc] = dau1996apreproc_1Ch(insig ,fs,fc);
    try
        [outsig2 fc] = dau1996apreproc_1Ch(insig2,fs,fc);
        [outsig3 fc] = dau1996apreproc_1Ch(insig3,fs,fc);
        [outsig4 fc] = dau1996apreproc_1Ch(insig4,fs,fc);
    end
end

if strcmp(model,'dau1996') % Overshoot limit
    [outsig  fc] = dau1996preproc_1Ch(insig,fs,fc);
    try
        [outsig2 fc] = dau1996preproc_1Ch(insig2,fs,fc);
        [outsig3 fc] = dau1996preproc_1Ch(insig3,fs,fc);
        [outsig4 fc] = dau1996preproc_1Ch(insig4,fs,fc);
    end
end

if strcmp(model,'dau1997') 
    [outsig fc mf] = dau1997preproc_1Ch(insig,fs,fc);
end

if strcmp(model,'jepsen2008') 
    fc2look = fc;
    [outsig  fc  mfc] = jepsen2008preproc(insig,fs);
    try
        [outsig2 fc mfc] = jepsen2008preproc(insig2,fs);
        [outsig3 fc mfc] = jepsen2008preproc(insig3,fs);
        [outsig4 fc mfc] = jepsen2008preproc(insig4,fs);
    end
    [xx,idx] = max(find(fc<=fc2look));
    fc = fc(idx);
    outsig = outsig{idx}(:,:);
end

opts    = Ensure_field(opts,'bAddNoise',1);

if opts.bAddNoise == 1
    opts    = Ensure_field(opts,'sigma',0);
else
    opts    = ef(opts,'sigma',0);
end

bAddNoise   = opts.bAddNoise;
sigma       = opts.sigma;
mu          = 0;

t = ( 1:size(outsig,1) )/fs;

N = size(outsig,1);
M = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Internal noise:

if bAddNoise
    for i = 1:M
        
        yn = normrnd(mu,sigma,N,1);
        outsig  = outsig  + yn;
        try
            outsig2 = outsig2 + yn;
            outsig3 = outsig3 + yn;
            outsig4 = outsig4 + yn;
        end
        
    end
    
    fprintf( 'Internal noise with mean = %.0f, std = %.2f [MU]', mean(yn), std(yn) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outsigs = outsig;
try
    outsigs = [outsigs outsig2];
end
try
	outsigs = [outsigs outsig3];
end
try
    outsigs = [outsigs outsig4];
end
    
if nargout == 0
    
    figure; 
    plot(outsig), hold on; 
    plot(outsig2,'g'); 
    plot(outsig3,'r')
    
    legend('1','2','3')
    title(sprintf('Channel tuned to fc = %.2f [Hz]',fc));
    
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
