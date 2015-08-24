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
%       Tested cross-platform: Yes
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

opts = Ensure_field(opts,'fmin',fc);
opts = Ensure_field(opts,'fmax',fc);
fmin = opts.fmin;
fmax = opts.fmax;

insig   = insigs(:,1);
try
    insig2  = insigs(:,2);
    insig3  = insigs(:,3);
    insig4  = insigs(:,4);
end
    
if strcmp(model,'dau1996a') % No overshoot limit
    
    if fmin ~= fc | fmax ~= fc
        error('Single-channel model being used, update this script and re-run the analysis');
    end
    
    [outsig  fc] = dau1996apreproc_1Ch(insig ,fs,fc);
    try
        [outsig2 fc] = dau1996apreproc_1Ch(insig2,fs,fc);
        [outsig3 fc] = dau1996apreproc_1Ch(insig3,fs,fc);
        [outsig4 fc] = dau1996apreproc_1Ch(insig4,fs,fc);
    end
end

if strcmp(model,'dau1996') % Overshoot limit
    
    if fmin ~= fc | fmax ~= fc
        bSingleChannel = 0;
        bMultiChannel = 1;
    else
        bSingleChannel = 1;
        bMultiChannel = 0;
    end
    
    if bSingleChannel == 1
        [outsig  fc] = dau1996preproc_1Ch(insig,fs,fc);
        try
            [outsig2 fc] = dau1996preproc_1Ch(insig2,fs,fc);
            [outsig3 fc] = dau1996preproc_1Ch(insig3,fs,fc);
            [outsig4 fc] = dau1996preproc_1Ch(insig4,fs,fc);
        end
    end
    
    if bMultiChannel
        [outsig  fc] = dau1996preproc(insig,fs);
        try
            [outsig2 fc] = dau1996preproc(insig2,fs);
            [outsig3 fc] = dau1996preproc(insig3,fs);
            [outsig4 fc] = dau1996preproc(insig4,fs);
        end
        
        [xx,idx] = find(fc>=fmin & fc<=fmax);
        fc = fc(idx);
        
        outsig = outsig(:,idx);
        try
            outsig2 = outsig2(:,idx);
            outsig3 = outsig3(:,idx);
            outsig4 = outsig4(:,idx);
        end
        
    end
end

if strcmp(model,'dau1997') 
    
    if fmin ~= fc | fmax ~= fc
        error('Single-channel model being used, update this script and re-run the analysis');
    end
    
    [outsig fc mf] = dau1997preproc_1Ch(insig,fs,fc);
end

if strcmp(model,'jepsen2008') 
    
    % fc2look = fc;
    [outsig  fc  mfc] = jepsen2008preproc(insig,fs);
    [outsig  fc  mfc] = jepsen2008preproc(insig,fs,'resample_intrep');
    try
        [outsig2 fc mfc] = jepsen2008preproc(insig2,fs,'resample_intrep');
        [outsig3 fc mfc] = jepsen2008preproc(insig3,fs,'resample_intrep');
        [outsig4 fc mfc] = jepsen2008preproc(insig4,fs,'resample_intrep');
    end

    [xx,idx] = find(fc>=fmin & fc<=fmax);
    fc = fc(idx);
    outsigtmp = [];
    for i = 1:length(idx)
        outsigtmp = [outsigtmp outsig{idx(i)}(:,:)];
        
    end
    outsig = outsigtmp;
    
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
