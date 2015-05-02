function [outsigs, fc, t, opts] = Get_internal_representations_deterministic(insigs,fs,model,opts)
% function [outsigs, fc, t, opts] = Get_internal_representations_deterministic(insigs,fs,model,opts)
%
% 1. Description:
%       LvNoise - level of the internal noise
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
% Last update on: 30/04/2015 % Update this date manually
% Last use on   : 30/04/2015 % Update this date manually
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
insig2  = insigs(:,2);
insig3  = insigs(:,3);
try
    insig4  = insigs(:,4);
end
    
if strcmp(model,'dau1996a') % No overshoot limit
    [outsig  fc] = dau1996apreproc_1Ch(insig ,fs,fc);
    [outsig2 fc] = dau1996apreproc_1Ch(insig2,fs,fc);
    [outsig3 fc] = dau1996apreproc_1Ch(insig3,fs,fc);
    try
        [outsig4 fc] = dau1996apreproc_1Ch(insig4,fs,fc);
    end
end

if strcmp(model,'dau1996') % Overshoot limit
    [outsig  fc] = dau1996preproc_1Ch(insig,fs,fc);
    [outsig2 fc] = dau1996preproc_1Ch(insig2,fs,fc);
    [outsig3 fc] = dau1996preproc_1Ch(insig3,fs,fc);
    try
        [outsig4 fc] = dau1996preproc_1Ch(insig4,fs,fc);
    end
end

% if strcmp(model,'dau1997') 
%     [outsig fc mf] = dau1997preproc_1Ch(insig,fs,fc);
% end

opts    = Ensure_field(opts,'bAddNoise',1);
opts    = Ensure_field(opts,'sigma',0);

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
        outsig2 = outsig2 + yn;
        outsig3 = outsig3 + yn;
        try
            outsig4 = outsig4 + yn;
        end
        
    end
    
    fprintf( 'Internal noise with mean = %.0f, std = %.2f [MU]', mean(yn), std(yn) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outsigs = [outsig outsig2 outsig3];
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
