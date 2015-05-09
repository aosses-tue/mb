function [outsig, fc, t, opts] = Get_internal_representations(insig,fs,model,opts)
% function [outsig, fc, t, opts] = Get_internal_representations(insig,fs,model,opts)
%
% 1. Description:
%       LvNoise - level of the internal noise
% 
% 2. Stand-alone example:
%       [insig t]   = Create_sin; % default sine tone
%       fs          = 44100;
%       var         = 0;
%       model       = 'dau1996';
%       Get_internal_representations(insig,fs,model,var);
% 
% 3. Additional info:
%       Tested cross-platform: No
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 22/10/2014
% Last update on: 29/04/2015 % Update this date manually
% Last use on   : 08/05/2015 % Update this date manually
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

if strcmp(model,'dau1996a') % No overshoot limit
    [outsig fc outint] = dau1996apreproc(insig,fs);
end

if strcmp(model,'dau1996') % Overshoot limit
    [outsig fc outint] = dau1996preproc(insig,fs);
end

if strcmp(model,'dau1997') 
    [outsig fc mf outint] = dau1997preproc(insig,fs);
end

% if strcmp(model,'jepsen2008') 
%     [outsig fc] = jepsen2008preproc(insig,fs);
% end

opts    = Ensure_field(opts,'bAddNoise',1);
if opts.bAddNoise == 1
    opts    = Ensure_field(opts,'sigma',10);
end

bAddNoise   = opts.bAddNoise;
sigma       = opts.sigma;
mu          = 0;

t = ( 1:size(outsig,1) )/fs;

N = size(outsig,1);
M = size(outsig,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Internal noise:
% var = max(var,1e-5);
% yn = wgn(N,1, To_dB(var) );

if bAddNoise
    for i = 1:M
        
        yn = normrnd(mu,sigma,N,1);
        var(i) = std(yn);
        IntNoise(:,i) = yn;
    end
    
    fprintf('Internal noise with RMS of %.2f [dB]\n', mean( rmsdb(IntNoise) ));
    outsig = outsig + IntNoise;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout == 0
    idx2avg = size(outsig,1)/10; % the last tenth of the array
else
    idx2avg = size(outsig,1); 
end
   
try
    tmp = outsig(end-idx2avg+1:end,:);
    tmp = rms(tmp);
    idx(1) = find(  tmp == max( (tmp) )  );
    tmp(idx(1)) = 0;

    idx(2) = find(  tmp == max( (tmp) )  );
    tmp(idx(2)) = 0;

    idx(3) = find(  tmp == max( (tmp) )  );

    opts.idx = idx;
catch
    warning('Look at this script...');
end

opts.outint = outint;

if nargout == 0
    
    figure; 
    plot(t,outsig(:,idx)); grid on
    legend(     num2str(fc(idx(1))), ...
                num2str(fc(idx(2))), ...
                num2str(fc(idx(3))) );
	title('Three frequency bands with more energy')
    
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
