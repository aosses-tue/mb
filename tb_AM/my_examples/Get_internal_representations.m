function [outsig, fc, t, opts] = Get_internal_representations(insig,fs,model,var)
% function [outsig, fc, t, opts] = Get_internal_representations(insig,fs,model,var)
%
% 1. Description:
%       var - sigma of the internal noise
%       opts.idx - indexes of bands with larger average values (largest, second largest, third largest)
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
% Last update on: 30/10/2014 % Update this date manually
% Last use on   : 30/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4
    var = 1e-5;
end

if nargin < 3
    model = 'dau1996';
end

if nargin < 2
    error('Specify the sample rate of insig');
end

if strcmp(model,'dau1996a') % No overshoot limit
    [outsig fc] = dau1996apreproc(insig,fs);
end

if strcmp(model,'dau1996') % Overshoot limit
    [outsig fc] = dau1996preproc(insig,fs);
end

opts = [];

t = ( 1:size(outsig,1) )/fs;

N = size(outsig,1);
M = size(outsig,2);

var = max(var,1e-5);

yn = wgn(N,1, To_dB(var) );

outsig = outsig + repmat(yn,1,M);


idx2avg = size(outsig,1)/10; % the last tenth of the array
tmp = outsig(end-idx2avg+1:end,:);
tmp = rms(tmp);
idx(1) = find(  tmp == max( (tmp) )  );
tmp(idx(1)) = 0;

idx(2) = find(  tmp == max( (tmp) )  );
tmp(idx(2)) = 0;

idx(3) = find(  tmp == max( (tmp) )  );

opts.idx = idx;

if nargout == 0
    
    figure; 
    plot(t,outsig(:,idx)); grid on
    legend(     num2str(fc(idx(1))), ...
                num2str(fc(idx(2))), ...
                num2str(fc(idx(3))) );
    
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
