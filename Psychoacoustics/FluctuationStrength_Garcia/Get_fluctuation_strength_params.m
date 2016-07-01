function params = Get_fluctuation_strength_params(N,fs,dataset)
% function params = Get_fluctuation_strength_params(N,fs,dataset)
% 
% Creates the file 'params.mat' (if the file already exists it is deleted
% first) containing all the required parameters for the fluctuation
% strength model.
% 
%   - gzi was deleted on 07/01/2015
%   - Ndataset = 0; is the final fitting, as presented in the thesis
% 
% Author: Rodrigo Garcia
% Original name: Create_params (renamed when copied from Rodrigo's repository)
% Modified by: Alejandro Osses V.
% Last modified: 07/01/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    dataset = 0; % 0 = Approved version
end
if nargin < 2
    fs = 44100;
end
if nargin < 1
    N = 2*fs; % 2 seconds
end

params         = struct;
params.fs      = fs;
params.N       = N;
params.Chno    = 47;
params.debug   = 'none';

switch dataset
    case 0
        params.window_type = 'cosine';
        params.filterbank = 'terhardt'; 
        % params.p_g     = 1; warning('pg temporarily to 1')
        params.p_g     = 1; 
        params.p_m     = 1.7;  
        params.p_k     = 1.7; % warning('Temporal value') 
        params.a0_in_time = 1;
        params.a0_in_freq = ~params.a0_in_time;
        
        params.cal     = 0.2490; % on 15/06/2016
        params.bIdle = 1; % v5
        warning('loading FS, v5, set bIdle back to 0 for v4');
        
    case 1
        params.window_type = 'blackman';
        params.filterbank = 'terhardt'; 
        params.p_g     = 2;
        params.p_m     = 0.25;
        params.p_k     = 0.375;
        params.cal     = 0.15/1.116;
        params.a0_in_time = 0;
        params.a0_in_freq = ~params.a0_in_time;
        
    case 90
        params.window_type = 'cosine';
        params.filterbank = 'erb';
        params.p_g     = 1;   
        params.p_m     = 1.7;  
        params.p_k     = 1.7; 
        params.cal     = 0.0419;
        params.a0_in_time = 0;
        params.a0_in_freq = ~params.a0_in_time;
        
    case 99
        params.window_type = 'cosine';
        params.filterbank = 'terhardt'; 
        params.a0_in_time = 0;
        params.a0_in_freq = ~params.a0_in_time;
        
end

params.dataset = dataset;

params.Hweight = Get_Hweight_fluctuation(fs);
params.gzi     = Calculate_gzi(params.Chno,dataset);