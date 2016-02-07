function FS = Get_fluctuation_strength(filename,model_par)
% function FS = Get_fluctuation_strength(filename,model_par)
% 
% Returns fluctuation strength value using debug model.
% 
% Inputs:
% filename: The stimulus filename.
% m_p: The model parameters struct.
% 
% Outputs:
% fs: The fluctuation strength value.
% 
% Author: Rodrigo Garcia
%
% Example:
%   
%   filename = 'AM-fm_0.00-fc_1000-md_40-SPL_70-w_25-fs_44100-N_529200.wav';
%   [x fs] = Wavread(filename);
%   load('spec-20160203-1932.mat');
%   load('m_p-20160203-1933.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate stimulus parameters
[insig fs] = Wavread(filename);
N = model_par.N;
if length(insig)>N
    insig  = insig(1:N);
end
FS     = FluctuationStrength_Garcia(insig,fs,N,model_par);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
