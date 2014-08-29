function [b,a]=ch_int_tp(f_abt);
% function [b,a]=ch_int_tp(f_abt);
%
% 1. Description:
%       Models temporal integration of loudness with simple 1st order low 
%       pass filter [butter 8Hz] 
%       Filter designed to match experimental data on pure tone loudness 
%       integration and fluctuation strength. 
%       
%       AO: original Chalupper's function had as input parameter
%               input 1: n      - signal to be smoothed
%               input 2: f_abt  - frequency related to hop-size
%               output 1: n_t   - smoothed (filtered) signal
% 
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.1.2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Altered by MFFM : Matt Flax is flatmax for the Psy-Sound project
% Jan. 2007

if nargin == 0
   f_abt = 500; 
end   

[b,a]=butter(1,8/(f_abt/2));
%n_t=filter(b,a,n);
%y=find(n_t <0);
%n_t(y) =0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
