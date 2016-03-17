function [y, returnData] = ch_rms_tep(fout, f_abt, fs)
% function [y, returnData] = ch_rms_tep(fout,f_abt,fs)
%
% Calculation of short term RMS levels with auditory temporal window 
% (Plack&Moore)
%
% erd=4ms, 
% w=-51 dB
% f_abt (envelope), default f_abt=500ms
% fs(time signal), default fs = 44100 Hz
%
%       y - short term RMS levels
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author        : Josef Chalupper (josef.chalupper@siemens.com)
% Created on    : 12/12/2000
% Edited on     : 06/01/2007 (new version with comments and examples)
% Downloaded on : 07/08/2014 (approx.)
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 29/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Altered by flatmax is Matt Flax for the Psy-Sound project
% Jan. 2007

if nargin <2
   f_abt=500;
   fs=44100;
end

if size(fout,2) < size(fout,1)
	fout = fout';
end

[t_pa,w,t_sb,t_sa,t_pb] = staticParamDLM;
[h, t, erd] = tep_window(t_pb, t_pa, t_sb, t_sa, w, fs);
h = fliplr(h.^2)'; %due to convolution and intensity; 

dauer   = erd*fs;
step    = round(1/f_abt*fs);
wlen    = length(h);
n_steps = floor((length(fout)-wlen)/step)+1;

rms   = [];
power = fout'.^2;
for i=0:(n_steps-1)
  rms(i+1) = sum(power(i*step+1:i*step+wlen) .* h)/dauer;
end
returnData = fout(n_steps*step:end);

y=sqrt(rms);
j=find(y==0);  
y(j) = realmin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end