function [Fc Pyyw extra] = One_third_OB_analysis(x,fs, BandsPerOctave,info)
% function [Fc Pyyw extra] = One_third_OB_analysis(x,Fs, BandsPerOctave,info)
%
% 1. Description: one-third octave band (OB) analysis for an audio file
%       Each column corresponds to a temporal serie
% 
%       BandsPerOctave = 1 - OB analysis
%                        3 - one-third OB analysis
%       Pyyw - band power
% 
% % Example 1 (loading from an audio file):
%       filename = 'Choice.wav';
%       One_third_OB_analysis(filename);
% 
% % Example 2 (reading x vector and Fs):
%       BandsPerOctave = 1;
%       [Fc Py] = One_third_OB_analysis(x,Fs,BandsPerOctave);
% 
% % See file ../MATLAB_svn/Meas/Experiments/experiment_report_20140305_OB_analysis.m
%
% Programmed by Alejandro Osses, TUe
% Created on : 06/03/2014
% Last update: 09/05/2014
% Last used  : 09/05/2014 (Linux)
%              03/06/2014 (Win7)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
disp([mfilename '.m: try to use OB_analysis.m script instead'])

try
    filename = x;
    [x, fs] = wavread(filename);
    if length(x) > fs
        x = x(1:fs); 
        warning('Signal truncated to just 1 second for analysis')
    end
catch
    if nargin < 2
        error('You have to specify a sampling frequency fs [Hz]')
    end
end

if nargin < 4
    info = [];
end

info = Ensure_field(info,'f0min',0);
info = Ensure_field(info,'f0max', fs/2);

if nargin < 3
    BandsPerOctave = 3;
end

extra.rms = rmsdb(x);

disp(['RMS value: ' num2str(extra.rms)]);

N   = 8;           % Filter Order
F0  = 1000;       % Center Frequency (Hz)
f   = fdesign.octave(BandsPerOctave,'Class 1','N,F0',N,F0,fs);

F0  = validfrequencies(f);

F0 = F0(find(F0 < fs/2)); % Extra check
F0 = F0(find(F0 < info.f0max));
F0 = F0(find(F0 > info.f0min));

Nfc = length(F0);

for i=1:Nfc,
    f.F0 = F0(i);
    Hd3(i) = design(f,'butter');
end

Nx = length(x);

yw  = zeros(Nx,Nfc);
Pyyw = nan(Nfc,1);
for i=1:Nfc,
    yw(:,i) = filter(Hd3(i),x);
    
    Pyyw(i) = rmsdb(yw(:,i));
    
end

for i=size(yw,2):-1:1
    if sum(isnan(yw(:,i))) ~= 0
        disp('We are going to delete one column (there are NaN values)')
        yw(:,i) = [];
        F0(i) = [];
        Nfc = Nfc - 1;
        Pyyw(i) = [];
        
    end
end

Fc = F0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end