function [Pyyw Fc extra] = OB_analysis(x,fs, BandsPerOctave,info)
% function [Pyyw Fc extra] = OB_analysis(x,fs, BandsPerOctave,info)
%
% 1. Description: one-third octave band (OB) analysis for an audio file
%       Each column corresponds to a temporal serie.
%       Optimised for vector processing (compared to One_third_OB_analysis.m)
% 
%       BandsPerOctave = 1 - OB analysis
%                        3 - one-third OB analysis
%       Pyyw - band power
% 
% % Example 1 (reading x vector and Fs):
%       BandsPerOctave = 1;
%       [Py Fc] = OB_analysis(x,Fs,BandsPerOctave);
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 3/6/2014
% Last update: 3/6/2014 % Update this date manually
% Last used: 3/6/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
Ny = size(x,2);
Pyyw = nan(Nfc,Ny);
for i=1:Nfc,
    yw = filter(Hd3(i),x);
    Get_filter_specs(Hd3(i));
    Pyyw(i,:) = rmsdb(yw);
end

Fc = F0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end