function [h ha] = Get_waveforms_and_F0_praat(f1,f2,options,stPlot)
% function [h ha] = Get_waveforms_and_F0_praat(f1,f2,options,stPlot)
%
% 1. Description:
% 
% 2. Stand-alone example:
% 
% 3. Additional info:
%   Tested cross-platform: Yes
% 
% Inspired in script: VoD_read_aligned;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Original file name: VoD_read_aligned
% Created on    : 21/01/2015
% Last update on: 18/02/2015 % Update this date manually
% Last use on   : 18/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    stPlot = [];
end 

if nargin < 3
    options = [];
end 

options = Ensure_field(options,'type',4);
options = Ensure_field(options,'label1','input-1');
options = Ensure_field(options,'label2','input-2');

type = options.type;

[y1 fs] = Wavread(f1);
[y2 fs] = Wavread(f2);
L = min( length(y1), length(y2));

options = ef(options,'tanalysis',[0 (0:L-1)/fs]);

y1 = y1(1:L);
y2 = y2(1:L);
t = (0:L-1)/fs + options.tanalysis(1);

f1praat = Delete_extension(f1,'wav');
f2praat = Delete_extension(f2,'wav');

try
    [tfm f0m] = Get_F0_praat_from_txt( [f1praat '.txt.praat'] );
    disp('Praat file found')
catch
    [tfm, f0m] = Get_F0_AC_praat([f1praat '.wav'],[f1praat '.txt.praat'],type);
end

try
    [tfp f0p] = Get_F0_praat_from_txt( [f2praat '.txt.praat'] );
    disp('Praat file found')
catch
    [tfp, f0p] = Get_F0_AC_praat([f2praat '.wav'],[f2praat '.txt.praat'],type);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = [];
[h(end+1) ha(1:2)] = Plot_f0_waveform(t,y1,y2,options,stPlot);
[h(end+1) ha(3:4)] = Plot_fundamental_frequency(tfm,f0m,tfp,f0p,options,stPlot);

linkaxes(ha,'x');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end