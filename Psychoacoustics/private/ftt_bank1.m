function [output, S] = ch_ftt_bank1(x, S, f_abt,fs)
% function [output, S] = ch_ftt_bank1(x, S, f_abt,fs)
%
% 1. Description
%       Applies filterbank (as calculated in makefilterbank1.m) and envelope
%       extraction (with auditory temporal window). It requires ch_rms_tep.m 
%       and ch_tep_window.m
%
% Author: Josef Chalupper (josef.chalupper@siemens.com)
% original version: 12.12.2000
% new version (with comments and examples): 6.1.2007
% Edited by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Altered by Matt Flax is flatmax for the Psy-Sound project
% Jan. 2007

fcoefs = S.fcoefs;

if nargin<3
  f_abt = 500;
  fs    = 44100;
end

if size(fcoefs,2) ~= 10
  error('fcoefs parameter sind falsch.');
end

if size(x,2) < size(x,1)
  x = x';
end

x = x.*2*10^.5; %full scale = 107 dB SPL

A0  = fcoefs(:,1);
A11 = fcoefs(:,2);
A12 = fcoefs(:,3);
A13 = fcoefs(:,4);
A14 = fcoefs(:,5);
A2  = fcoefs(:,6);
B0  = fcoefs(:,7);
B1  = fcoefs(:,8);
B2  = fcoefs(:,9);
gain= fcoefs(:,10);

len = size(gain, 1);
if ~isfield(S, 'y4')
  S.y4 = cell(1, len);
end

for chan = 1:len
  [y1, S.Zf1{chan}] = filter(...
      [A0(chan)/gain(chan) A11(chan)/gain(chan) A2(chan)/gain(chan)], ...
      [B0(chan) B1(chan) B2(chan)], x, S.Zf1{chan});
  
  [y2, S.Zf2{chan}] = filter( [A0(chan) A12(chan) A2(chan)], ...
                              [B0(chan) B1(chan)  B2(chan)], ...
                              y1, S.Zf2{chan});
  
  [y3, S.Zf3{chan}] = filter([A0(chan) A13(chan) A2(chan)], ...
                             [B0(chan) B1(chan)  B2(chan)], ...
                             y2, S.Zf3{chan});
  
  [y4, S.Zf4{chan}] = filter([A0(chan) A14(chan) A2(chan)], ...
                             [B0(chan) B1(chan)  B2(chan)], ...
                             y3, S.Zf4{chan});
  
  S.y4{chan} = [S.y4{chan} y4];

  [temp, S.y4{chan}] = ch_rms_tep(S.y4{chan}, f_abt, fs);
  y4 = 20*log10(temp/2e-5);
  if chan==1
    output = zeros(length(y4),size(gain,1));
  end
  output(:,chan) = y4';
end

output = fliplr(output);

end