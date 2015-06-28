function l7S690_20150621_RA_Computational
% function l7S690_20150621_RA_Computational
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 21/06/2015
% Last update on: 21/06/2015 % Update this date manually
% Last use on   : 21/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

bDiary = 0;
Diary(mfilename,bDiary);

dir_files = 'D:\Documenten-TUe\09-Training+activities\2015-Q3-Architectural-acoustics\Experiments-2\model\txt-exported-4-MATLAB\';
% fname = [dir_files 'RAdBcommas.mat']; load(fname);
% fname = [dir_files 'RAdBcommas2.mat']; load(fname); RAdBcommas = RAdBcommas2; delete RAdBcommas2;
% fname = [dir_files 'RAdBcommas3.mat']; load(fname); RAdBcommas = RAdBcommas3; delete RAdBcommas3;
fname = [dir_files 'RAdBcommas4.mat']; load(fname); RAdBcommas = RAdBcommas4; delete RAdBcommas4;

f = RAdBcommas(:,1);
lvl = RAdBcommas(:,2:end);

pref = 2e-5;
pr = pref*10.^(lvl/20);

pravg = transpose( mean(transpose(pr)) );
lvlavg = 20*log10( pravg/pref );

idx = find( f<=98 & f>=97 );
idxinf = idx(1:2);
idxsup = idx(2:3);
lvlref  = lvlavg(idx(2));
f0      = f( idx(2) );
lvlm3dB = lvlref - 3;
finf    = interp1(lvlavg(idxinf),f(idxinf),lvlm3dB);

fsup    = interp1(lvlavg(idxsup),f(idxsup),lvlm3dB);
w0      = 2*pi*f0;
winf    = 2*pi*finf;
wsup    = 2*pi*fsup;
Dw      = wsup - winf;
Q       = w0/Dw;
delta   = w0/(Q*2);

T = 3*log(10)/delta;    % 5.5911 s                        
                        % 3.3699 s with Z = 100000 Pa s/m
                        % 4.6879 s with Z = 200000 Pa s/m
                        % 4.9508 s with Z = 250000 Pa s/m

figure;
plot(f,lvlavg), grid on
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
hold on
plot(f(idx),lvlavg(idx),'LineWidth',2)

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
