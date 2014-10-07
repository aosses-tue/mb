function r20140919_update(options)
% function r20140919_update(options)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       options.bSave = 0;
%       r20140930_FIA_TUe(options);
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 17/09/2014
% Last update on: 17/09/2014 % Update this date manually
% Last use on   : 17/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Voice of the Dragon

opts = Get_VoD_params;
% Voice of the dragon
ceff = 310;
L = 0.7;
n = 2:5;
fn = n*ceff/(2*L);

var2latex([n' fn' opts.mf']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KUL

if nargin == 0
    options = [];
end

close all

bCreateWhitenoise = 0; % Whitenoise for LIST-f

options.dest_folder      = [Get_TUe_paths('outputs') 'tmp-Audio'  delim];
options.dest_folder_fig  = [Get_TUe_paths('outputs') 'tmp-Figure' delim];
options.nAnalyser = 10;

bDiary = 0;
Diary(mfilename,bDiary);

tic

tmp = Get_TUe_subpaths('db_speechmaterials');
root_folder = tmp.allfiles_LISTf;
filename = [root_folder 'wdz2.wav']; % 'Elke zaterdag ga ik naar de markt'
Wavread(filename);
filenoise1 = [root_folder 'whitenoise-LISTf.wav'];
filenoise2 = [root_folder 'wivineruis.wav'];

options.CalMethod = 5;
options.bPlot = 0;

if bCreateWhitenoise
    [x1 fs1] = Wavread([tmp.allfiles_LISTf 'wivineruis.wav']);
    x1RMS = rmsdb(x1);
    [x2 fs2] = Wavread([tmp.allfiles_PB 'whitenoise.wav']); % originally -10.78 dBFS RMS
    x2 = setdbspl(x2,x1RMS+100);
    rmsdb(x2)
    Wavwrite(x2,fs1,[tmp.allfiles_LISTf 'whitenoise-LISTf']);
end

out_file1 = PsySoundCL(filename,options);
out_noise1 = PsySoundCL(filenoise1,options);
out_noise2 = PsySoundCL(filenoise2,options);
f = out_file1.f;

figure;
semilogx(   f, out_file1.DataSpecOneThirdAvg, 'bo--', ...
            f, out_noise1.DataSpecOneThirdAvg-5, 'kx--', ...
            f, out_noise2.DataSpecOneThirdAvg-5,'r>-');
grid on
legend('wdz2','white noise, 60 dB SPL','SSN, 60 dB SPL')
xlabel('Frequency [Hz]')
ylabel('Sound Pressure Level [dB]')

% Saveas(gcf, 'figure','epsc') % As sent to TF, e-mail 18-09-2014

disp('');

toc

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
