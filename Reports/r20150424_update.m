function r20150424_update
% function r20150424_update
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 22/04/2015
% Last update on: 22/04/2015 % Update this date manually
% Last use on   : 22/04/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

[insig,fs] = greasy;
insig = resample(insig,44100,fs);
fs = 44100;
[outsig, fc, mfc] = jepsen2008preproc(insig, fs);

N = length(insig);
M = length(fc);
IntRep.Y = zeros(N,M,length(mfc)); % memory allocation
for i = 1:M
    P = size( outsig{i,1}, 2);
    IntRep.Y(1:N,i,1:P) = outsig{i,1}(1:N,1:P);
end

PlotIntRepImage(fc,mfc,IntRep);








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path1 = [Get_TUe_paths('db_voice_of_dragon') '02-Wav-files'           delim '2015-02-wav-files' delim '02-calibrated' delim 'meas-ac-4-dist-ane.wav'];
path2 = [Get_TUe_paths('db_voice_of_dragon') '03-Wav-files-predicted' delim '2015-02-wav-files' delim '02-calibrated' delim 'model-ac-4-dist-ane.wav'];

file1 = path1;
file2 = path2;

[x1 fs1] = Wavread(file1);
[x2 fs2] = Wavread(file2);

if fs1 == fs2
    fs = fs1;
end

ti = 0.124;
tf = 1.202;

Ni = round(ti*fs);
Nf = round(tf*fs);

x1 = x1(Ni:Nf);
x2 = x2(Ni:Nf);

model = 'casp';
pmode = 'IntRepImage'; %'all';
DemoMain(model,pmode,x1);

DemoMain(model,pmode,x2);

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])

function PlotIntRepImage(freq1,freq2,IntRep)

if nargin < 3, error('not enough input arguments'); end

% calculation of rms across time in every auditory and modulation filter
IntRepImage = squeeze(sqrt(mean((IntRep.Y).^2)));

IntRepImage(1,1)=300; % 220
warning('Temporal variable assigned by AO');

if length(freq1) > 1 && length(freq2) > 1
    figH7=figure;
    set(figH7,'Name','internal representation (energy is color encoded) across f_c and f_m','NumberTitle','off')
    imagesc(IntRepImage')
    colorbar
    axis xy
    set(gca,'XTick',1:3:length(freq1),'XTickLabel',freq1(1:3:end))
    set(gca,'YTick',1:length(freq2),'YTickLabel',freq2)
    xlabel('auditory filters f_c in Hz')
    ylabel('modulation filters f_m in Hz')
    title('internal representation (energy is color encoded) across f_c and f_m','Fontweight','bold')
    %     colorbar
else
    error('illegal plot modus: not available with one auditory or one modulation filter')
end


