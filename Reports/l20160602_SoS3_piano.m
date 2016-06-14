function l20160602_SoS3_piano
% function l20160602_SoS3_piano
%
% 1. Description:
%       Define below (around Line 28-30) the attenuation curve you want to apply
%       by specifying Rx in dB.
%       The sound that will be attenuated is the one defined by 'file' (adjust 
%       the file name accordingly).
% 
% 2. Stand-alone example:
%       l20160602_SoS3_piano;      
% 
% 3. Additional info:
%       Tested cross-platform: No
%       Functions needed: From_dB, Wavread.
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 22/04/2016
% Last update on: 22/04/2016 
% Last use on   : 02/06/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir        = 'D:\Documenten-TUe\01-Text\70-Presentaties-TUe\20160424-SoS3\';
diraudio   = [dir 'Audio' delim];
file       = [diraudio 'T8-sample-Klavier_orig.wav']; % file = [dir 'T3-sample-Klavier_orig.wav'];

FontSize = 14;

% Write down manually the corresponding attenuation curve defined by Rx:
fc = [63  80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150];
Rx = [22  22  23  27  30  34  38  45  51  56  60  59   61   67   71   70   68   67];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(1:length(fc),Rx,'o--','LineWidth',2), grid on, hold on
ha = gca;
h = gcf;
set(ha,'XTick',1:length(fc));
set(ha,'XTickLabel',fc);
xlabel('Frequency [Hz]');
ylabel('Reduction [dB]');

set(h,'PaperPositionMode', 'auto')
Poscur = get(h,'Position'); % 403 246 560 420
set(h,'Position',[Poscur(1) Poscur(2) 2*Poscur(3) Poscur(4)]);
set(ha,'FontSize',FontSize);
leg1 = 'Rx';
legend(leg1,'Location','SouthEast')
% Saveas(h,[dirfigures 'Rx'],'emf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[insig fs] = Wavread(file);
B = il_rx_piano(fc,Rx,fs); 
N = length(B); % filter order

outsig    = filter(B   ,1,insig);

fileout = [diraudio 'T8-sample-Klavier-out-Att-real.wav'];
Wavwrite(outsig   ,fs,fileout);

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = il_rx_piano(f,Rx,fs)

f = [0 f fs/2];
Rx = -1*[Rx(1) Rx Rx(end)];

N = 2^12;

B = fir2(N,f/(fs/2),From_dB(Rx));
[H1,Fn]=freqz(B,1);

% figure;
% semilogx(Fn, 20*log10(abs([H1])));
% 
% xlabel('Frequency [Hz]')
% legend([num2str(N) ' taps']);
% title('FIR filter to be used as approximation to isolation curve')

disp('')