function l20160422_SoS3_piano
% function l20160422_SoS3_piano
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 22/04/2016
% Last update on: 22/04/2016 
% Last use on   : 02/06/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir        = 'D:\Documenten-TUe\01-Text\70-Presentaties-TUe\20160424-SoS3\';
diraudio   = [dir 'Audio' delim];
dirfigures = [dir 'Figures' delim];
file       = [diraudio 'T8-sample-Klavier_orig.wav']; % file = [dir 'T3-sample-Klavier_orig.wav'];

FontSize = 14;

close all
fc = [63  80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150];
Rx = [22  22  23  27  30  34  38  45  51  56  60  59   61   67   71   70   68   67];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(1:length(fc),Rx,'o--','LineWidth',2), grid on, hold on
ha = gca;
h = gcf;
set(ha,'XTick',1:length(fc));
set(ha,'XTickLabel',fc);
Xlabel('Frequency [Hz]',FontSize);
Ylabel('Reduction [dB]',FontSize);

set(h,'PaperPositionMode', 'auto')
Poscur = get(h,'Position'); % 403 246 560 420
set(h,'Position',[Poscur(1) Poscur(2) 2*Poscur(3) Poscur(4)]);
set(ha,'FontSize',FontSize);
leg1 = 'Rx';
legend(leg1,'Location','SouthEast')
Saveas(h,[dirfigures 'Rx'],'emf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[insig fs] = Wavread(file);

insig = From_dB(10)*insig;

B = il_rx_piano(fc,Rx,fs); 
Bm20 = il_rx_piano(fc,Rx-20,fs); 
N = length(B);

outsigm20 = filter(Bm20,1,insig);
outsig    = filter(B   ,1,insig);

filein      = [diraudio 'T8-sample-Klavier.wav'];
fileoutReal = [diraudio 'T8-sample-Klavier-out-Att-real.wav'];
Wavwrite( insig,fs,filein);
Wavwrite(outsig   ,fs,fileoutReal);
Wavwrite(outsigm20,fs,[diraudio 'T8-sample-Klavier-out-Att-Rxm20.wav']);
%%%
option = [];
option.nAnalyser = 10;
option.CalMethod = 1;
option.bPlot = 0;
out = PsySoundCL(filein,option);
OB_1 = out.DataSpecOneThirdAvg(5:22);
plot(1:length(fc),OB_1,'rs-','LineWidth',2);
leg2 = 'piano/tenant';
legend({leg1,leg2},'Location','SouthEast')
Saveas(h,[dirfigures 'Rx+piano'],'emf');
%%%

option = [];
option.nAnalyser = 10;
option.CalMethod = 1;
option.bPlot = 0;
out = PsySoundCL(fileoutReal,option);
OB_2 = out.DataSpecOneThirdAvg(5:22);
plot(1:length(fc),OB_2,'k>-.','LineWidth',2);
leg3 = 'piano/neighbour';
legend({leg1,leg2,leg3},'Location','SouthEast')
ylim([20 80])
Saveas(h,[dirfigures 'Rx+piano+Att'],'emf');

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = il_rx_piano(f,Rx,fs)

f = [0 f fs/2];
Rx = -1*[Rx(1) Rx Rx(end)];

N = 2^12;

B = fir2(2^12,f/(fs/2),From_dB(Rx));
[H1,Fn]=freqz(B,1);

% figure;
% semilogx(Fn, 20*log10(abs([H1])));
% 
% xlabel('Frequency [Hz]')
% legend([num2str(N) ' taps']);
% title('FIR filter to be used as approximation to isolation curve')

disp('')