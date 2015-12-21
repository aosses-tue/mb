function r=icra5_noise(insig,fs,method,effort_type,opts)
% function r=icra5_noise(insig,fs,method,effort_type,opts)
%
% 1. Description:
%       Imitate the making of the icra5 noise (Dreschler et al, 2005, Int J Aud)
%
% 2. Stand-alone example:
%       file = [Get_TUe_data_paths('db_speechmaterials') 'Spanish' delim 'Matrix' delim '00131.wav'];
%       [insig fs] = Wavread(file);
%       method = 1;
%       effort_type = 'raised';
%       outsig = icra5_noise(insig,fs,method,effort_type);
%       lvl = rmsdb(insig);
%       outsig = setdbspl(outsig,lvl+100);
%       sound(insig,fs);
%       pause(length(insig)/fs*1.2);
%       sound(outsig,fs);
% 
%       file = [Get_TUe_data_paths('piano') '01-Chabassier' delim 'SONS' delim 'Cd5' delim 'pressionexpe.wav'];
%       [insig fs] = Wavread(file);
%       outsig = icra5_noise(insig,fs);
%       lvl = rmsdb(insig);
%       outsig = setdbspl(outsig,lvl+100);
%       sound(insig,fs);
%       pause(length(insig)/fs*1.2);
%       sound(outsig,fs);
% 
% 3. See also:
%       r20151211_update_ICRA;
% 
% Author: Tom Francart 2010-2013, KU Leuven
% Modified by: Alejandro Osses, HTI, TU/e 2015
% Last used on: 07/12/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDebug = 0;

if nargin < 5
    opts = [];
end

if nargin < 4
    effort_type = 'normal'; % 'raised', 'loud','shout'
end

if nargin < 3
    method = 1;
end

if nargin < 2
    fs=44100;
end

opts = ef(opts,'bPlot',0);
opts = Ensure_field(opts,'Gender','male');
bPlot  = opts.bPlot;
Gender = opts.Gender;

%%%
if bDebug
    calcMTF(insig,fs);
end
%%%

RMSin   = rmsdb(insig) + 100;

[B1,B2,B3,A1,A2,A3]=il_getfilters(fs,bPlot);

r       = zeros(length(insig),3);
cross   = [0 800 2400 fs/2]; % same as used in il_getfilters

BW(1)   = cross(2)-cross(1);
BW(2)   = cross(3)-cross(2);
BW(3)   = cross(4)-cross(3);

r(:,1)=filter(B1,A1,insig);
r(:,2)=filter(B2,A2,insig);
r(:,3)=filter(B3,A3,insig);

%%%
r(:,1)=il_schroeder(r(:,1));
r(:,2)=il_schroeder(r(:,2));
r(:,3)=il_schroeder(r(:,3));

%%%
if bDebug
    calcMTF(sum(r,2),fs);
end
%%%

% interimGlobalRMS = rmsdb(sum(r,2))+100;
interimRMS          = rmsdb(r)+100;
% interimRMS_per_Hz   = interimRMS - 10*log10(BW); 

if bPlot == 1
    K = 8192/2;
    figure;
    subplot(2,1,1)
    freqfft2(r,K,fs);
    legend('LP','BP','HP')
end

r(:,1)=filter(B1,A1,r(:,1));
r(:,2)=filter(B2,A2,r(:,2));
r(:,3)=filter(B3,A3,r(:,3));

%%%
% calcMTF(sum(r,2),fs);
%%%

ra = r;
for i=1:3
    switch method
        case 1
            r(:,i)=r(:,i)/rms(r(:,i))*sqrt(BW(i)); % 1/rms(r(:,i)) = makes RMS equal to 1
        case 2
            r(:,i) = setdbspl( r(:,i),interimRMS(i) );
    end
end

%%%
if bDebug
    calcMTF(sum(r,2),fs);
end
%%%

finalRMS          = rmsdb(r)+100;
finalRMS_per_Hz   = finalRMS - 10*log10(BW); 

if bPlot == 1
    subplot(2,1,2)
    freqfft2(r,K,fs);
    legend('LP','BP','HP')
end

r=sum(r,2);
if method == 1
   r = setdbspl(r,RMSin);
end

%%%
if bDebug
    calcMTF(sum(r,2),fs);
    title('5')
end
%%%

if bPlot == 1
    figure;
    freqfft2(r,K,fs);
    title('Signals added together')
end

RMSbefore = rmsdb(r)+100;
%%% Phase randomising:
r       = il_randomize_phase(r); % it can increase or decrease the level

%%%
if bDebug
    calcMTF(sum(r,2),fs);
    title('6')
end
%%%

%%% Gender filter:
switch Gender
    case 'male'
        B = malespectrum_filter(fs,'ansi',bPlot); % effort filter not being applied inside (only for plot)
    case 'female'
        B = femalespectrum_filter(fs,'ansi',bPlot);
end

N2pad   = round(length(B)/2);
r       = [r; zeros(N2pad,1)]; 

r       = filter(B,1,r);
r       = r(N2pad+1:end);
RMSafter= rmsdb(r)+100;
r       = gaindb(r,RMSbefore-RMSafter); % compensate decrease or increase in level after phase randomisation

%%% Effort filter
if ~strcmp(effort_type,'normal');
    B   = il_Get_effort_filter(fs,effort_type,bPlot);
    N2pad = round(length(B)/2);
    r   = [r; zeros(N2pad,1)]; 
    r   = filter(B,1,r);
    r   = r(N2pad+1:end);
end

%%%
if bDebug
    calcMTF(sum(r,2),fs);
    title('7')
end
%%%

disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B1,B2,B3,A1,A2,A3]=il_getfilters(fs,bPlot)

cross=[800 2400];

% order=100/20;         % butterworth: slope=20*order
if fs <= 22050          % tested for fs = 20000 Hz
   order=floor(84/6);   % butterworth: A = 6 dB/Oct * Order
else                    % tested for fs = 48000 Hz
   order=floor(66/6);   % at 72 starts a little bit of ripple
end

[B1,A1]=butter(order         , cross(1)/fs*2);
[B2,A2]=butter(round(order/2), cross/fs*2);
[B3,A3]=butter(order         , cross(2)/fs*2,'high');

if (bPlot)
    [H1,F]=freqz(B1,A1,[],fs);
    [H2,F]=freqz(B2,A2,[],fs);
    [H3,F]=freqz(B3,A3,[],fs);
    figure
    semilogx(F,20*log10([abs(H1) abs(H2) abs(H3)]));
    title(sprintf('fs=%.0f [Hz]',fs))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=il_schroeder(d)
r=d.*sign(rand(length(d),1)-0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=il_randomize_phase(d)
N=512;      % window length
hop=512/8;
window=hamming(N);
r=zeros(length(d),1);
for i=0:round(length(d)/hop)-1-N/hop
   from = i*hop+1;
   to   = i*hop+N;
   stuk = fft(d(from:to).*window,N);
   stuk = stuk.*exp(j*rand(N,1)*2*pi);
   r(from:to)=r(from:to)+real(ifft(stuk));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = il_Get_effort_filter(fs,effort_type,bPlot)

cross = [100 500];
abs_slope = [12 22];

[E  f] = SpeechSptr(effort_type);
[En f] = SpeechSptr('normal');

% F=logspace(0, log10(fs/2));
F = [0 f fs/2];

% R=zeros(1,length(F)); % response, memory allocation

ref=log2(cross(1));
i = 1;
while (F(i)<cross(1)) & F(i)<f(1)
    Ediff_dB(i)=10^((ref-log2(F(i)))*(-abs_slope(1))/20);
    i=i+1;
end
Ediff_dB= [Ediff_dB E-En];

ref=log2(cross(2));
Ediff_dB(end+1)=10^((log2(F(end))-ref)*(-abs_slope(2))/20);

Ediff   = From_dB(Ediff_dB)
B       = fir2(10000,F/fs*2, Ediff);

if (bPlot)
    
    figure
    semilogx(F,20*log10(Ediff),'-xr');
    title( sprintf('Effort filter: %s',effort_type) );
    grid on;
    hold on
    [resp,Fresp]=freqz(B,1,[],fs);
    plot(Fresp,20*log10(abs(resp)), '-ob');
end

