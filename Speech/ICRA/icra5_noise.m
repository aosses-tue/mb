function r=icra5_noise(data,fs,method)
% function r=icra5_noise(data,fs,method)
%
% 1. Description:
%       Imitate the making of the icra5 noise (Dreschler et al, 2005, Int J Aud)
%
% 2. Stand-alone example:
%       file = [Get_TUe_data_paths('db_speechmaterials') 'Spanish' delim 'Matrix' delim '00131.wav'];
%       [insig fs] = Wavread(file);
%       outsig = icra5_noise(insig);
%       lvl = rmsdb(insig);
%       outsig = setdbspl(outsig,lvl+100);
%       sound(insig,fs);
%       pause(length(insig)/fs*1.2);
%       sound(outsig,fs);
% 
%       file = [Get_TUe_data_paths('piano') '01-Chabassier' delim 'SONS' delim 'Cd5' delim 'pressionexpe.wav'];
%       [insig fs] = Wavread(file);
%       outsig = icra5_noise(insig);
%       lvl = rmsdb(insig);
%       outsig = setdbspl(outsig,lvl+100);
%       sound(insig,fs);
%       pause(length(insig)/fs*1.2);
%       sound(outsig,fs);
% 
% 3. See also:
%       r20151211_update_ICRA
% 
% Author: Tom Francart 2010-2013, KU Leuven
% Modified by: Alejandro Osses, HTI, TU/e 2015
% Last used on: 07/12/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    method = 1;
end

if fs < 2
    fs=44100;
end

RMSin   = rmsdb(data) + 100;

[B1,B2,B3,A1,A2,A3]=il_getfilters(fs);

r       = zeros(length(data),3);
cross   = [0 800 2400 fs/2]; % same as used in il_getfilters

BW(1)   = cross(2)-cross(1);
BW(2)   = cross(3)-cross(2);
BW(3)   = cross(4)-cross(3);

r(:,1)=filter(B1,A1,data);
r(:,2)=filter(B2,A2,data);
r(:,3)=filter(B3,A3,data);

r(:,1)=il_schroeder(r(:,1));
r(:,2)=il_schroeder(r(:,2));
r(:,3)=il_schroeder(r(:,3));

interimRMS          = rmsdb(r)+100;
interimRMS_per_Hz   = interimRMS - 10*log10(BW); 

bDebug = 1;
if bDebug == 1
    K = 8192/2;
    figure;
    subplot(2,1,1)
    freqfft2(r,K,fs);
    legend('LP','BP','HP')
end
% rms=meanrms(r);

r(:,1)=filter(B1,A1,r(:,1));
r(:,2)=filter(B2,A2,r(:,2));
r(:,3)=filter(B3,A3,r(:,3));

for i=1:3
    switch method
        case 1
            r(:,i)=r(:,i)/meanrms(r(:,i))*sqrt(BW(i));
        case 2
            r(:,i) = setdbspl( r(:,i),interimRMS(i) );
    end
end

finalRMS          = rmsdb(r)+100;
finalRMS_per_Hz   = finalRMS - 10*log10(BW); 

if bDebug == 1
    subplot(2,1,2)
    freqfft2(r,K,fs);
    legend('LP','BP','HP')
end

r=sum(r,2);
if method == 1
   r = setdbspl(r,RMSin);
end

if bDebug == 1
    figure;
    freqfft2(r,K,fs);
    title('Signals added together')
end

RMSbefore = rmsdb(r)+100;
r       = il_randomize_phase(r); % it can increase or decrease the level
B       = malespectrum_filter(fs); % effort filter not being applied inside (only for plot)

N2pad   = round(length(B)/2);
r       = [r; zeros(N2pad,1)]; 

r       = filter(B,1,r);
r       = r(N2pad+1:end);
RMSafter= rmsdb(r)+100;
r       = gaindb(r,RMSbefore-RMSafter);

disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B1,B2,B3,A1,A2,A3]=il_getfilters(fs)

cross=[800 2400];

% order=100/20;          % butterworth: slope=20*order
order=floor(84/6);       % butterworth: A = 6 dB/Oct * Order

[B1,A1]=butter(order  , cross(1)/fs*2);
[B2,A2]=butter(order/2, cross/fs*2);
[B3,A3]=butter(order  , cross(2)/fs*2,'high');

doplot=1;
if (doplot)
    [H1,F]=freqz(B1,A1,[],fs);
    [H2,F]=freqz(B2,A2,[],fs);
    [H3,F]=freqz(B3,A3,[],fs);
    figure
    semilogx(F,20*log10([abs(H1) abs(H2) abs(H3)]));
end

function r=il_schroeder(d)
r=d.*sign(rand(length(d),1)-0.5);

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
