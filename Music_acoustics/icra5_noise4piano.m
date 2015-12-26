function r=icra5_noise4piano(insig,fs,method)
% function r=icra5_noise4piano(insig,fs,method)
%
% 1. Description:
%       Imitate the making of the icra5 noise (Dreschler et al, 2005, Int J Aud)
%
% 2. Stand-alone example:
%       file = [Get_TUe_data_paths('db_speechmaterials') 'Spanish' delim 'Matrix' delim '00131.wav'];
%       [insig fs] = Wavread(file);
%       outsig = icra5_noise4piano(insig);
%       lvl = rmsdb(insig);
%       outsig = setdbspl(outsig,lvl+100);
%       sound(insig,fs);
%       pause(length(insig)/fs*1.2);
%       sound(outsig,fs);
% 
%       file = [Get_TUe_data_paths('piano') '01-Chabassier' delim 'SONS' delim 'Cd5' delim 'pressionexpe.wav'];
%       [insig fs] = Wavread(file);
%       outsig = icra5_noise4piano(insig);
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
% First modified on: 09/12/2015
% Last used on: 23/12/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    method = 3; % Method 3 = with ERB bands, 'approved' (nice first attempt) by AK on 22/12/2015
end

if nargin < 2
    fs=44100;
end

switch method
    case 2
        type    = 1;
        cross   = Get_icra_cutoff_freqs(type);
        [B,A]   = il_getfilters(fs,cross);
        Nbands  = length(cross)+1; 
        r       = zeros(length(insig),Nbands);
        for i = 1:Nbands
            r(:,i)  = filter(B(i,:),A(i,:),insig);
        end
    case 3
        [r, fc] = auditoryfilterbank(insig,fs);
        Nbands = length(fc);
end

interimRMS = rmsdb(r)+100;

for i = 1:Nbands
    r(:,i)  = il_schroeder(r(:,i));
end

switch method
    case 2
        for i = 1:Nbands
            r(:,i)  = filter(B(i,:),A(i,:),r(:,i));
        end
    case 3
        for i = 1:Nbands
            rtmp = auditoryfilterbank_one_freq(r(:,i),fs,fc(i));
            r(:,i) = rtmp;
        end
        
end

interimGlobalRMS    = rmsdb(sum(r,2))+100;
interimRMS2         = rmsdb(r)+100;

bDebug = 0;
if bDebug == 1
    K = 8192/2;
    figure;
    subplot(2,1,1)
    freqfft2(r,K,fs);
    legend('LP','BP','HP')
end

% currentRMS = rmsdb(r)+100; % figure; plot(interimRMS,'ro'); hold on; plot(currentRMS,'s--')

for i=1:Nbands
    r(:,i) = setdbspl( r(:,i),interimRMS(i) );
end

interimRMS3 = rmsdb(r)+100;
r=sum(r,2);

RMSbefore = rmsdb(r)+100;
r         = il_randomise_phase(r); % it can increase or decrease the level

% B       = malespectrum_filter(fs); % effort filter not being applied inside (only for plot)
% N2pad   = round(length(B)/2);
% r       = [r; zeros(N2pad,1)]; 
% r       = filter(B,1,r);
% r       = r(N2pad+1:end);
RMSafter  = rmsdb(r)+100;
r         = gaindb(r,RMSbefore-RMSafter); % compensate decrease or increase in level after phase randomisation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B,A]=il_getfilters(fs,cross)

if fs <= 22050          % tested for fs = 20000 Hz
   order=floor(84/6);   % butterworth: A = 6 dB/Oct * Order
else                    % tested for fs = 48000 Hz
   order=floor(60/6);   % at 72 starts a little bit of ripple
end

Nbands = length(cross)+1;

B = zeros(Nbands,order+1);
A = zeros(Nbands,order+1);

[B(1,:),A(1,:)]=butter(order  , cross(1)/fs*2);
for i = 2:Nbands-1
    [B(i,:),A(i,:)]=butter(order/2, cross(i-1:i)/fs*2);
end
[B(Nbands,:),A(Nbands,:)]=butter(order, cross(Nbands-1)/fs*2,'high');

doplot=1;
if (doplot)
    
    for i = 1:Nbands
        [H(:,i),F]=freqz(B(i,:),A(i,:),[],fs);
    end
    figure
    semilogx(F,20*log10([abs(H)]));
end

function r=il_schroeder(d)
r=d.*sign(rand(length(d),1)-0.5);

function r=il_randomise_phase(d)
N=512;      % window length
hop=512/8;
window=hamming(N);
r=zeros(length(d),1);
for i=0:round(length(d)/hop)-1-N/hop
   from=i*hop+1;
   to=i*hop+N;
   stuk=fft(d(from:to).*window,N);
   stuk=stuk.*exp(j*rand(N,1)*2*pi);
   r(from:to)=r(from:to)+real(ifft(stuk));
end
