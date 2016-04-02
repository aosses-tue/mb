function [r info4pede pede_sample] = icra_noise4piano(insig,fs)
% function [r info4pede pede_sample] = icra5_noise4piano(insig,fs)
%
% 1. Description:
%       Imitate the making of the icra5 noise (Dreschler et al, 2005, Int J Aud)
%       The processing introduced here is similar to that of icra5_noise4piano
%       using method = 3
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
%       r20160316_update_Antoine.m, icra5_noise4piano.m
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 21/03/2016
% Last update on: 21/03/2016 
% Last use on   : 21/03/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    fs=44100;
end

[r, fc] = auditoryfilterbank(insig,fs);
Nbands = length(fc);

interimRMStot = rmsdb(insig)+100;
interimRMS = rmsdb(r)+100;

for i = 1:Nbands
    r(:,i)  = il_schroeder(r(:,i));
end

for i = 1:Nbands
    rtmp = auditoryfilterbank_one_freq(r(:,i),fs,fc(i));
    r(:,i) = rtmp;
end   

% for i=1:Nbands
%     r(:,i) = setdbspl( r(:,i),interimRMS(i) );
% end

r         = sum(r,2);
r         = il_randomise_phase(r); % it can increase or decrease the level

ri = Apply_IIR_Butter(r ,fs,audtofreq( 2.5),'high',4);
ro = Apply_IIR_Butter(ri,fs,audtofreq(33.5),'low',8);

r         = setdbspl( r,interimRMStot );
% r = gaindb(r,RMSbefore-RMSafter); % compensate decrease or increase in level after phase randomisation
info4pede.RMS    = interimRMS - max(interimRMS);
info4pede.RMStot = rmsdb(r);
info4pede.RMSmax = rmsdb(prctile(abs(r),99));
info4pede.fc     = fc;
info4pede.fs     = fs;
info4pede.Length = length(r);

if nargout >= 3
    
    SNR_dB = 0;
    pede_sample = icra_noise4piano_pedestal(r,fs,info4pede.RMS,SNR_dB);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r=il_schroeder(d)
r=d.*sign(rand(length(d),1)-0.5);

function r = il_randomise_phase(d)
N      = 512; % window length
hop    = 512/8;
window = hamming(N);
r      = zeros(length(d),1);
for i=0:round(length(d)/hop)-1-N/hop
   from = i*hop+1;
   to   = i*hop+N;
   stuk = fft(d(from:to).*window,N);
   stuk = stuk.*exp(j*2*pi*rand(N,1));
   r(from:to)=r(from:to)+real(ifft(stuk));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
