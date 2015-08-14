function y = r20150814_update(x)
% function y = r20150814_update(x)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 12/08/2015
% Last update on: 12/08/2015 
% Last use on   : 12/08/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bDiary = 0;
Diary(mfilename,bDiary);


bDoEnveAMT = 0;
bDoHartmann = 0;
bDoTFMF = 1;
bDoTFDRNL = 1;

dire = [Get_TUe_paths('outputs') 'AMTControl-examples' delim]; % 'D:\Output\AMTControl-examples\'
dBFS = 100;

if bDoEnveAMT
    
    f1 = [dire 'randomnoise-BBN_SPL-70.wav'];
    [insig1 fs] = Wavread(f1);

    fc = 1100;
    BW = 200;
    fcut = 1000;

    SpectrumLvl = 40;
    lvl1 = SpectrumLvl + 10*log10(fs/2);
    lvl2 = SpectrumLvl + 10*log10(BW);

    windowtype = 'hanning';
    K  = length(insig1)/2;
    [xx y1dB f] = freqfft2(insig1,K,fs,windowtype,dBFS);

    insig1 = setdbspl(insig1,lvl1,dBFS);

    insig2 = Set_Fourier_coeff_to_zero(insig1,fs,fc-BW/2,fc+BW/2);
    insig2 = setdbspl(insig2,lvl2,dBFS);
    [xx y2dB f] = freqfft2(insig2,K,fs,windowtype,dBFS);

    figure; 
    plot(f,y1dB,'b',f,y2dB,'r'); grid on
    ylim([0 70])

    [b_highest,a_highest] = butter(2,fcut/(fs/2));

    [xx yenvdB1] = envfreqfft(insig1,K,fs,'hanning',dBFS,fcut);

    % Ntimes = 100;
    % 
    % ytmp = [];
    % for i = 1:Ntimes
    %     % [xx y4dB f] = freqfft2(insig4,K,fs,windowtype,dBFS);
    %     outsigtmp = Randomise_insig(insig1);
    %     yenv    = abs(hilbert(outsigtmp));
    %     insig4 = filter(b_highest, a_highest,yenv);
    %     [xx yy f] = freqfft2(insig4,K,fs,windowtype,dBFS);
    %     ytmp = [ytmp yy];
    % end
    % yenvdB1 = transpose( dbmean(transpose(ytmp)) );

    % ytmp = [];
    % for i = 1:Ntimes
    %     % [xx y4dB f] = freqfft2(insig4,K,fs,windowtype,dBFS);
    %     outsigtmp = Randomise_insig(insig2);
    %     yenv    = abs(hilbert(outsigtmp));
    %     insig5 = filter(b_highest, a_highest,yenv);
    %     [xx yy f] = freqfft2(insig5,K,fs,windowtype,dBFS);
    %     ytmp = [ytmp yy];
    % end
    % yenvdB2 = transpose( dbmean(transpose(ytmp)) );
    [xx yenvdB2] = envfreqfft(insig2,K,fs,'hanning',dBFS,fcut);

    figure; 
    plot(f,yenvdB1,'b',f,yenvdB2,'r'); grid on
    ylim([0 70])
    xlabel('Envelope frequency [Hz]')
    ylabel('Amplitude')

end

if bDoHartmann
    
    % Envelope, Ex 1 (Hartmann2005, pp 418):
    f = 980:10:1020;
    A = [1/4 1/2 1 1/2 1/4];
    phi = -pi/2;
    dur = 120e-3;
    fs = 44100;
    
    y = zeros(dur*fs,1);
    for i = 1:5
        ytmp = A(i) * Create_sin_phase(f(i),phi,dur,fs);
        y = y+ytmp;
    end
    
    t = (1:length(y)) /fs;
    env = abs( 1 + 0.5*cos(40*pi*t) + cos(20*pi*t) );
    
    figure;
    plot(t*1000, env,'r',t*1000,y,'b'); hold on
    plot(t*1000,-env,'r')
    xlabel('Time [ms]')
    ylabel('Amplitude')
    grid on
    legend('Envelope')
    
    % Useless envelope, Ex 3:
    f = [100 1000];
    A = [1 1/4];
    phi = 0;
    dur = 20e-3; % 2 100-Hz periods
    fs = 44100;
    y = zeros(dur*fs,1);
    for i = 1:2
        ytmp = A(i) * Create_sin_phase(f(i),phi,dur,fs);
        y = y+ytmp;
    end
    
    t = (1:length(y)) /fs;
    env = sqrt( 17/16 + 0.5*cos(2*pi*900*t) );
    
    figure;
    plot(t*1000, env,'r',t*1000,y,'b'); hold on
    plot(t*1000,-env,'r')
    xlabel('Time [ms]')
    ylabel('Amplitude')
    grid on
    legend('Envelope')
    
    SPL = 70;
    insig = setdbspl(y,SPL,dBFS);
    filename = sprintf('%ssine-%.0f-plus-%.0f-Hz-%.0f-dB.wav',dire,f(1),f(2),SPL);
    Wavwrite(insig,fs,filename);
    
end

fmax2plot = 1000;
N = 8192*8;

if bDoTFMF
    
    % To do: compare outsig with TFs
    insig = [zeros(N/2-1,1); 1; zeros(N/2,1)];
    fs = 44100;
    fc = 4000;
    [outsig,mfc,params] = modfilterbank_debug(insig,fs,fc);
    outsig = outsig{1};
    
    K = N/2;
    [x xdB f] = freqfft2(insig,K,fs);
    
    h(:,1) = freqz(params.b_mf1, params.a_mf1, K);
        
    ydB = [];
    for i = 2:length(mfc)
        
        exp1 = sprintf('h(:,%.0f) = freqz(params.b_mf%.0f, params.a_mf%.0f, K);',i,i,i);
        eval(exp1);
        
    end
    figure;
    plot(f,To_dB(abs(h))); grid on;
    xlim([0 fmax2plot])
    ylim([-18 5])
    xlabel('Modulation frequency [Hz]')
    ylabel('Attenuation [dB]')
    title('Only modulation filterbank')
    
    figure;
    hhigh = freqz(params.b_highest,params.a_highest,K);
    hhigh = repmat(hhigh,1,length(mfc));
    plot(f,To_dB(abs(h.*hhigh))); grid on; hold on
    xlim([0 fmax2plot])
    ylim([-18 5])
    xlabel('Modulation frequency [Hz]')
    ylabel('Attenuation [dB]')
    plot(f,To_dB(hhigh(:,1)),'k--','LineWidth',2)
     
end

if bDoTFDRNL
    
    
    
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
